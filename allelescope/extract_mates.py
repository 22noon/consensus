import pysam
import sys
from collections import defaultdict


#######  Indel tagger code ############

def summarise_indels(indels):
    if not indels:
        return 0, '', 0, 0
    parts   = []
    max_del = max_ins = net = 0
    for typ, length in indels:
        parts.append(f"{typ}{length}")
        if typ == 'D':
            max_del  = max(max_del, length)
            net     -= length
        else:
            max_ins  = max(max_ins, length)
            net     += length
    return net, ','.join(parts), max_del, max_ins


def get_allele_at_base(read, target_pos):
    """
    Returns a descriptive allele string for what the read shows
    at target_pos (0-based reference coordinate).

    Returns one of:
      'D<len>'    — deletion e.g. 'D7'
      'I<seq>'    — insertion e.g. 'IATT' (actual inserted bases)
      'R<base>'   — reference match e.g. 'RA'
      'M<base>'   — mismatch/SNP e.g. 'MT'
      'S'         — soft clipped at this position
      '.'         — no coverage / unmapped
    """
    if read.is_unmapped or read.cigartuples is None:
        return '.'

    query_seq = read.query_sequence
    query_pos = 0       # index into query_sequence (includes soft clip)
    ref_pos   = read.reference_start

    for op, length in read.cigartuples:

        # ── Soft clip ────────────────────────────────────────────────
        if op == 4:
            # Soft clip consumes query but not reference.
            # Check if target_pos falls under the clipped region.
            # For leading clip: ref_pos hasn't moved yet so compare
            # to read start; for trailing: ref_pos is past target.
            for i in range(length):
                # We can only say "this base is soft-clipped"
                # Map soft clip bases to ref approximately:
                # leading clip sits just before reference_start,
                # so we don't assign them reference positions —
                # just advance query and move on.
                query_pos += 1
            continue

        # ── Hard clip / padding ──────────────────────────────────────
        if op in (5, 6):
            continue

        # ── Insertion ───────────────────────────────────────────────
        if op == 1:
            # Insertion is anchored between ref_pos-1 and ref_pos.
            # Convention: insertion belongs to the ref base just before it.
            if ref_pos - 1 == target_pos or ref_pos == target_pos:
                ins_seq = query_seq[query_pos: query_pos + length]
                return f"I{ins_seq}"
            query_pos += length
            continue

        # ── Deletion / intron skip ───────────────────────────────────
        if op in (2, 3):
            if ref_pos <= target_pos < ref_pos + length:
                return f"D{length}"
            ref_pos += length
            continue

        # ── Alignment match / mismatch / sequence match ──────────────
        # op 0=M, 7== (EQ), 8=X
        if op in (0, 7, 8):
            if ref_pos <= target_pos < ref_pos + length:
                offset = target_pos - ref_pos
                base   = query_seq[query_pos + offset].upper()
                # For EQ (7) we know it's a match; for X (8) mismatch;
                # for M (0) we need to compare to reference — but we
                # don't have the reference here, so use base quality
                # or the MD tag to discriminate. See note below.
                if op == 7:
                    return f"R{base}"
                elif op == 8:
                    return f"M{base}"
                else:
                    ref_base = get_ref_base_from_md(read, target_pos)
                    if ref_base is None:
                        return f"R{base}"   # can't tell, assume ref
                    elif base == ref_base:
                        return f"R{base}"
                    else:
                        return f"M{base}"
            ref_pos   += length
            query_pos += length
            continue

    return '.'   # target_pos not covered by this read


def get_ref_base_from_md(read, target_pos):
    """
    Parse the MD tag to find the reference base at target_pos.
    MD looks like: '10A5^TT3' meaning:
      10 matches, then ref=A (mismatch), 5 matches, then deletion of TT, 3 matches.
    Returns the uppercase reference base, or None if MD tag absent/unparseable.
    """
    try:
        read.get_tag("MD")
    except KeyError:
        return None


    # Walk CIGAR in parallel with MD to track reference base at target
    # MD only covers aligned (non-clipped) bases, so we need CIGAR too
    # to correctly advance query_pos through soft clips first.

    # Precompute: offset into the MD-covered positions
    # Simpler approach: use pysam's get_aligned_pairs with seq
    try:
        for qpos, rpos, ref_base in read.get_aligned_pairs(with_seq=True, matches_only=False):
            if rpos is None:
                continue    # insertion in query — no ref pos
            if rpos == target_pos:
                if ref_base is None:
                    return None
                return ref_base.upper()
    except (ValueError, TypeError):
        return None

    return None


def get_indels_at_base(read, target_pos):
    """
    Returns list of (type, length) for all indels overlapping target_pos.
    Kept for the XI/XD/XN tags.
    """
    if read.cigartuples is None:
        return []
    indels  = []
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op == 1:
            if ref_pos == target_pos or ref_pos - 1 == target_pos:
                indels.append(('I', length))
        elif op == 2:
            if ref_pos <= target_pos < ref_pos + length:
                indels.append(('D', length))
            ref_pos += length
        elif op == 3:
            ref_pos += length
        elif op in (0, 7, 8):
            ref_pos += length
    return indels


###########  End of indel tagger code, start of mate extraction code ############

Tag_bits = {
    'OVERLAP':          0x1,
    'NON_SPANNING_MATE': 0x2,
    'SPLITX':           0x4,
    'SPLIT':            0x8,
    'PROPER':           0x10,
    'IMPROPER':         0x20,
    'READ1':            0x40,
    'READ2':            0x80,
    'SECONDARY':        0x100,
    'QCFAIL':           0x200,
    'DUPLICATE':        0x400,
    'SUPPLEMENTARY':    0x800,
}


def parse_sa_tag(sa_string):
    """Parse SA tag to get locations of other alignments."""
    parts = []
    for entry in sa_string.rstrip(';').split(';'):
        rname, pos, strand, cigar, mapq, nm = entry.split(',')
        parts.append({
            'rname':  rname,
            'pos':    int(pos),
            'strand': strand,
            'cigar':  cigar,
            'mapq':   int(mapq),
            'nm':     int(nm),
        })
    return parts


def calculate_softclip_lengths(read):
    cigar_tuples = read.cigartuples
    left_sc, right_sc = 0, 0
    if cigar_tuples:
        left_sc  = cigar_tuples[0][1]  if cigar_tuples[0][0]  == 4 else 0
        right_sc = cigar_tuples[-1][1] if cigar_tuples[-1][0] == 4 else 0
    read.set_tag("SL", left_sc,  value_type="i")
    read.set_tag("SR", right_sc, value_type="i")
    return read


def get_pair_orientation(read):
    """Determine the orientation of a paired read (F1R2, R2F1, etc.)."""
    if read.is_unmapped or read.mate_is_unmapped:
        return "UNMAPPED"
    is_read1    = read.is_read1
    read_reverse = read.is_reverse
    mate_reverse = read.mate_is_reverse
    if is_read1:
        read1_strand = 'R1' if read_reverse else 'F1'
        read2_strand = 'R2' if mate_reverse else 'F2'
        return f"{read1_strand}{read2_strand}"
    else:
        read2_strand = 'R2' if read_reverse else 'F2'
        read1_strand = 'R1' if mate_reverse else 'F1'
        return f"{read1_strand}{read2_strand}"


def tag_read(read, XO):
    if read.has_tag('SA'):
        XO |= Tag_bits['SPLIT']
        #return read
    orientation = get_pair_orientation(read)
    if read.is_paired:
        read.set_tag('XR', orientation)
        if orientation != "F1R2" or not read.is_proper_pair:
            XO |= Tag_bits['IMPROPER']
        else:
            XO |= Tag_bits['PROPER']
    read.set_tag('XO', XO)
    return read

"""

Extract every alignment (split or not) overlapping a genomic region of
interest into ``out_bam``, and the mates / split-alignment partners of
those overlapping reads into ``mates_bam``.

Assumes the following are already available in scope (imported from the
existing project modules, same as the original implementation):

    Tag_bits                -> dict of RG-flag bit values
    tag_read(read, XO)      -> sets the XO tag on ``read`` and returns it
    calculate_softclip_lengths(read)
    parse_sa_tag(sa_string)
    get_allele_at_base(read, target_pos)
    get_indels_at_base(read, target_pos)
    summarise_indels(indels)
"""

def extract_mates(
    src_bam: str,
    out_bam: str,
    mates_bam: str,
    chrom: str,
    start: int,
    end: int,
) -> defaultdict:
    """
    Extract overlapping reads and their mates for a genomic region.

    OVERLAP semantics
    ------------------
    A given alignment record is tagged OVERLAP only if BOTH of the
    following hold:
      1. It physically crosses the region of interest itself 
      2. Both mates of the fragment (R1 and R2) have an alignment
         that crosses the region of interest. This is meant to capture
         "independent sequencing evidence from both strands of the
         insert" -- a lone read (or only one mate) crossing the region
         is not enough

    Any alignment (of either mate, split or not) that does not itself
    cross the region but belongs to a fragment where at least one other
    alignment does, is fetched separately and written to ``mates_bam``,
    tagged NON_SPANNING_MATE (paired-end mate) or SPLITX (split-read
    partner) as appropriate.

    Args:
        src_bam:   Input BAM path (e.g. alignment.markdup.bam)
        out_bam:   Output BAM for overlapping reads at the locus
        mates_bam: Output BAM for fetched mates / split partners
        chrom:     Chromosome name
        start:     1-based start position
        end:       1-based end position

    Returns:
        allele_counts: dict mapping allele string -> read count (one
        observation per physical read/fragment, deterministic).
    """
    allele_counts = defaultdict(int)

    bam = pysam.AlignmentFile(src_bam, "rb")
    header = bam.header.to_dict()

    if 'RG' not in header:
        header['RG'] = []
    rg_ids = [
        'NON_SPANNING_MATE',
        'F1R2', 'R2F1', 'F2R1', 'R1F2', 'F1F2', 'R1R2',
        'SPLIT', 'SPLITX', 'OVERLAP', 'IMPROPER', 'UNPAIRED', 'PROPER',
    ]
    existing_ids = {rg.get('ID') for rg in header['RG']}
    for rg_id in rg_ids:
        if rg_id not in existing_ids:
            header['RG'].append({'ID': rg_id, 'SM': 'sample', 'PL': 'ILLUMINA'})
            existing_ids.add(rg_id)

    overlap_bam = pysam.AlignmentFile(out_bam,   "wb", header=header)
    mate        = pysam.AlignmentFile(mates_bam, "wb", header=header)

    # ``start``/``end`` are part of the public, 1-based API. The tag
    # helpers and pysam both operate on 0-based, half-open coordinates.
    region_start0 = start - 1
    region_end0   = end
    target_pos    = start - 1

    def overlaps_region(read) -> bool:
        return read.reference_start < region_end0 and read.reference_end > region_start0

    written_records = set()               # (query_name, reference_start) already in overlap_bam
    XA_dict          = {}                 # query_name -> chosen allele call, propagated onto mates
    pending          = defaultdict(list)  # query_name -> [(chrom, pos, kind), ...]
    covering         = defaultdict(list)  # query_name -> [(is_supplementary, reference_start, allele), ...]

    # ---- Pass 1a: which mates (R1/R2) have a PRIMARY alignment crossing the region ----
    print(f"Scanning region {chrom}:{start}-{end} for overlapping alignments >(split or not)<...")
    buffered_reads = list(bam.fetch(chrom, region_start0, region_end0))

    overlap_mate = defaultdict(set)  # query_name -> subset of {1, 2}
    for read in buffered_reads:
        if read.is_unmapped:
            continue
        # Only primary alignments count as independent sequencing evidence --
        # a supplementary/split segment crossing the region is the same
        # strand/fragment as its primary, not a second, independent mate.
        #if read.is_supplementary or read.is_secondary:
        #    continue
        if read.is_read1:
            overlap_mate[read.query_name].add(1)
        elif read.is_read2:
            overlap_mate[read.query_name].add(2)
        # else: unpaired / single-end -- no mate distinction possible,
        # so it can never satisfy the "both mates cross" condition below.

    def crosses_as_pair(qname: str) -> bool:
        return overlap_mate[qname] == {1, 2}

    # ---- Pass 1b: write every alignment overlapping the region, split or not ----
    print(f"Writing overlapping alignments to {out_bam}...")
    for read in buffered_reads:
        if read.is_unmapped:
            continue
#        if read.query_name != "M02063:376:000000000-DYVFW:1:1101:22715:19471":
#            continue

        XO = read.get_tag('XO') if read.has_tag('XO') else 0
        if overlaps_region(read) and crosses_as_pair(read.query_name):
            XO |= Tag_bits['OVERLAP']
        read = tag_read(read, XO)

        allele = get_allele_at_base(read, target_pos)
        indels = get_indels_at_base(read, target_pos)
        signed_len, summary_str, max_del, max_ins = summarise_indels(indels)

        # Defer the allele_counts decision -- other alignments of this same
        # read may still be coming later in fetch order (e.g. a duplication
        # where both split segments independently cover target_pos).
        if read.reference_start <= target_pos < read.reference_end:
            covering[read.query_name].append(
                (read.is_supplementary, read.reference_start, allele)
            )

        read.set_tag("XA", allele,      value_type="Z")
        read.set_tag("XI", signed_len,  value_type="i")
        read.set_tag("XD", max_del,     value_type="i")
        read.set_tag("XN", max_ins,     value_type="i")
        read.set_tag("XS", summary_str, value_type="Z")

        overlap_bam.write(calculate_softclip_lengths(read))
        written_records.add((read.query_name, read.reference_start))

        if read.has_tag('SA'):
            for loc in parse_sa_tag(read.get_tag('SA')):
                # SA tag positions are 1-based per the SAM spec.
                pending[read.query_name].append((loc['rname'], loc['pos'] - 1, 'SPLIT'))

        if read.is_paired and not read.mate_is_unmapped and read.next_reference_name is not None:
            pending[read.query_name].append(
                (read.next_reference_name, read.next_reference_start, 'MATE')
            )

    overlap_bam.close()
    # ---- Deterministically resolve one allele call per physical read ----
    # Preference order: non-supplementary (primary) alignment first; among
    # remaining ties, leftmost reference_start. Independent of fetch order.
    for qname, candidates in covering.items():
        _, _, allele = min(candidates, key=lambda c: (c[0], c[1]))
        allele_counts[allele] += 1
        XA_dict[qname] = allele

    # ---- Build the set of "other alignments" that still need fetching ----
    # Anything whose target (chrom, pos) matches an already-written record
    # for the same read is dropped here (e.g. the "both segments already
    # overlap" case needs no separate mate fetch).
    to_check      = defaultdict(list)
    scan_by_chrom = defaultdict(set)
    for qname, targets in pending.items():
        for c, p, kind in targets:
            if c is None or p is None or p < 0:
                continue
            if (qname, p) in written_records:
                continue
            to_check[qname].append((c, p, kind))
            scan_by_chrom[c].add(p)

    # ---- Pass 2: fetch and write the outstanding mates / split partners ----
    WINDOW = 300
    n_targets = sum(len(v) for v in to_check.values())
    print(f"Fetching {n_targets} additional alignment(s) for {len(to_check)} read(s) "
          f"and writing to {mates_bam}...")

    for scan_chrom, positions in scan_by_chrom.items():
        for pos in sorted(positions, reverse=True):
            try:
                for cand in bam.fetch(scan_chrom, max(0, pos - WINDOW), pos + WINDOW):
                    if cand.is_unmapped:
                        continue
                    key = (cand.query_name, cand.reference_start)
                    if key in written_records:
                        continue
                    if cand.query_name not in to_check:
                        continue

                    remaining = to_check[cand.query_name]
                    for i, (exp_chrom, exp_pos, kind) in enumerate(remaining):
                        if cand.reference_name == exp_chrom and cand.reference_start == exp_pos:
                            XO = cand.get_tag('XO') if cand.has_tag('XO') else 0
                            XO |= Tag_bits['NON_SPANNING_MATE'] if kind == 'MATE' else Tag_bits['SPLITX']
                            cand = tag_read(cand, XO)
                            cand.set_tag("XA", XA_dict.get(cand.query_name, ""), value_type="Z")
                            mate.write(calculate_softclip_lengths(cand))
                            written_records.add(key)
                            del remaining[i]
                            if not remaining:
                                del to_check[cand.query_name]
                            break
            except pysam.utils.SamtoolsError:
                pass

    mate.close()
    bam.close()

    n_missing = sum(len(v) for v in to_check.values())
    print(f"Done processing mates. {n_missing} expected alignment(s) across "
          f"{len(to_check)} read(s) not found.")

    return allele_counts


def _write_indel_stats(allele_counts: defaultdict, out_bam: str) -> None:
    """Write indel statistics file alongside the output BAM."""
    indelstat_output = f"{out_bam}.indelstats.txt"
    total = sum(allele_counts.values())
    with open(indelstat_output, 'w') as f:
        for allele, count in sorted(allele_counts.items(), key=lambda kv: kv[1], reverse=True):
            if   allele[0] == 'R': allele_type = 'ref'
            elif allele[0] == 'D': allele_type = 'del'
            elif allele[0] == 'I': allele_type = 'ins'
            elif allele[0] == 'M': allele_type = 'snp'
            else:                  allele_type = 'other'
            f.write(f"{allele}\t{count}\t{count/total:.4f}\t{allele_type}\n")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python extract_mates.py <input_bam> <extracted_bam> <mate_bam> <chrom> <start> <end>")
        sys.exit(1)

    input_bam     = sys.argv[1]
    extracted_bam = sys.argv[2]
    mate_bam      = sys.argv[3]
    chrom         = sys.argv[4]
    start         = int(sys.argv[5])
    end           = int(sys.argv[6])

    allele_counts = extract_mates(input_bam, extracted_bam, mate_bam, chrom, start, end)
    _write_indel_stats(allele_counts, extracted_bam)
