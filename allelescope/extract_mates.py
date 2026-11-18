import hashlib
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


def get_allele_at_base(read, target_pos, ref_base=None):
    """
    Returns a descriptive allele string for what the read shows
    at target_pos (0-based reference coordinate).

    ``ref_base`` is the reference base at target_pos. It is the same for
    every read, so callers should look it up once; without it the base is
    parsed from each read's MD tag, which is slower and wrong when the BAM
    has no MD tags.

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
                    if ref_base is None:
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
    'F1R2':             0x1000,
    'R1F2':             0x2000,
    'F1F2':             0x4000,
    'R1R2':             0x8000,
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


# Orientations of a normal inward-facing (FR) pair. Which mate lies on the
# forward strand depends only on which strand of the fragment was read first,
# so both are equally "proper".
PROPER_ORIENTATIONS = ("F1R2", "R1F2")


def tag_read(read, XO):
    if read.has_tag('SA'):
        XO |= Tag_bits['SPLIT']
        #return read
    orientation = get_pair_orientation(read)
    if read.is_paired:
        read.set_tag('XR', orientation)
        if orientation in Tag_bits:
            XO |= Tag_bits[orientation]
        if read.is_proper_pair and orientation in PROPER_ORIENTATIONS:
            XO |= Tag_bits['PROPER']
        else:
            XO |= Tag_bits['IMPROPER']
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

# Targets on one sequence closer together than this (bp) share a fetch.
# Reading a short gap of alignments costs less than an extra index seek.
MATE_FETCH_MERGE_GAP = 1000


def _reference_base(ref_fasta, chrom, pos0):
    """Upper-case reference base at a 0-based position, or None if unavailable."""
    if not ref_fasta:
        return None
    try:
        with pysam.FastaFile(str(ref_fasta)) as fasta:
            base = fasta.fetch(chrom, pos0, pos0 + 1).upper()
    except (OSError, ValueError, KeyError, IndexError):
        return None
    return base if len(base) == 1 and base in "ACGTN" else None


def _base_quality_at(read, target_pos):
    """Base quality of the read at reference position target_pos (0-based),
    or None where the read has no base there (deletion, skip)."""
    quals = read.query_qualities
    if quals is None:
        return None
    qpos, rpos = 0, read.reference_start
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8):
            if rpos <= target_pos < rpos + length:
                return quals[qpos + target_pos - rpos]
            qpos += length
            rpos += length
        elif op in (1, 4):
            qpos += length
        elif op in (2, 3):
            if rpos <= target_pos < rpos + length:
                return None
            rpos += length
    return None


# ---- Read evidence per allele ---------------------------------------------
# Distances in bp. A base this close to either end of the aligned part of a
# read counts as "near the end"; a soft clip of at least CLIP_MIN_BP whose
# edge is this close to the site counts as "clipped near the site".
NEAR_END_BP = 10
CLIP_NEAR_BP = 10
CLIP_MIN_BP = 5


def _site_query_pos(read, target_pos):
    """
    (query index, base quality) of the read at reference position target_pos
    (0-based). Inside a deletion the index is that of the first base after it
    and the quality is None. (None, None) when the read does not reach there.
    """
    quals = read.query_qualities
    qpos, rpos = 0, read.reference_start
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8):
            if rpos <= target_pos < rpos + length:
                q = qpos + target_pos - rpos
                return q, (quals[q] if quals is not None else None)
            qpos += length
            rpos += length
        elif op in (1, 4):
            qpos += length
        elif op in (2, 3):
            if rpos <= target_pos < rpos + length:
                return qpos, None
            rpos += length
    return None, None


def _allele_edit_length(allele):
    """How much the allele itself adds to the read's NM (edit distance)."""
    if not allele:
        return 0
    kind = allele[0]
    if kind == "M":
        return 1
    if kind == "D":
        try:
            return int(allele[1:])
        except ValueError:
            return 0
    if kind == "I":
        return len(allele) - 1
    return 0


def _read_evidence(read, target_pos, allele):
    """
    What one read says about how trustworthy its call at target_pos is.
    Must be called before the read's tags are rewritten (XA in particular,
    which bwa uses for alternative hits and this module reuses for alleles).

    Returns a tuple:
      (is_reverse, mapq, base_quality, end_distance, other_edits,
       ambiguous, chimeric, clipped_near, insert_size, f1r2, duplicate,
       aligned_start, aligned_end)
    aligned_start and aligned_end are 1-based and inclusive: the first and
    last reference bases of the aligned part (used to find primer sequence).
    other_edits is the NM tag minus the allele's own edits (None without NM);
    insert_size is None unless both mates map to the same sequence;
    f1r2 is None for unpaired reads.
    """
    q, bq = _site_query_pos(read, target_pos)
    end_distance = None
    if q is not None:
        qas, qae = read.query_alignment_start, read.query_alignment_end
        end_distance = max(0, min(q - qas, qae - 1 - q))

    other_edits = None
    if read.has_tag("NM"):
        other_edits = max(0, int(read.get_tag("NM")) - _allele_edit_length(allele))

    # Mapped elsewhere about as well: bwa lists alternative hits in XA, and
    # reports the best other score in XS (an integer, unlike this module's XS).
    ambiguous = read.mapping_quality == 0 or read.is_secondary or read.has_tag("XA")
    if not ambiguous and read.has_tag("AS") and read.has_tag("XS"):
        xs = read.get_tag("XS")
        if isinstance(xs, int) and xs >= read.get_tag("AS"):
            ambiguous = True

    chimeric = read.is_supplementary or read.has_tag("SA")

    clipped_near = False
    cig = read.cigartuples or ()

    def soft_clip(ops):
        n = 0
        for op, length in ops:          # clips sit at the ends: H then S
            if op == 4:
                n += length
            elif op != 5:
                break
        return n

    if soft_clip(cig) >= CLIP_MIN_BP and target_pos - read.reference_start <= CLIP_NEAR_BP:
        clipped_near = True
    if soft_clip(reversed(cig)) >= CLIP_MIN_BP and (read.reference_end - 1) - target_pos <= CLIP_NEAR_BP:
        clipped_near = True

    insert_size = None
    f1r2 = None
    if read.is_paired:
        f1r2 = read.is_read1 != read.is_reverse
        if (not read.mate_is_unmapped and read.template_length
                and read.reference_id == read.next_reference_id):
            insert_size = abs(read.template_length)

    return (read.is_reverse, read.mapping_quality, bq, end_distance, other_edits,
            ambiguous, chimeric, clipped_near, insert_size, f1r2, read.is_duplicate,
            read.reference_start + 1, read.reference_end)


def _new_allele_stats():
    return {"fwd": 0, "rev": 0, "bq_sum": 0, "bq_n": 0, "mq_sum": 0,
            "end_dist": [], "nm_sum": 0, "nm_n": 0, "mq0": 0, "ambiguous": 0,
            "chimeric": 0, "clipped": 0, "isize": [], "f1r2": 0, "f2r1": 0, "dup": 0,
            "starts": defaultdict(int), "ends": defaultdict(int)}


def _add_evidence(st, ev):
    (rev, mq, bq, end_distance, other_edits, ambiguous, chimeric,
     clipped_near, insert_size, f1r2, dup, aligned_start, aligned_end) = ev
    st["starts"][aligned_start] += 1
    st["ends"][aligned_end] += 1
    st["rev" if rev else "fwd"] += 1
    st["mq_sum"] += mq
    if mq == 0:
        st["mq0"] += 1
    if bq is not None:
        st["bq_sum"] += bq
        st["bq_n"] += 1
    if end_distance is not None:
        st["end_dist"].append(end_distance)
    if other_edits is not None:
        st["nm_sum"] += other_edits
        st["nm_n"] += 1
    st["ambiguous"] += ambiguous
    st["chimeric"] += chimeric
    st["clipped"] += clipped_near
    if insert_size is not None:
        st["isize"].append(insert_size)
    if f1r2 is not None:
        st["f1r2" if f1r2 else "f2r1"] += 1
    st["dup"] += dup


def _median(values):
    if not values:
        return None
    s = sorted(values)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


def allele_evidence_summary(allele_stats):
    """
    Per-allele summary written to <bam>.evidence.json and sent to the browser.
    Counts are numbers of reads (one per fragment, as for the allele counts).
    """
    out = {}
    for allele, st in allele_stats.items():
        n = st["fwd"] + st["rev"]
        dist = st["end_dist"]
        out[allele] = {
            "n": n,
            "end_median": _median(dist),
            "near_end": sum(1 for d in dist if d < NEAR_END_BP),
            "end_n": len(dist),
            "nm_mean": round(st["nm_sum"] / st["nm_n"], 2) if st["nm_n"] else None,
            "mq0": st["mq0"],
            "ambiguous": st["ambiguous"],
            "chimeric": st["chimeric"],
            "clipped": st["clipped"],
            "isize_median": _median(st["isize"]),
            "isize_n": len(st["isize"]),
            "f1r2": st["f1r2"],
            "f2r1": st["f2r1"],
            "dup": st["dup"],
            # Where the reads' aligned parts start and end: {position: reads}.
            # With a primer scheme these show which reads have primer
            # sequence at the site (annotation.Annotation.primer_reads).
            "starts": {str(k): v for k, v in sorted(st["starts"].items())},
            "ends": {str(k): v for k, v in sorted(st["ends"].items())},
        }
    return out


def _name_rank(qname: str) -> bytes:
    """Fixed pseudo-random sort key for a read name (same on every run)."""
    return hashlib.blake2b(qname.encode("utf-8", "surrogateescape"), digest_size=8,
                           person=b"allelescope").digest()


def _choose_read_names(names, max_pairs):
    """The read names to keep, or None to keep all."""
    if not max_pairs or max_pairs <= 0 or len(names) <= max_pairs:
        return None
    return set(sorted(names, key=_name_rank)[:max_pairs])


def _merge_positions(positions, gap):
    """
    Merge 0-based start positions into half-open fetch intervals.

    Positions less than ``gap`` bp apart end up in the same interval; each
    interval ends one base after its last position, which is enough for
    ``fetch`` to return the records that start there.
    """
    intervals = []
    for pos in sorted(positions):
        if intervals and pos - (intervals[-1][1] - 1) <= gap:
            intervals[-1][1] = pos + 1
        else:
            intervals.append([pos, pos + 1])
    return [tuple(iv) for iv in intervals]


def extract_mates(
    src_bam: str,
    out_bam: str,
    mates_bam: str,
    chrom: str,
    start: int,
    end: int,
    ref_fasta: str = None,
    max_pairs: int = None,
    sample_info: dict = None,
    allele_stats: dict = None,
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
        max_pairs: keep at most this many read pairs (fragments) crossing
                   the site, chosen by a fixed hash of the read name; their
                   mates and split partners are kept with them. Allele
                   counts still use every read. None keeps everything.
        sample_info: optional dict, filled with "pairs_total" and
                   "pairs_kept".
        allele_stats: optional dict, filled per allele with strand counts,
                   base and mapping quality, and the read evidence summarised
                   by allele_evidence_summary(), over the same reads as the
                   allele counts (one read per fragment).
        ref_fasta: indexed reference FASTA. Used to look up the reference
                   base at the site once; without it each read's MD tag
                   is parsed instead.

    Returns:
        allele_counts: dict mapping allele string -> read count (one
        observation per physical read/fragment, deterministic).
    """
    allele_counts = defaultdict(int)

    site_ref_base = _reference_base(ref_fasta, chrom, start - 1)

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
    covering         = defaultdict(list)  # query_name -> [(is_supplementary, reference_start, allele, evidence), ...]

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
    # ---- Optional subsampling of the read pairs crossing the site ----
    # Allele counts use every read; only the chosen read pairs are tagged,
    # written and have their mates fetched. The choice is a fixed hash of the
    # read name, so both mates and all split segments of a fragment are kept
    # or dropped together, and the same reads are chosen every time.
    names = {read.query_name for read in buffered_reads if not read.is_unmapped}
    keep = _choose_read_names(names, max_pairs)
    if sample_info is not None:
        sample_info["pairs_total"] = len(names)
        sample_info["pairs_kept"] = len(names) if keep is None else len(keep)

    tagged_reads = []
    for read in buffered_reads:
        if read.is_unmapped:
            continue
        if keep is not None and read.query_name not in keep:
            # Dropped from the display: only its allele call is needed.
            if read.reference_start <= target_pos < read.reference_end:
                allele = get_allele_at_base(read, target_pos, site_ref_base)
                covering[read.query_name].append(
                    (read.is_supplementary, read.reference_start, allele,
                     _read_evidence(read, target_pos, allele) if allele_stats is not None else None)
                )
            continue

        XO = read.get_tag('XO') if read.has_tag('XO') else 0
        if overlaps_region(read) and crosses_as_pair(read.query_name):
            XO |= Tag_bits['OVERLAP']
        read = tag_read(read, XO)

        allele = get_allele_at_base(read, target_pos, site_ref_base)
        indels = get_indels_at_base(read, target_pos)
        signed_len, summary_str, max_del, max_ins = summarise_indels(indels)

        # Defer the allele_counts decision -- other alignments of this same
        # read may still be coming later in fetch order (e.g. a duplication
        # where both split segments independently cover target_pos).
        if read.reference_start <= target_pos < read.reference_end:
            covering[read.query_name].append(
                (read.is_supplementary, read.reference_start, allele,
                 _read_evidence(read, target_pos, allele) if allele_stats is not None else None)
            )

        read.set_tag("XA", allele,      value_type="Z")
        read.set_tag("XI", signed_len,  value_type="i")
        read.set_tag("XD", max_del,     value_type="i")
        read.set_tag("XN", max_ins,     value_type="i")
        read.set_tag("XS", summary_str, value_type="Z")
        tagged_reads.append(read)

    for read in tagged_reads:
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
        _, _, allele, evidence = min(candidates, key=lambda c: (c[0], c[1]))
        allele_counts[allele] += 1
        XA_dict[qname] = allele
        if allele_stats is not None:
            if allele not in allele_stats:
                allele_stats[allele] = _new_allele_stats()
            # When both mates cross the site, the site is only "near a read
            # end" if it is near the end of each of them: use the read where
            # it lies furthest from an end.
            dists = [c[3][3] for c in candidates
                     if c[2] == allele and c[3] is not None and c[3][3] is not None]
            if len(dists) > 1 and max(dists) != evidence[3]:
                evidence = evidence[:3] + (max(dists),) + evidence[4:]
            _add_evidence(allele_stats[allele], evidence)

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
    # Every target is an exact (chrom, start) taken from RNEXT/PNEXT or the
    # SA tag, and a record is only accepted if its reference_start equals it.
    # Around a deep site the targets are thousands of positions within a few
    # hundred bp, so fetching a window per target decoded the same reads over
    # and over (cost ~ number of targets x depth). Instead the targets on each
    # sequence are merged into a few intervals, each record is read once, and
    # its start position (an int) is checked before its name is decoded.
    n_targets = sum(len(v) for v in to_check.values())
    print(f"Fetching {n_targets} additional alignment(s) for {len(to_check)} read(s) "
          f"and writing to {mates_bam}...")

    for scan_chrom, positions in scan_by_chrom.items():
        for fetch_start, fetch_end in _merge_positions(positions, MATE_FETCH_MERGE_GAP):
            if not to_check:
                break
            try:
                candidates = bam.fetch(scan_chrom, fetch_start, fetch_end)
                for cand in candidates:
                    if cand.reference_start not in positions or cand.is_unmapped:
                        continue
                    qname = cand.query_name
                    remaining = to_check.get(qname)
                    if remaining is None:
                        continue
                    key = (qname, cand.reference_start)
                    if key in written_records:
                        continue

                    for i, (exp_chrom, exp_pos, kind) in enumerate(remaining):
                        if cand.reference_start == exp_pos and cand.reference_name == exp_chrom:
                            XO = cand.get_tag('XO') if cand.has_tag('XO') else 0
                            XO |= Tag_bits['NON_SPANNING_MATE'] if kind == 'MATE' else Tag_bits['SPLITX']
                            cand = tag_read(cand, XO)
                            cand.set_tag("XA", XA_dict.get(qname, ""), value_type="Z")
                            mate.write(calculate_softclip_lengths(cand))
                            written_records.add(key)
                            del remaining[i]
                            if not remaining:
                                del to_check[qname]
                            break
            except (ValueError, pysam.utils.SamtoolsError):
                # Sequence named in RNEXT/SA but absent from the BAM header
                # (ValueError), or an htslib read error: skip these targets.
                pass

    mate.close()
    bam.close()

    n_missing = sum(len(v) for v in to_check.values())
    print(f"Done processing mates. {n_missing} expected alignment(s) across "
          f"{len(to_check)} read(s) not found.")

    return allele_counts


def _write_indel_stats(allele_counts: defaultdict, out_bam: str, allele_stats: dict = None) -> None:
    """
    Write the allele table alongside the output BAM, one allele per line:
        allele  count  fraction  type  [fwd  rev  mean_BQ  mean_MAPQ]
    The bracketed columns are written when allele_stats is given; mean_BQ is
    empty for alleles without a base at the site (deletions).
    """
    indelstat_output = f"{out_bam}.indelstats.txt"
    total = sum(allele_counts.values())
    with open(indelstat_output, 'w', encoding='utf-8') as f:
        for allele, count in sorted(allele_counts.items(), key=lambda kv: kv[1], reverse=True):
            if   allele[0] == 'R': allele_type = 'ref'
            elif allele[0] == 'D': allele_type = 'del'
            elif allele[0] == 'I': allele_type = 'ins'
            elif allele[0] == 'M': allele_type = 'snp'
            else:                  allele_type = 'other'
            line = f"{allele}\t{count}\t{count/total:.4f}\t{allele_type}"
            if allele_stats is not None and allele in allele_stats:
                st = allele_stats[allele]
                bq = f"{st['bq_sum'] / st['bq_n']:.1f}" if st["bq_n"] else ""
                mq = f"{st['mq_sum'] / max(1, st['fwd'] + st['rev']):.1f}"
                line += f"\t{st['fwd']}\t{st['rev']}\t{bq}\t{mq}"
            f.write(line + "\n")



def _write_evidence(allele_stats: dict, out_bam: str) -> None:
    """Write allele_evidence_summary() to <out_bam>.evidence.json."""
    import json
    with open(f"{out_bam}.evidence.json", "w", encoding="utf-8") as f:
        json.dump(allele_evidence_summary(allele_stats), f)

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
