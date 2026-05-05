import pysam
from collections import defaultdict


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
                    # op == 0: use MD tag to determine match/mismatch
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
        md = read.get_tag("MD")
    except KeyError:
        return None

    import re
    ref_pos   = read.reference_start
    query_pos = 0

    # Walk CIGAR in parallel with MD to track reference base at target
    # MD only covers aligned (non-clipped) bases, so we need CIGAR too
    # to correctly advance query_pos through soft clips first.

    # Precompute: offset into the MD-covered positions
    # Simpler approach: use pysam's get_aligned_pairs with seq
    try:
        for qpos, rpos, ref_base in read.get_aligned_pairs(
                with_seq=True, matches_only=False):
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
    indels    = []
    ref_pos   = read.reference_start
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


def tag_bam(input_bam, output_bam, chrom, target_pos):
    """
    Writes per-read tags at target_pos (0-based):

      XA:Z  — allele string: D7 | IATT | RA | MT | S | .
      XI:i  — signed net indel length (+ ins, - del, 0 = ref/SNP)
      XD:i  — largest deletion length
      XN:i  — largest insertion length
      XS:Z  — full indel summary e.g. 'D7,I3'
    """
    allele_counts = defaultdict(int)

    with pysam.AlignmentFile(input_bam, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:

        for read in infile.fetch(chrom, target_pos, target_pos + 1):

            allele = get_allele_at_base(read, target_pos)
            indels = get_indels_at_base(read, target_pos)
            signed_len, summary_str, max_del, max_ins = summarise_indels(indels)

            read.set_tag("XA", allele,     value_type="Z")
            read.set_tag("XI", signed_len, value_type="i")
            read.set_tag("XD", max_del,    value_type="i")
            read.set_tag("XN", max_ins,    value_type="i")
            read.set_tag("XS", summary_str,value_type="Z")

            allele_counts[allele] += 1
            outfile.write(read)

    sorted_bam = output_bam.replace('.bam', '.sorted.bam')
    pysam.sort("-o", sorted_bam, output_bam)
    pysam.index(sorted_bam)

    print(f"\nAllele pileup at {chrom}:{target_pos + 1}")
    print(f"{'Allele':<12} {'Count':>6}  {'Bar'}")
    print("-" * 40)
    total = sum(allele_counts.values())
    for allele, count in sorted(allele_counts.items(),
                                 key=lambda x: -x[1]):
        pct = count / total * 100
        bar = '█' * int(pct / 2)
        print(f"{allele:<12} {count:>6}  {bar}  {pct:.1f}%")

    return sorted_bam
tag_bam("Datasets/RPV/alignment.markdup.bam", "tagged.bam", "MN632612.1", 4568)

