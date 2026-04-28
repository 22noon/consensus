import pysam
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




###########  End of indel tagger code, start of mate extraction code ############
Tag_bits = {
    'OVERLAP': 0x1,
    'NON_SPANNING_MATE': 0x2,
    'SPLITX': 0x4,
    'SPLIT': 0x8,
    'PROPER': 0x10,
    'IMPROPER': 0x20,
    'READ1': 0x40,
    'READ2': 0x80,
    'SECONDARY': 0x100,
    'QCFAIL': 0x200,
    'DUPLICATE': 0x400,
    'SUPPLEMENTARY': 0x800
}

def parse_sa_tag(sa_string):
#Parse SA tag to get locations of other alignments
    parts = []
    for entry in sa_string.rstrip(';').split(';'):
        rname, pos, strand, cigar, mapq, nm = entry.split(',')
        parts.append({
            'rname': rname,
            'pos': int(pos),
            'strand': strand,
            'cigar': cigar,
            'mapq': int(mapq),
            'nm': int(nm)
        })
    return parts

def calculate_softclip_lengths(read):
    cigar_tuples = read.cigartuples
    left_sc, right_sc = 0, 0
    if  cigar_tuples:
        left_sc = cigar_tuples[0][1] if cigar_tuples[0][0] == 4 else 0
        right_sc = cigar_tuples[-1][1] if cigar_tuples[-1][0] == 4 else 0
    read.set_tag("SL", left_sc, value_type="i")  # Left soft-clip
    read.set_tag("SR", right_sc, value_type="i")  # Right soft-clip
    return read


# Usage
def process_mates_optimized(input_bam, extracted_bam, mate_bam, chrom, start,end):
    bam = pysam.AlignmentFile(input_bam, "rb")
    header = bam.header.to_dict()
    # Add read group to header
    if 'RG' not in header:
        header['RG'] = []
    rg_ids = [
        'NON_SPANNING_MATE',
        'F1R2','R2F1','F2R1','R1F2','F1F2','R1R2',
        'SPLIT','SPLITX','OVERLAP','IMPROPER','UNPAIRED','PROPER'
    ]
    existing_ids = {rg.get('ID') for rg in header['RG']}
    for rg_id in rg_ids:
        if rg_id not in existing_ids:
            header['RG'].append({'ID': rg_id, 'SM': 'sample', 'PL': 'ILLUMINA'})
            existing_ids.add(rg_id)
    overlap_bam = pysam.AlignmentFile(extracted_bam, "wb", header=header)
    mate = pysam.AlignmentFile(mate_bam, "wb", header=header)
    
    # Pass 1: Collect all reads in region and query mate positions
    print(f"Scanning region {chrom}:{start}-{end} for overlapping reads and their mates...")
    mate_regions = defaultdict()   # unique (chrom, pos) tuples to query
    overlap_count = defaultdict(int)
    overlaps = set()
    split_reads = defaultdict(set) # read_name -> set of (chrom, pos) for split alignments 
    
    for read in bam.fetch(chrom, start-1, end):
        if read.is_unmapped:
            continue
        if not read.has_tag('SA'): #Read is not split
           overlap_count[read.query_name] += 1
    overlaps = {r for r, count in overlap_count.items() if count > 1}

#Second pass: Write overlapping reads to extracted_bam and collect mate and split positions for non-overlapping reads
    print(f"Writing overlapping reads to {extracted_bam} and collecting mate positions...")
    for read in bam.fetch(chrom, start-1, end):
        if read.is_unmapped:
            continue

        XO = read.get_tag('XO') if read.has_tag('XO') else 0 
        if read.has_tag('SA'): #Read is split
            other_locations = parse_sa_tag(read.get_tag('SA'))
            for loc in other_locations:
                if read.reference_name != loc['rname']:
                    continue
                split_reads[read.query_name].add((loc['rname'], loc['pos']))
        else:
            mate_chrom = read.next_reference_name
            mate_pos = read.next_reference_start
            if read.query_name in overlaps:
                XO |= Tag_bits['OVERLAP']
            if read.query_name not in mate_regions.keys():
                mate_regions[read.query_name]=(mate_chrom, mate_pos)
            else:
                del mate_regions[read.query_name]
        read = tag_read(read,XO)

        # Indel tags
        allele = get_allele_at_base(read, start-1)
        indels = get_indels_at_base(read, start-1)
        signed_len, summary_str, max_del, max_ins = summarise_indels(indels)

        read.set_tag("XA", allele,     value_type="Z")
        read.set_tag("XI", signed_len, value_type="i")
        read.set_tag("XD", max_del,    value_type="i")
        read.set_tag("XN", max_ins,    value_type="i")
        read.set_tag("XS", summary_str,value_type="Z")
# Then add two tags
        overlap_bam.write(calculate_softclip_lengths(read)) 
    
    overlap_bam.close()
    scan_positions = defaultdict(set)
    scan_positions_count = defaultdict(int)
    mates_to_check = set(mate_regions.keys())
    splits_to_check = set(split_reads.keys())

    # Populate scan_positions and counts from mate_regions and write mates_to_check to a file
    for read_name, mate_info in mate_regions.items():
        mate_chrom, mate_pos = mate_info
        if mate_chrom is None or mate_pos is None or mate_pos < 0:
            continue
        scan_positions[mate_chrom].add(mate_pos)
        scan_positions_count[(mate_chrom, mate_pos)] += 1

    for read_name, split_info in split_reads.items():
        for mate_chrom, mate_pos in split_info:
            if mate_chrom is None or mate_pos is None or mate_pos < 0:
                print(f"Skipping split read {read_name} with invalid mate info: {mate_chrom}:{mate_pos}")
                continue
            scan_positions[mate_chrom].add(mate_pos)
            scan_positions_count[(mate_chrom, mate_pos)] += 1


    print(f"Fetching mates for {len(scan_positions_count)} positions and writing to {mate_bam}...")
    for mate_chrom in scan_positions.keys(): 
        for mate_pos in sorted(scan_positions[mate_chrom], reverse=True):
            #print(f"\tFetching mates for {mate_chrom}:{mate_pos} with {scan_positions_count[(mate_chrom, mate_pos)]} reads to check")
            window = 300 
            if scan_positions_count[(mate_chrom, mate_pos)] == 0:
                continue

            try:
                for mate_read in bam.fetch(mate_chrom, max(0, mate_pos - window), mate_pos + window):
                    XO = mate_read.get_tag('XO') if mate_read.has_tag('XO') else 0 
                    if mate_read.query_name in mates_to_check:
                        expected_chrom, expected_pos = mate_regions.get(mate_read.query_name, (None, None))
                        if expected_chrom is None or expected_pos is None:
                            print(f"Skipping read {mate_read.query_name} with missing mate info")
                        else:
                            if mate_read.reference_name != expected_chrom or mate_read.reference_start != expected_pos:
                                continue
                        scan_positions_count[mate_regions[mate_read.query_name]] -= 1
                        XO |= Tag_bits['NON_SPANNING_MATE']
                        mate_read = tag_read(mate_read,XO)
                        mate.write(calculate_softclip_lengths(mate_read))
                        mates_to_check.remove(mate_read.query_name)
                    if mate_read.query_name in splits_to_check:
                        split_info = split_reads.get(mate_read.query_name, set())
                        #print(f"Checking split info for {mate_read.query_name}: {split_info} against {mate_read.reference_name}:{mate_read.reference_start +1}")
                        if (mate_read.reference_name, mate_read.reference_start +1 ) in split_info:
                            mate_read.set_tag('XO', XO | Tag_bits['SPLITX'])
                            #mate_read = tag_read(mate_read)
                            mate.write(calculate_softclip_lengths(mate_read))
                            splits_to_check.remove(mate_read.query_name)
            except pysam.utils.SamtoolsError:
                pass
        
    mate.close()
    print(f"Done processing mates. Remaining reads without found mates: {len(mates_to_check)}")

def get_pair_orientation(read):
    """
    Determine the orientation of a paired read.
    Returns: F1R2, R2F1, F2R1, R1F2, F1F2, R1R2, etc.
    """
    if read.is_unmapped or read.mate_is_unmapped:
        return "UNMAPPED"
    
    # Check if this is read1 or read2
    is_read1 = read.is_read1
    
    # Check strand orientation
    read_reverse = read.is_reverse
    mate_reverse = read.mate_is_reverse
    
    # Determine orientation string
    if is_read1:
        # This is read 1
        read1_strand = 'R1' if read_reverse else 'F1'
        read2_strand = 'R2' if mate_reverse else 'F2'
        return f"{read1_strand}{read2_strand}"
    else:
        # This is read 2
        read2_strand = 'R2' if read_reverse else 'F2'
        read1_strand = 'R1' if mate_reverse else 'F1'
        return f"{read1_strand}{read2_strand}"

def tag_read(read,XO):
    if read.has_tag('SA'):
        XO |= Tag_bits['SPLIT']
        return read
    orientation = get_pair_orientation(read)
    if read.is_paired:
        read.set_tag('XR', orientation)
        if orientation != "F1R2" or not read.is_proper_pair:
            XO |= Tag_bits['IMPROPER']
        else:
            XO |= Tag_bits['PROPER']
    #else:
    #    read.set_tag('RG', 'UNPAIRED')
    read.set_tag('XO', XO)
    return read

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 7:
        print("Usage: python extract_mates.py <input_bam> <extracted_bam> <mate_bam> <chrom> <start> <end>")
        sys.exit(1)
    
    input_bam = sys.argv[1]
    extracted_bam = sys.argv[2]
    mate_bam = sys.argv[3]
    chrom = sys.argv[4]
    start = int(sys.argv[5])
    end = int(sys.argv[6])

    process_mates_optimized(input_bam, extracted_bam, mate_bam, chrom, start, end)
