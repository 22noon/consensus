import pysam
from collections import defaultdict

def process_mates_optimized(input_bam, extracted_bam, mate_bam, chrom, start,end):
    bam = pysam.AlignmentFile(input_bam, "rb")
    header = bam.header.to_dict()
    # Add read group to header
    if 'RG' not in header:
        header['RG'] = []
    if 'NON_SPANNING_MATE' not in [rg['ID'] for rg in header['RG']]:
        header['RG'].append({
            'ID': 'NON_SPANNING_MATE',
            'SM': 'sample',
            'PL': 'ILLUMINA'
        })
    if 'OVERLAP' not in [rg['ID'] for rg in header['RG']]:
        header['RG'].append({
            'ID': 'OVERLAP',
            'SM': 'sample',
            'PL': 'ILLUMINA'
        })
    if 'IMPROPER' not in [rg['ID'] for rg in header['RG']]:
        header['RG'].append({
            'ID': 'IMPROPER',
            'SM': 'sample',
            'PL': 'ILLUMINA'
        })
    if 'UNPAIRED' not in [rg['ID'] for rg in header['RG']]:
        header['RG'].append({
            'ID': 'UNPAIRED',
            'SM': 'sample',
            'PL': 'ILLUMINA'
        })
    if 'PROPER' not in [rg['ID'] for rg in header['RG']]:
        header['RG'].append({
            'ID': 'PROPER',
            'SM': 'sample',
            'PL': 'ILLUMINA'
        })
    overlap_bam = pysam.AlignmentFile(extracted_bam, "wb", header=header)
    mate = pysam.AlignmentFile(mate_bam, "wb", header=header)
    
    # Pass 1: Collect all reads in region and query mate positions
    print(f"Scanning region {chrom}:{start}-{end} for overlapping reads and their mates...")
    mate_regions = defaultdict(set)   # unique (chrom, pos) tuples to query
    overlap_count = defaultdict(int)
    overlaps = set()
    
    for read in bam.fetch(chrom, start-1, end):
        if read.is_unmapped:
            continue
        overlap_count[read.query_name] += 1
    overlaps = {r for r, count in overlap_count.items() if count > 1}

    # bam header is copied to extracted_bam
    print(f"Writing overlapping reads to {extracted_bam} and collecting mate positions...")
    for read in bam.fetch(chrom, start-1, end):
        if read.is_unmapped:
            continue
        mate_chrom = read.next_reference_name
        mate_pos = read.next_reference_start
        if read.query_name in overlaps:
            read.set_tag('RG', 'OVERLAP')
        if read.query_name not in mate_regions.keys():
            mate_regions[read.query_name]=(mate_chrom, mate_pos)
        else:
            #read.set_tag('RG', 'NON_SPANNING_MATE')
            del mate_regions[read.query_name]
        # Mark improperly paired reads with RG=IMPROPER (overrides any previous RG tag)
        read = tag_read(read)
        overlap_bam.write(read) 
    
    overlap_bam.close()
    scan_positions = defaultdict(set)
    scan_positions_count = defaultdict(int)
    reads_to_check = set(mate_regions.keys())
    # Populate scan_positions and counts from mate_regions and write reads_to_check to a file
    for read_name, mate_info in mate_regions.items():
        try:
            mate_chrom, mate_pos = mate_info
        except Exception:
            continue
        if mate_chrom is None or mate_pos is None or mate_pos < 0:
            continue
        scan_positions[mate_chrom].add(mate_pos)
        scan_positions_count[(mate_chrom, mate_pos)] += 1


    print(f"Fetching mates for {len(scan_positions_count)} positions and writing to {mate_bam}...")
    for mate_chrom in scan_positions.keys(): 
        for mate_pos in sorted(scan_positions[mate_chrom], reverse=True):
            print(f"\tFetching mates for {mate_chrom}:{mate_pos} with {scan_positions_count[(mate_chrom, mate_pos)]} reads to check")
            window = 300 
            if scan_positions_count[(mate_chrom, mate_pos)] == 0:
                continue
            #print(f"Fetching mates for {mate_chrom}:{mate_pos}")

            try:
                for mate_read in bam.fetch(mate_chrom, max(0, mate_pos - window), mate_pos + window):
                    if mate_read.query_name in reads_to_check:
                        expected_chrom, expected_pos = mate_regions.get(mate_read.query_name, (None, None))
                        #if mate_read.query_name == "M02063:330:000000000-LV7WG:1:2107:16849:15250":
                        #    print(f"Found mate read: {mate_read} position: {mate_read.reference_name}:{mate_read.reference_start} next_position: {mate_read.next_reference_name}:{mate_read.next_reference_start}")
                        #    print(f"expected position: {mate_regions[mate_read.query_name]} scan count: {scan_positions_count[(mate_chrom, mate_pos)]}")
                        if expected_chrom is None or expected_pos is None:
                            print(f"Skipping read {mate_read.query_name} with missing mate info")
                        else:
                            if mate_read.reference_name != expected_chrom or mate_read.reference_start != expected_pos:
                                continue
                        scan_positions_count[mate_regions[mate_read.query_name]] -= 1
                        #scan_positions_count[(mate_chrom, mate_pos)] -= 1
                        mate_read.set_tag('RG', 'NON_SPANNING_MATE')
                        mate_read = tag_read(mate_read)
                        mate.write(mate_read)
                        reads_to_check.remove(mate_read.query_name)
                        #print(f"Written mate for read {mate_read.query_name} at {mate_chrom}:{mate_pos} with {scan_positions_count[(mate_chrom, mate_pos)]} remaining")
            except pysam.utils.SamtoolsError:
                pass
        
    mate.close()
    print(f"Done processing mates. Remaining reads without found mates: {len(reads_to_check)}")

def tag_read(read):
    if read.is_paired:
        if not read.is_proper_pair:
            read.set_tag('RG', 'IMPROPER')
        else:
            read.set_tag('RG', 'PROPER')
    else:
        read.set_tag('RG', 'UNPAIRED')
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
