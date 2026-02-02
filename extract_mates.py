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
    overlap_bam = pysam.AlignmentFile(extracted_bam, "wb", header=header)
    mate = pysam.AlignmentFile(mate_bam, "wb", header=header)
    
    # Pass 1: Collect all reads in region and query mate positions
    mate_regions = defaultdict(set)   # unique (chrom, pos) tuples to query
    
    # bam header is copied to extracted_bam
    for read in bam.fetch(chrom, start-1, end):
        if read.is_unmapped:
            #print(f"Skipping unmapped read {read.query_name}")
            continue
        #print(f"mapped read {read.query_name} at {read.reference_name}:{read.reference_start} with mate at {read.next_reference_name}:{read.next_reference_start}")
        mate_chrom = read.next_reference_name
        mate_pos = read.next_reference_start
        if read.query_name not in mate_regions.keys():
            mate_regions[read.query_name]=(mate_chrom, mate_pos)
        else:
            read.set_tag('RG', 'NON_SPANNING_MATE')
            del mate_regions[read.query_name]
        overlap_bam.write(read) 
    
    overlap_bam.close()
    scan_positions = defaultdict(set)
    scan_positions_count = defaultdict(int)
    reads_to_check = set(mate_regions.keys())
    for _, (mate_chrom, mate_pos) in mate_regions.items():
        scan_positions[mate_chrom].add(mate_pos)
        scan_positions_count[(mate_chrom, mate_pos)] += 1
        #print (f"Need to find mate for read at {mate_chrom}:{mate_pos} read={R} count={scan_positions_count[(mate_chrom, mate_pos)]}")

    # Check each mate position
    for mate_chrom in scan_positions.keys(): 
        for mate_pos in sorted(scan_positions[mate_chrom]):
            window = 1
            #print(f"Fetching mates for {mate_chrom}:{mate_pos}")

            try:
                for mate_read in bam.fetch(mate_chrom, max(0, mate_pos - window), mate_pos + window):
                    if mate_read.query_name in reads_to_check:
                        scan_positions_count[(mate_chrom, mate_pos)] -= 1
                        mate_read.set_tag('RG', 'NON_SPANNING_MATE')
                        mate.write(mate_read)
                        reads_to_check.remove(mate_read.query_name)
                        #print(f"Written mate for read {mate_read.query_name} at {mate_chrom}:{mate_pos} with {scan_positions_count[(mate_chrom, mate_pos)]} remaining")
            except pysam.utils.SamtoolsError:
                pass
        
    mate.close()

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