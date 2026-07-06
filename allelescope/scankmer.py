from collections import defaultdict

# INPUTS
#read = "GCCCAAACCCTCCCAACGAACCCCCCCCCCTCCCCCCCGAACTCAACCCACGG"
#read = "CCGTGGGTTGAGTTCGGGGGGGAGGGGGGGGGGTTCGTTGGGAGGGTTTGGGC"
read = "CACCACTCCCCTCCTGCCACACCGGCCCGCAAAGGCCCACAGCAGGGCAAGACGCCCGACACCAGCATCAACCCATTCGAACACCAGGGCCCAAACCCTCCCAACGAACCCCCCCCCCTCCCCCCCGAACTCAACCCACGGAGCACCCAGC"
k = 19
ref = open("region.fa").read().upper().replace("\n", "")

# Build k-mer index for reference
ref_index = defaultdict(int)
for i in range(len(ref) - k + 1):
    kmer = ref[i:i+k]
    ref_index[kmer] += 1

# Scan read
results = []
for i in range(len(read) - k + 1):
    kmer = read[i:i+k]
    count = ref_index[kmer]
    results.append((i, kmer, count))

# Print
for pos, kmer, count in results:
    print(f"{pos:3d} {kmer}  count={count}")
