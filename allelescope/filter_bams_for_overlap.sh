#!/bin/bash
#VCF=$1
VCF=calls.norm.vcf.gz
 mkdir -p bam
# samtools collate -u -@ 10 -O alignment.bam tmp_prefix | \
# samtools fixmate -@ 10 -m -u - - | \
# samtools sort -@ 10 -u -O bam - | \
# samtools markdup -@ 10 - alignment.markdup.bam

# samtools index alignment.markdup.bam

# For each variant, extract reads overlapping that position
while IFS=$'\t' read -r CHROM POS; do
    echo "Processing ${CHROM}:${POS}"
    
    # samtools already filters to reads overlapping this position
    # path of the current running script (absolute)
    SCRIPT_PATH=$(dirname "$(readlink -f "${BASH_SOURCE[0]:-$0}")")
    python $SCRIPT_PATH/extract_mates.py alignment.markdup.bam bam/variant_${CHROM}_${POS}.bam bam/variant_${CHROM}_${POS}.mates.bam "${CHROM}" "${POS}" "${POS}" 
    samtools sort -@ 10 -o bam/variant_${CHROM}_${POS}.mates.bam bam/variant_${CHROM}_${POS}.mates.bam
    samtools merge -f -c -o bam/variant_${CHROM}_${POS}.bam.t  bam/variant_${CHROM}_${POS}.bam bam/variant_${CHROM}_${POS}.mates.bam 
    mv bam/variant_${CHROM}_${POS}.bam.t bam/variant_${CHROM}_${POS}.bam 
    samtools index bam/variant_${CHROM}_${POS}.bam 
    #rm bam/variant_${CHROM}_${POS}.mates.bam 

    
done < <(bcftools query -f '%CHROM\t%POS\n' $VCF)

echo "Done creating per-variant BAM files"
