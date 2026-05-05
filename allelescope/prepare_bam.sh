#!/bin/bash
BAM=$1
REF=$2
VARIANT_LIST=$3
mkdir -p bam
samtools collate -u -@ 10 -O $BAM tmp_prefix | \
samtools fixmate -@ 10 -m -u - - | \
samtools sort -@ 10 -u -O bam - | \
samtools markdup -@ 10 - alignment.markdup.bam
samtools index alignment.markdup.bam
samtools faidx $REF
tabix -p vcf $VARIANT_LIST 

