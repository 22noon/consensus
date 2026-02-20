CHROM=$1
POS=$2
echo "Processing ${CHROM}:${POS}"
if [ -f bam/variant_${CHROM}_${POS}.bam ]; then
    exit 0
fi
SCRIPT_PATH=$(dirname "$(readlink -f "${BASH_SOURCE[0]:-$0}")")
python $SCRIPT_PATH/extract_mates.py alignment.markdup.bam bam/variant_${CHROM}_${POS}.bam bam/variant_${CHROM}_${POS}.mates.bam "${CHROM}" "${POS}" "${POS}" 
if [ ! -s bam/variant_${CHROM}_${POS}.mates.bam ]; then
    echo "No mates found for ${CHROM}:${POS}, creating empty BAM"
else
    echo "Mates found for ${CHROM}:${POS}, sorting BAM"
    samtools sort -@ 10 -o bam/variant_${CHROM}_${POS}.mates.bam bam/variant_${CHROM}_${POS}.mates.bam
    samtools merge -f -c -o bam/variant_${CHROM}_${POS}.bam.t  bam/variant_${CHROM}_${POS}.bam bam/variant_${CHROM}_${POS}.mates.bam 
    rm bam/variant_${CHROM}_${POS}.mates.bam 
    mv bam/variant_${CHROM}_${POS}.bam.t bam/variant_${CHROM}_${POS}.bam 

samtools index bam/variant_${CHROM}_${POS}.bam 
count=$(samtools view -c bam/variant_${CHROM}_${POS}.bam )
if [ $count -eq 0 ]; then
    ERRORLEVEL=1
else
    ERRORLEVEL=0
fi
exit $ERRORLEVEL 