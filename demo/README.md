# Demo dataset

A small synthetic dataset for trying Allelescope and for following the user
manual. It is the data shown in the manual's screenshots.

Everything here is simulated: the reference is random sequence and the reads
are 150 bp pairs written directly against it, so nothing has biological
meaning. The whole dataset is about 250 KB, and a browser instance builds from
it in a few seconds.

## Files

| File | Contents |
| --- | --- |
| `demo_ref.fa`, `.fai` | Reference with three sequences: `demo_seg1` (6,000 bp), `demo_seg2` (3,000 bp), `demo_seg3` (1,500 bp) |
| `demo_reads.bam`, `.bai` | 4,900 coordinate-sorted paired-end reads, with `MD` and `NM` tags |
| `demo_calls.vcf.gz`, `.tbi` | Five called variants on `demo_seg1` |
| `simulate_demo_data.py` | The script that generated the data |

## Try it

From the repository root, with Allelescope and samtools installed (see the
[main README](../allelescope/README.md#installation)):

```bash
build_variant_browser \
  -b demo/demo_reads.bam \
  -r demo/demo_ref.fa \
  -v demo/demo_calls.vcf.gz \
  -o Browsers/demo \
  -t 2

start_server --data-dir Browsers --port 8000
```

Then open <http://localhost:8000/demo/>. Build into a folder outside `demo/`,
as above, so the generated files are not mixed with the dataset.

## What to look for

All positions are on `demo_seg1`. The ones marked *not in VCF* can be opened
with **Jump to position**.

| Position | In VCF | What was simulated | Try |
| --- | --- | --- | --- |
| 1200 | yes | SNP C>A in nearly all reads | A clear homozygous-looking call |
| 2450 | yes | SNP A>C in about 30% of read pairs | Exclude **Duplicates** and **mates** and compare the allele balance |
| 2800 | not in VCF | Two low-frequency alternative bases (C about 9%, G about 7%) | Add the position; use the allele panel to filter to the minor alleles |
| 3300 | yes | 7-base deletion in about 40% of read pairs | Allele panel compares at 3301; filter to `DEL 7` |
| 3700 | not in VCF | Breakpoint: split reads continue at `demo_seg2:801` | Include **split reads**, then open `demo_seg2:700-1100` |
| 4100 | yes | 3-base insertion (`GAT`) in about 25% of read pairs | Allele panel compares at 4101; filter to `INS GAT` |
| 5000 | yes | SNP T>A carried only by reverse-strand reads | Switch **Strand** between Plus and Minus: a strand-bias artefact |
| 5150 to 5800 | not in VCF | A few pairs with both mates on the forward strand | Add 5300, tick **View as Pairs**, include orientation **F1F2** |

Across the whole dataset about 8% of read pairs are duplicated (the build step
marks them) and about 5% of reads have a soft-clipped end. `demo_seg2` and
`demo_seg3` have no variants; `demo_seg2` receives the other half of the split
reads.

## Regenerating the data

The script uses a fixed random seed, so it reproduces these files exactly.
It needs `pysam`, `samtools`, `bgzip` and `tabix`:

```bash
python demo/simulate_demo_data.py /tmp/demo_out
samtools calmd -b /tmp/demo_out/demo_reads.bam /tmp/demo_out/demo_ref.fa > demo/demo_reads.bam
samtools index demo/demo_reads.bam
cp /tmp/demo_out/demo_ref.fa demo/ && samtools faidx demo/demo_ref.fa
bgzip -c /tmp/demo_out/demo_calls.vcf > demo/demo_calls.vcf.gz
tabix -f -p vcf demo/demo_calls.vcf.gz
```
