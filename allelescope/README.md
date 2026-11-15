# Allelescope Variant Browser

Allelescope is an interactive, browser-based variant viewer. Each **browser instance** is a self-contained directory of processed alignment, reference, and variant files alongside a  viewer script. A lightweight file server serves one or more instances from a shared **data directory**, and users navigate to each instance in their web browser. The speciality of Allelescope is it allows deeper analysis of reads supporting a variant, offering advanced filtering capabilities that cannot be replicated in browsers like IGV for technical reasons.

Allelescope requires minimally a sorted bam file and a reference. If a VCF file is available variant locations in it can automatically be added to the browser. 

For instructions on using the generated viewer, see the [Allelescope user guide](USER_GUIDE.md).

To try Allelescope without your own data, use the small synthetic dataset in [`demo/`](../demo/README.md); it builds a browser instance in a few seconds.

---
## Installation
Allelescope can be obtained from PyPI. You can simply type 
```bash
pip install allelescope
```
and the package will be installed. The package requires Python 3.12 or above.

## Workflow overview

```mermaid

graph TD
    subgraph "Raw inputs"
        A[sample.bam]
        B[reference.fa]
        C[calls.vcf.gz]
    end
    A --> D
    B --> D
    C --> D
    D[build_variant_browser] --> E[Browsers/BTV/]
    E --> F[alignment.markdup.bam]
    E --> G[reference.fa]
    E --> H[calls.vcf.gz]
    E --> I[interactive_variants.html]
    J[start_server --data-dir Browsers --port 8000] --> K[localhost:8000/BTV]
    J --> L[localhost:8000/RPV]
```

---

## VCF preparation with bcftools

For best results, variant files should be prepared with `bcftools` before being passed to the build script. This ensures the tags that Allelescope relies on are present and correctly formatted.

```bash
# Normalise, left-align indels, and fill in missing tags
bcftools norm -f reference.fa -m -any raw_calls.vcf.gz \
  | bcftools +fill-tags -- -t all \
  | bcftools view -o calls.norm.vcf.gz -O z

bcftools index -t calls.norm.vcf.gz
```

This will create the normalized variant file `calls.norm.vcf.gz`. This can be passed in the step below to the build script with `-v` option. (This vcf preparation was tested with bcftools version 1.19.)

---

## Building a Browser For a Dataset 

`build_variant_browser` handles the data preparation step: it takes raw input files, processes and indexes them, and produces the files needed for the viewer — all placed into the instance output directory. The script is invoked with the following syntax:


```
build_variant_browser -b BAM -r FASTA [-v VCF] [-g GENOME_NAME] [-o OUT_DIR] [-t THREADS]
                      [--bam-processing auto|full|skip] [--link]
```
The detailed description of flags are shown below.

| Flag | Required | Default | Description |
|---|---|---|---|
| `-b` | Yes | — | Input BAM file |
| `-r` | Yes | — | Reference FASTA file |
| `-v` | No | — | Variant file (VCF or bgzipped VCF.gz, must be sorted) |
| `-g` | No | auto-inferred | Genome name passed to `generate_viewer.py` |
| `-o` | No | `./out_dir` | Output directory (created if it does not exist) |
| `-t` | No | `10` | Number of threads for samtools |
| `--bam-processing` | No | `auto` | `auto`: run only the BAM processing steps the file still needs (see below); `full`: always run every step; `skip`: use a coordinate-sorted BAM as it is |
| `--link` | No | off | When the BAM needs no processing, hard-link it into the output directory instead of copying it. Only works on the same file system; otherwise the BAM is copied |
| `-h` | — | — | Print help and exit |


If `-g` is not supplied, the genome name is resolved in this order:

1. First `CHROM` value in the variant file (first non-header data row)
2. First reference sequence name from the BAM `@SQ` header

If neither can be determined, the script exits with an error and asks you to supply `-g` explicitly.


### Dependencies
The following tools are assumed to be installed in the current software environment for the pre-processing work of `build_variant_browser`.

| Tool | Required | Notes |
|---|---|---|
| `samtools` | Always | collate, fixmate, sort, markdup, index, faidx |
| `tabix` | Only with `-v` | indexes the VCF |
| `bgzip` | Only with `-v` | used if the supplied VCF is not already gzipped |

All three are available through [htslib](https://www.htslib.org/download/) or via conda/micromamba. The package was developed using samtools 1.19.2 and htslib 1.19, so it is advised to be cautious using very old versions of samtools that might lack some features. The following can instructs how to build an environment to run Allelescope.

```bash
# Using conda
conda create -n allelescope -c conda-forge -c bioconda samtools=1.19.2 htslib=1.19 
conda activate allelescope

# Using micromamba
micromamba create -n allelescope -c conda-forge -c bioconda samtools=1.19.2 htslib=1.19 
micromamba activate allelescope
```

if you want a compatible 'bcftools' installed in the same environment run:

```bash
conda install bcftools=1.19
```
or the equivalent command for the environment.

If you already have a suitable active environment, you can also install directly:

```bash
conda install -c bioconda samtools htslib
```



---

## Examples

**Minimal — BAM and reference only:**
```bash
build_variant_browser \
  -b sample.bam \
  -r reference.fa
```

**With variants:**
```bash
build_variant_browser \
  -b sample.bam \
  -r reference.fa \
  -v calls.norm.vcf.gz
```

**Full options:**
```bash
build_variant_browser.sh \
  -b sample.bam \
  -r reference.fa \
  -v calls.norm.vcf.gz \
  -g D388-WT_OR813926 \
  -o /data/browser_output \
  -t 16
```

---

## What the script does

1. **Validates** that all supplied input files exist and required tools are available.
2. **Resolves the genome name** from the VCF or BAM header if not provided.
3. **Processes the BAM** and writes `alignment.markdup.bam` to the output directory. A raw BAM goes through a `collate → fixmate → sort → markdup` pipeline. With `--bam-processing auto` (the default), the script first reads the BAM header and the first 20,000 records, and skips the steps that have already been done:

   | The BAM | Steps run |
   |---|---|
   | is coordinate-sorted and an `@PG` header line shows a duplicate-marking program (samtools markdup, Picard or GATK MarkDuplicates, sambamba, biobambam, samblaster) | none: the BAM is copied (or hard-linked with `--link`) and its index is reused if it is up to date |
   | has duplicates marked but is not coordinate-sorted | `sort` |
   | is coordinate-sorted and every sampled pair has the `MC` and `ms` tags that `samtools fixmate -m` adds | `markdup` |
   | anything else, including single-end BAMs without duplicate marking | `collate → fixmate → sort → markdup` |

   `collate` leaves no trace of its own: it only groups reads by name so that `fixmate` can see both mates, so it is run only together with `fixmate`. Duplicate flags without a recognised program in the header are not trusted, because they may come from an earlier, unrelated run; use `--bam-processing skip` to keep them. Duplicates marked by another tool are kept as that tool marked them.
4. **Indexes the BAM** (`.bai`).
5. **Copies the reference FASTA** to the output directory and indexes it with `samtools faidx` (`.fai`).
6. **Copies and indexes the VCF** (if provided). If the file is not bgzipped, it is bgzipped automatically before indexing with `tabix` (`.tbi`).
7. **Invokes `generate_viewer.py`** if it is found in the same directory as the script. If not found, the equivalent command is printed so you can run it manually.

---

## Output directory layout

```
out_dir/
├── alignment.markdup.bam       # Duplicate-marked BAM
├── alignment.markdup.bam.bai   # BAM index
├── reference.fa                # Copy of the reference FASTA
├── reference.fa.fai            # FASTA index
├── calls.norm.vcf.gz           # Copy of the variant file (bgzipped)
├── calls.norm.vcf.gz.tbi       # Tabix index
└── interactive_variants.html   # Generated viewer (entry point)
```

> The variant files are only present if `-v` was supplied.

---

## Starting the file server

The file server takes a **data directory** (the root from which instances are served) and an optional **port**. It does not need to be restarted when new instances are added.

```bash
start_server --data-dir /path/to/Browsers --port 8000
```

| Parameter | Default | Description |
|---|---|---|
| `--data-dir` | current directory | Root directory containing one or more browser instances |
| `--port` | `8000` | Port to listen on |

---

## Accessing browser instances

Once the server is running, navigate to:

```
http://localhost:<port>/<relative-path>
```

where `<relative-path>` is the path from the data directory to the specific instance directory.

## Accessing from a local machine via SSH tunnel

If the file server is running on a remote host, you can forward the server port to your local machine with SSH and then open the browser locally.

On your local machine, run:

```bash
ssh -L <port>:localhost:<port> user@remote.host
```

Then open the URL in Chrome, Edge, or Firefox:

```
http://localhost:<port>/<relative-path>
```

For example, if the remote server uses a port 9000, adjust the tunnel like:

```bash
ssh -L 9000:localhost:8000 user@remote.host
```

and then open:

```
http://localhost:9000/<relative-path>
```

### Multi-instance example

Suppose you have two samples prepared under a shared `Browsers/` directory:

```
Browsers/
├── BTV/
│   ├── alignment.markdup.bam
│   ├── reference.fa
│   └── interactive_variants.html
└── RPV/
    ├── alignment.markdup.bam
    ├── reference.fa
    └── interactive_variants.html
```

Start the server with `Browsers/` as the data directory:

```bash
start_server --data-dir Browsers --port 8000
```

Then open each instance in your browser:

| Sample | URL |
|---|---|
| BTV | http://localhost:8000/BTV |
| RPV | http://localhost:8000/RPV |

To build the two instances in the first place:

```bash
build_variant_browser.sh -b BTV.bam -r BTV_ref.fa -v BTV_calls.vcf.gz -o Browsers/BTV
build_variant_browser.sh -b RPV.bam -r RPV_ref.fa -v RPV_calls.vcf.gz -o Browsers/RPV
```

---

## Notes

- The script exits immediately on any error (`set -euo pipefail`). If a step fails, check the last printed stage heading to identify where it stopped.
- A plain, unsorted VCF will cause `tabix` to fail. Sort your VCF first with `bcftools sort` or `vcf-sort` before passing it to this script.
