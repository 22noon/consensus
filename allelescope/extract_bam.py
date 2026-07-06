import os
import sys
from pathlib import Path
import pysam
from baq_filter import baq_filter  
from extract_mates import extract_mates


def process_variant(
    chrom: str,
    pos: str,
    ref: str,
    bam_path: str,
    source_bam: str = "alignment.markdup.bam",
    threads: int = 10,
    window: int = 2,
    primary_only: bool = True,
    add_indel_tags: bool = True,
) -> int:
    """
    Process a BAM file for a specific genomic variant.

    Args:
        chrom: Chromosome name (e.g. "chr1")
        pos: Position as string (e.g. "12345")
        ref: Path to reference FASTA file
        bam_path: Base directory containing the BAM files
        source_bam: Source BAM filename (default: "alignment.markdup.bam")
        threads: Number of threads for samtools sort (default: 10)
        window: Window size for BAQ filter (default: 2)
        primary_only: Only include primary alignments (default: True)
        add_indel_tags: Add indel tags during annotation (default: True)

    Returns:
        0 on success (reads found), 1 if output BAM is empty
    """
    print(f"Processing {chrom}:{pos}")

    base = Path(bam_path)
    bam_dir = base / "bam"
    bam_dir.mkdir(parents=True, exist_ok=True)

    out_bam       = bam_dir / f"variant_{chrom}_{pos}.bam"
    mates_bam     = bam_dir / f"variant_{chrom}_{pos}.mates.bam"
    merged_tmp    = bam_dir / f"variant_{chrom}_{pos}.bam.t"
    calmd_bam     = bam_dir / f"variant_{chrom}_{pos}.calmd.bam"
    ref           = base / f"{ref}"

    # --- Skip if already done ---
    if out_bam.exists():
        return 0

    src_bam = str(base / source_bam)

    # --- Extract reads & mates (your existing extract_mates logic) ---
    extract_mates(src_bam, str(out_bam), str(mates_bam), chrom, pos, pos)

    if not mates_bam.exists() or mates_bam.stat().st_size == 0:
        print(f"No mates found for {chrom}:{pos}, creating empty BAM")
    else:
        print(f"Mates found for {chrom}:{pos}, sorting BAM")

        # samtools sort -@ 10 -o mates.bam mates.bam
        pysam.sort("-@", str(threads), "-o", str(mates_bam), str(mates_bam))

        # samtools merge -f -c -o merged.bam.t out.bam mates.bam
        pysam.merge("-f", "-c", "-o", str(merged_tmp), str(out_bam), str(mates_bam))
        mates_bam.unlink()

        # samtools calmd -b -r merged.bam.t ref.fa > calmd.bam
        calmd_output = pysam.calmd("-b", "-r", str(merged_tmp), ref)
        calmd_bam.write_bytes(calmd_output.encode()
                              if isinstance(calmd_output, str) else calmd_output)
        pysam.index(str(calmd_bam))
        merged_tmp.unlink()

        # Your baq_filter logic — call as a Python function instead of subprocess
        baq_filter(
            chrom=chrom,
            pos=int(pos),
            in_bam=str(calmd_bam),
            out_bam=str(out_bam),
            window=window,
            primary_only=primary_only,
            add_indel_tags=add_indel_tags,
        )

    # --- Index final BAM ---
    pysam.index(str(out_bam))

    # --- Count reads ---
    with pysam.AlignmentFile(str(out_bam), "rb") as bam:
        count = bam.count()

    return 0 if count > 0 else 1

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: process_variant.py CHROM POS REF BAM_PATH")
        sys.exit(1)

    exit_code = process_variant(
        chrom=sys.argv[1],
        pos=sys.argv[2],
        ref=sys.argv[3],
        bam_path=sys.argv[4],
    )
    sys.exit(exit_code)
