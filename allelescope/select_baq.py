#!/usr/bin/env python3
"""
List read names that cross a locus and PASS BAQ-adjusted baseQ >= threshold.

Defaults approximate `bcftools mpileup -Q 13` (BAQ on, no -A):
- stepper="samtools"
- ignore_orphans=True
- collapse overlapping mates (ignore_overlaps=True)
- filters: unmapped, secondary, supplementary, QC-fail, duplicates are excluded

Usage:
  python reads_crossing_site_BAQpass.py \
      -b input.bam \
      -f reference.fa \
      -c MN632612.1 \
      -p 4572 \
      -q 13 \
      -o reads_baqpass.txt \
      --print-details
"""

import argparse
import sys
import pysam

def parse_args():
    ap = argparse.ArgumentParser(
        description="Read names crossing a locus with BAQ-adjusted baseQ >= threshold"
    )
    ap.add_argument("-b", "--bam", required=True, help="Input BAM/CRAM (indexed)")
    ap.add_argument("-f", "--fasta", required=True, help="Reference FASTA (indexed .fai)")
    ap.add_argument("-c", "--chrom", required=True, help="Chromosome/contig (e.g., MN632612.1)")
    ap.add_argument("-p", "--pos", type=int, required=True, help="1-based genomic position (e.g., 4572)")
    ap.add_argument("-q", "--min-baq", type=int, default=13,
                    help="Minimum BAQ-adjusted base quality to keep (default: 13)")
    ap.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    ap.add_argument("--max-depth", type=int, default=800000,
                    help="Max per-position depth to avoid downsampling (default: 800000)")
    ap.add_argument("--stepper", choices=["samtools", "all"], default="samtools",
                    help="Pileup stepper (default: samtools). Use 'all' for permissive behavior.")
    ap.add_argument("--include-orphans", action="store_true",
                    help="Include orphaned/discordant pairs (like bcftools -A). Default: off.")
    ap.add_argument("--no-overlap-collapse", action="store_true",
                    help="Do NOT collapse overlapping mates (default collapses them).")
    ap.add_argument("--print-details", action="store_true",
                    help="Print details (chrom, pos, name, BAQ, base, strand, MAPQ) instead of names only.")
    ap.add_argument("--require-proper-pair", action="store_true",
                    help="Keep only properly paired reads (optional).")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="Minimum mapping quality to keep a read (default: 0)")
    return ap.parse_args()

def main():
    args = parse_args()

    # Open reference and BAM
    try:
        ref = pysam.FastaFile(args.fasta)
    except Exception as e:
        sys.stderr.write(f"[error] Failed to open FASTA: {e}\n")
        sys.exit(1)

    try:
        bam = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        sys.stderr.write(f"[error] Failed to open BAM/CRAM: {e}\n")
        sys.exit(1)

    region = f"{args.chrom}:{args.pos}-{args.pos}"
    out = sys.stdout if args.output == "-" else open(args.output, "w")

    if args.print_details:
        print("#chrom\tpos\tread_name\tbaq\tbase_read_orient\tbase_ref_orient\tstrand\tmapq", file=out)

    # Build pileup kwargs safely; some pysam versions vary in supported kwargs
    pileup_kwargs = dict(
        region=region,
        stepper=args.stepper,           # "samtools" ≈ bcftools defaults; "all" is permissive
        fastafile=ref,
        compute_baq=True,               # BAQ ON
        min_base_quality=args.min_baq,  # filter on BAQ-adjusted base qualities
        max_depth=args.max_depth,
        truncate=True
    )

    # These kwargs may not exist in some versions; set if possible
    # Collapse overlapping mates unless explicitly disabled
    try:
        pileup_kwargs["ignore_overlaps"] = (not args.no_overlap_collapse)
    except Exception:
        pass
    # Ignore orphans by default (like no -A); include if requested
    try:
        pileup_kwargs["ignore_orphans"] = (not args.include_orphans)
    except Exception:
        pass

    pileup_iter = bam.pileup(**pileup_kwargs)

    # Deduplicate by read name so each template appears once
    seen = set()
    found_any_column = False

    # Complement table for optional reporting of base in reference-forward orientation
    comp_tab = str.maketrans("ACGTNacgtn", "TGCANtgcan")

    for col in pileup_iter:
        if col.reference_name != args.chrom or (col.pos + 1) != args.pos:
            continue
        found_any_column = True

        for pr in col.pileups:
            aln = pr.alignment

            # Manual per-read filters to mimic bcftools/samtools defaults
            if aln.is_unmapped:
                continue
            if aln.is_secondary:
                continue
            if aln.is_supplementary:
                continue
            if aln.is_qcfail:
                continue
            if aln.is_duplicate:
                continue
            if aln.mapping_quality < args.min_mapq:
                continue
            if args.require_proper_pair and not aln.is_proper_pair:
                continue

            # Skip deletions and reference skips at this position
            if pr.is_del or pr.is_refskip:
                continue

            qpos = pr.query_position
            if qpos is None:
                continue

            # BAQ-adjusted qualities are applied internally; min_base_quality already filtered most,
            # but we re-check to be explicit/robust
            baq_q = aln.query_qualities[qpos]
            if baq_q is None or baq_q < args.min_baq:
                continue

            rn = aln.query_name
            if rn in seen:
                continue
            seen.add(rn)

            if args.print_details:
                base_read_orient = aln.query_sequence[qpos]
                # Convert to reference-forward base for reporting (useful to inspect allele)
                base_ref_orient = base_read_orient.translate(comp_tab) if aln.is_reverse else base_read_orient
                strand = "-" if aln.is_reverse else "+"
                print(
                    f"{args.chrom}\t{args.pos}\t{rn}\t{baq_q}\t"
                    f"{base_read_orient}\t{base_ref_orient}\t{strand}\t{aln.mapping_quality}",
                    file=out
                )
            else:
                print(rn, file=out)

    if not found_any_column:
        sys.stderr.write(
            f"[warn] No pileup produced for {region}. "
            "Check contig name, indexing (.bai/.crai), and that reads overlap the site.\n"
        )

    if out is not sys.stdout:
        out.close()
    bam.close()
    ref.close()

if __name__ == "__main__":
    main()
