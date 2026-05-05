#!/usr/bin/env python3
import argparse
import statistics
import pysam

def _baq_at_qpos(qual_list, bq_str, qpos):
    """Compute BAQ-adjusted Phred at query index qpos using QUAL and BQ.
    BAQ_i = QUAL_i - (ord(BQ_i) - 64), floored at 0.
    """
    off = ord(bq_str[qpos]) - 64
    return max(0, int(qual_list[qpos]) - off)

def make_key(aln):
    """Key that uniquely identifies an alignment record across passes."""
    return (
        aln.query_name,
        aln.is_secondary,
        aln.is_supplementary,
        aln.reference_id,
        aln.reference_start,
        aln.cigarstring,
    )

def annotate_baq_with_neighbors_cli():
    ap = argparse.ArgumentParser(
        description="Annotate reads at a single site with BAQ-derived tags. "
                    "XB: BAQ at center or 0; XS: suspect by flanking BAQ; "
                    "XD/XR: deletion/refskip at center; XI/XL: indel lengths at center."
    )
    ap.add_argument("--chrom", required=True, help="Chromosome (e.g., chr7)")
    ap.add_argument("--pos", type=int, required=True, help="1-based position")
    ap.add_argument("--in-bam", required=True, help="Input BAM (with BQ tag from samtools calmd -r)")
    ap.add_argument("--out-bam", required=True, help="Output BAM with tags")
    # Fine-tuning knobs (defaults mirror earlier discussion)
    ap.add_argument("--window", type=int, default=2, help="Flank radius in bp (default: 2 -> ±2)")
    ap.add_argument("--primary-only", action="store_true", help="Ignore secondary/supplementary at pileup time")
    ap.add_argument("--min-mapq", type=int, default=None, help="Min MAPQ at pileup (mpileup -q analogue)")
    ap.add_argument("--min-baseq", type=int, default=None, help="Min baseQ at pileup (applies where base exists)")
    ap.add_argument("--ignore-overlaps", action="store_true", help="Avoid double-counting overlapping mates")
    ap.add_argument("--require-proper-pair", action="store_true", help="Ignore orphans at pileup")
    ap.add_argument("--max-depth", type=int, default=None, help="Max depth (mpileup -d analogue)")
    ap.add_argument("--neighbor-median-thr", type=int, default=5, help="XS: median(flank BAQ) < thr => suspect")
    ap.add_argument("--neighbor-low-thr", type=int, default=5, help="XS: flank BAQ < thr counts as 'low'")
    ap.add_argument("--neighbor-low-count-min", type=int, default=2, help="XS: min #low flanks to flag suspect")
    ap.add_argument("--add-indel-tags", action="store_true", help="Also emit XI/XL at center")
    args = ap.parse_args()

    chrom = args.chrom
    pos_1b = args.pos

    # Window positions (1-based)
    positions_1b = [p for p in range(pos_1b - args.window, pos_1b + args.window + 1) if p > 0]

    # Per-alignment info aggregated across positions
    info = {}
    def ensure_slot(key):
        if key not in info:
            info[key] = {
                "xb_center": 0,   # BAQ at center
                "flank_baqs": [], # BAQs at flanks
                "xd": 0,          # deletion spans center
                "xr": 0,          # refskip spans center
                "xi": 0,          # insertion length at center (>0)
                "xl": 0,          # deletion length at center (<0)
                "saw_center": False,
            }

    # Pileup kwargs approximating mpileup defaults
    pileup_kwargs = dict(stepper="samtools", truncate=True)
    if args.min_mapq is not None:
        pileup_kwargs["min_mapping_quality"] = int(args.min_mapq)
    if args.min_baseq is not None:
        pileup_kwargs["min_base_quality"] = int(args.min_baseq)
    if args.ignore_overlaps:
        pileup_kwargs["ignore_overlaps"] = True
    if args.require_proper_pair:
        pileup_kwargs["ignore_orphans"] = True
    if args.max_depth is not None:
        pileup_kwargs["max_depth"] = int(args.max_depth)

    # PASS 1: scan center and flanks
    with pysam.AlignmentFile(args.in_bam, "rb") as bam_in:
        for p1 in positions_1b:
            start0, end0 = p1 - 1, p1
            for col in bam_in.pileup(reference=chrom, start=start0, end=end0, **pileup_kwargs):
                if col.reference_name != chrom or (col.pos + 1) != p1:
                    continue

                is_center = (p1 == pos_1b)

                for pr in col.pileups:
                    aln = pr.alignment
                    if args.primary_only and (aln.is_secondary or aln.is_supplementary):
                        continue

                    key = make_key(aln)
                    ensure_slot(key)

                    # Refskip at this position
                    if pr.is_refskip:
                        if is_center:
                            info[key]["xr"] = 1
                            info[key]["saw_center"] = True
                        continue

                    # Deletion spanning this position
                    if pr.is_del:
                        if is_center:
                            info[key]["xd"] = 1
                            info[key]["saw_center"] = True
                            if args.add_indel_tags and pr.indel < 0:
                                info[key]["xl"] = int(pr.indel)
                        continue

                    # From here: aligned base exists
                    qpos = pr.query_position
                    if qpos is None:
                        continue

                    if is_center and args.add_indel_tags and pr.indel > 0:
                        info[key]["xi"] = int(pr.indel)

                    baq_here = None
                    if aln.query_qualities is not None and aln.has_tag("BQ"):
                        bq_str = aln.get_tag("BQ")
                        if qpos < len(aln.query_qualities) and qpos < len(bq_str):
                            baq_here = _baq_at_qpos(aln.query_qualities, bq_str, qpos)

                    if is_center:
                        info[key]["saw_center"] = True
                        if baq_here is not None:
                            info[key]["xb_center"] = int(baq_here)
                    else:
                        if baq_here is not None:
                            info[key]["flank_baqs"].append(int(baq_here))

        # Reset for PASS 2
        bam_in.reset()

        with pysam.AlignmentFile(args.out_bam, "wb", template=bam_in) as bam_out:
            for aln in bam_in.fetch(until_eof=True):
                key = make_key(aln)
                slot = info.get(key)

                xb = 0; xd = 0; xr = 0; xi = 0; xl = 0; xs = 0
                if slot is not None:
                    xb, xd, xr, xi, xl = (
                        int(slot["xb_center"]),
                        int(slot["xd"]),
                        int(slot["xr"]),
                        int(slot["xi"]),
                        int(slot["xl"]),
                    )

                    # XS decision:
                    # We flag suspect if center lacks a confident base (xb==0)
                    # AND flanking BAQ distribution is poor by thresholds.
                    flanks = slot["flank_baqs"]
                    if xb == 0 and flanks:
                        med = statistics.median(flanks)
                        low_ct = sum(1 for v in flanks if v < args.neighbor_low_thr)
                        xb = med # imputed BAQ at center for reference; XS=1 if suspect by neighbors 
                        if (med < args.neighbor_median_thr) or (low_ct >= args.neighbor_low_count_min):
                            xs = 1

                # Always write tags; mates that don't cross the site keep XB=0.
                aln.set_tag("XB", int(xb), value_type="i")
                aln.set_tag("XS", int(xs), value_type="i")
                #aln.set_tag("XD", int(xd), value_type="i")
                aln.set_tag("XR", int(xr), value_type="i")
                #if args.add_indel_tags:
                #    aln.set_tag("XI", int(xi), value_type="i")
                #    aln.set_tag("XL", int(xl), value_type="i")

                bam_out.write(aln)

if __name__ == "__main__":
    annotate_baq_with_neighbors_cli()