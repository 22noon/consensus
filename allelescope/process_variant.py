from .extract_mates import extract_mates, _write_indel_stats
from .baq_filter import baq_filter
from pathlib import Path
import pysam


def process_variant(
    serverpath,
    Chrom: str,
    Pos: str,
    Ref: str,
    source_bam: str = "alignment.markdup.bam",
    threads: int = 10,
    window: int = 2,
    primary_only: bool = True,
    add_indel_tags: bool = True,
) -> int:
    """
    Pure-Python replacement for processed_bam() / extract_bam.sh.
    Returns 0 on success, 1 if the output BAM is empty.
    """
    print(f"Processing BAM for {Chrom}:{Pos} ref={Ref}", flush=True)

    base     = Path(serverpath)
    bam_dir  = base / "bam"
    bam_dir.mkdir(parents=True, exist_ok=True)
    Ref      = base / f"{Ref}"

    out_bam    = bam_dir / f"variant_{Chrom}_{Pos}.bam"
    mates_bam  = bam_dir / f"variant_{Chrom}_{Pos}.mates.bam"
    merged_tmp = bam_dir / f"variant_{Chrom}_{Pos}.bam.t"
    calmd_bam  = bam_dir / f"variant_{Chrom}_{Pos}.calmd.bam"
    src_bam    = str(base / source_bam)

    try:
        position = int(Pos)
    except (TypeError, ValueError):
        position = 0

    if position <= 0:
        # The tab placeholders use position 0.  They have no variant to
        # extract reads around, but downstream code still expects both the
        # source and filtered BAM paths to be valid and indexed.
        with pysam.AlignmentFile(src_bam, "rb") as source:
            for bam_path in (out_bam, calmd_bam):
                print(f"Creating empty BAM {bam_path} for {Chrom}:{Pos}", flush=True)
                with pysam.AlignmentFile(str(bam_path), "wb", header=source.header):
                    pass
                pysam.index(str(bam_path))
        return 1

    # Already processed — skip
    if out_bam.exists():
        return 0

    # --- Step 1: extract overlapping reads + mates ---
    allele_counts = extract_mates(
        src_bam   = src_bam,
        out_bam   = str(out_bam),
        mates_bam = str(mates_bam),
        chrom     = Chrom,
        start     = int(Pos),
        end       = int(Pos),
    )
    _write_indel_stats(allele_counts, str(out_bam))

    if not mates_bam.exists() or mates_bam.stat().st_size == 0:
        print(f"No mates found for {Chrom}:{Pos}, creating empty BAM")
    else:
        print(f"Mates found for {Chrom}:{Pos}, sorting and merging...")

        # --- Step 2: sort mates ---
        pysam.sort("-@", str(threads), "-o", str(mates_bam), str(mates_bam))

        # --- Step 3: merge spanning reads + mates ---
        pysam.merge("-f", "-c", "-o", str(merged_tmp), str(out_bam), str(mates_bam))
        mates_bam.unlink()

        # --- Step 4: calmd (recalculate MD/NM tags against reference) ---
        #calmd_bytes = pysam.calmd("-b", "-r", str(merged_tmp), Ref)
        calmd_bytes = pysam.calmd("-b", "-r", str(merged_tmp), str(Ref))

        calmd_bam.write_bytes(
            calmd_bytes if isinstance(calmd_bytes, bytes) else calmd_bytes.encode()
        )
        pysam.index(str(calmd_bam))
        merged_tmp.unlink()

        # --- Step 5: BAQ filter + annotation ---
        baq_filter(
            chrom         = Chrom,
            pos           = int(Pos),
            in_bam        = str(calmd_bam),
            out_bam       = str(out_bam),
            window        = window,
            primary_only  = primary_only,
            add_indel_tags = add_indel_tags,
        )

    # --- Step 6: index and count ---
    pysam.index(str(out_bam))

    with pysam.AlignmentFile(str(out_bam), "rb") as bam:
        count = bam.count()

    return 0 if count > 0 else 1
