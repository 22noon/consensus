from .extract_mates import extract_mates, _write_indel_stats, _write_evidence
from pathlib import Path
import hashlib
import json
import re
import threading
import pysam

# Bump whenever the tags written to the per-site BAMs change meaning, so that
# files produced by an older version are rebuilt instead of being reused.
#   2: XO "proper" covers both F1R2 and R1F2 proper pairs
#   3: XA reference/mismatch decided against the FASTA rather than MD tags;
#      BQ (BAQ) only on reads that cross the site
PROCESSING_VERSION = "5"
_VERSION_FILE = ".allelescope_version"
_checked_dirs = set()
_check_lock = threading.Lock()


_UNSAFE_NAME_CHARS = re.compile(r"[^A-Za-z0-9._+-]")


def site_sample(value):
    """Parse a requested sample size: a positive int, or None for all reads."""
    try:
        n = int(str(value).strip())
    except (TypeError, ValueError):
        return None
    return n if n > 0 else None


def site_stem(chrom, pos, sample=None) -> str:
    """
    File-name stem for the per-site BAMs of chrom:pos, e.g. "variant_chr1_1200",
    or "variant_chr1_1200.s2000" for the version subsampled to 2,000 read pairs.

    Sequence names may contain characters that are not allowed in file names
    on Windows (GRCh38 has HLA contigs such as "HLA-A*01:01:01:01") or that
    would act as path separators. Those characters are replaced by "_" and a
    short hash of the original name is appended, so different names cannot
    collide. Names made only of letters, digits and "._+-" are used as they
    are, which keeps the file names of earlier versions.
    """
    chrom = str(chrom)
    safe = _UNSAFE_NAME_CHARS.sub("_", chrom)
    if safe != chrom or safe in ("", ".", ".."):
        digest = hashlib.sha1(chrom.encode("utf-8")).hexdigest()[:8]
        safe = f"{safe}-{digest}"
    stem = f"variant_{safe}_{int(pos)}"
    sample = site_sample(sample)
    return f"{stem}.s{sample}" if sample else stem


def site_info_path(serverpath, chrom, pos, sample=None) -> Path:
    return Path(serverpath) / "bam" / f"{site_stem(chrom, pos, sample)}.site.json"


def read_site_info(serverpath, chrom, pos, sample=None) -> dict:
    """Read-pair counts recorded when the site was prepared ({} if unknown)."""
    try:
        return json.loads(site_info_path(serverpath, chrom, pos, sample).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _calmd(in_bam, ref, out_bam, baq):
    """samtools calmd [-r] in_bam ref > out_bam (BAM output)."""
    args = ["-b", "-r"] if baq else ["-b"]
    data = pysam.calmd(*args, str(in_bam), str(ref))
    Path(out_bam).write_bytes(data if isinstance(data, bytes) else data.encode())


def ensure_current_bam_dir(serverpath) -> Path:
    """
    Make sure <instance>/bam only holds per-site BAMs written by this version.

    The per-site BAMs are a cache: they are rebuilt on demand from
    alignment.markdup.bam. If the directory was filled by an older version
    (no marker, or a different one) its contents are deleted so every site is
    re-extracted with the current tags. Checked once per directory per process,
    and again if the directory has been deleted since.
    """
    bam_dir = Path(serverpath) / "bam"
    with _check_lock:
        # Checked already, unless the folder has been deleted since (cache cleared).
        if bam_dir in _checked_dirs and bam_dir.is_dir():
            return bam_dir
        bam_dir.mkdir(parents=True, exist_ok=True)
        marker = bam_dir / _VERSION_FILE
        current = marker.read_text(encoding="utf-8").strip() if marker.exists() else None
        if current != PROCESSING_VERSION:
            stale = [p for p in bam_dir.iterdir() if p.is_file() and p.name != _VERSION_FILE]
            if stale:
                print(f"Per-site BAMs in {bam_dir} were written by an older version "
                      f"({current or 'unversioned'}); removing {len(stale)} file(s) so they are rebuilt.",
                      flush=True)
                for p in stale:
                    p.unlink()
            marker.write_text(PROCESSING_VERSION + "\n", encoding="utf-8")
        _checked_dirs.add(bam_dir)
    return bam_dir


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
    sample=None,
) -> int:
    """
    Pure-Python replacement for processed_bam() / extract_bam.sh.
    Returns 0 on success, 1 if the output BAM is empty.
    """
    print(f"Processing BAM for {Chrom}:{Pos} ref={Ref}", flush=True)

    base     = Path(serverpath)
    bam_dir  = ensure_current_bam_dir(base)
    Ref      = base / f"{Ref}"

    sample     = site_sample(sample)
    stem       = site_stem(Chrom, Pos, sample)
    out_bam    = bam_dir / f"{stem}.bam"
    mates_bam  = bam_dir / f"{stem}.mates.bam"
    merged_tmp = bam_dir / f"{stem}.bam.t"
    calmd_bam  = bam_dir / f"{stem}.calmd.bam"
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
    # The served file is the last one written, so it marks a finished site.
    if calmd_bam.exists() and out_bam.exists():
        return 0

    sample_info = {}
    allele_stats = {}
    allele_counts = extract_mates(
        src_bam   = src_bam,
        out_bam   = str(out_bam),
        mates_bam = str(mates_bam),
        chrom     = Chrom,
        start     = int(Pos),
        end       = int(Pos),
        ref_fasta = str(Ref),
        max_pairs = sample,
        sample_info = sample_info,
        allele_stats = allele_stats,
    )
    _write_indel_stats(allele_counts, str(out_bam), allele_stats)
    _write_evidence(allele_stats, str(out_bam))
    sample_info["sample"] = sample
    site_info_path(base, Chrom, Pos, sample).write_text(json.dumps(sample_info), encoding="utf-8")
    pysam.index(str(out_bam))

    # Recompute MD/NM against the reference for every record (samtools
    # calmd). Base alignment quality (BQ tag, calmd -r) is only computed for
    # the reads that cross the site: it is only ever evaluated there, and
    # BAQ is the costly part. Mates and split partners do not reach the site.
    site_part = bam_dir / f"{stem}.site.calmd.bam"
    _calmd(out_bam, Ref, site_part, baq=True)

    if mates_bam.exists() and mates_bam.stat().st_size > 0:
        print(f"Mates found for {Chrom}:{Pos}, sorting and merging...")
        mates_sorted = bam_dir / f"{stem}.mates.sorted.bam"
        mates_part   = bam_dir / f"{stem}.mates.calmd.bam"
        pysam.sort("-@", str(threads), "-o", str(mates_sorted), str(mates_bam))
        _calmd(mates_sorted, Ref, mates_part, baq=False)
        pysam.merge("-f", "-c", "-o", str(merged_tmp), str(site_part), str(mates_part))
        for tmp in (mates_bam, mates_sorted, mates_part, site_part):
            tmp.unlink(missing_ok=True)
    else:
        print(f"No mates found for {Chrom}:{Pos}")
        mates_bam.unlink(missing_ok=True)
        site_part.replace(merged_tmp)

    pysam.index(str(merged_tmp), str(calmd_bam) + ".bai")
    merged_tmp.replace(calmd_bam)

    with pysam.AlignmentFile(str(out_bam), "rb") as bam:
        count = bam.count()

    return 0 if count > 0 else 1
