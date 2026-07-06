#!/usr/bin/env python3
"""Build an Allelescope/variant-browser instance from BAM, reference FASTA and optional VCF.

Python replacement for build_variant_browser.sh.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Optional


class BuildError(RuntimeError):
    pass


class Colours:
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    YELLOW = "\033[1;33m"
    CYAN = "\033[0;36m"
    BOLD = "\033[1m"
    RESET = "\033[0m"


def info(message: str) -> None:
    print(f"{Colours.CYAN}[INFO]{Colours.RESET}  {message}")


def success(message: str) -> None:
    print(f"{Colours.GREEN}[OK]{Colours.RESET}    {message}")


def warn(message: str) -> None:
    print(f"{Colours.YELLOW}[WARN]{Colours.RESET}  {message}")


def step(message: str) -> None:
    print(f"\n{Colours.BOLD}▶ {message}{Colours.RESET}")


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise BuildError(f"{label} not found: {path}")


def require_command(command: str) -> None:
    if shutil.which(command) is None:
        raise BuildError(f"'{command}' not found in PATH")


def run(cmd: Iterable[str], *, stdout=None) -> subprocess.CompletedProcess:
    cmd = [str(c) for c in cmd]
    info("Running: " + " ".join(cmd))
    try:
        return subprocess.run(cmd, check=True, stdout=stdout)
    except subprocess.CalledProcessError as exc:
        raise BuildError(f"Command failed with exit code {exc.returncode}: {' '.join(cmd)}") from exc


def first_vcf_chrom(vcf_path: Path) -> str:
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    mode = "rt"
    try:
        with opener(vcf_path, mode) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields and fields[0]:
                    return fields[0]
    except OSError:
        # Handles cases where the file ends with .gz but is not valid gzip.
        return ""
    return ""


def first_bam_reference(bam_path: Path) -> str:
    proc = subprocess.run(
        ["samtools", "view", "-H", str(bam_path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    for line in proc.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        for field in line.split("\t"):
            if field.startswith("SN:"):
                return field[3:]
    return ""


def resolve_genome_name(genome_name: str, variants: Optional[Path], bam: Path) -> str:
    if genome_name:
        info(f"Genome name provided: {Colours.BOLD}{genome_name}{Colours.RESET}")
        return genome_name

    inferred = ""
    if variants is not None:
        inferred = first_vcf_chrom(variants)

    if not inferred:
        if variants is not None:
            warn("No variants found in VCF; falling back to first reference in BAM")
        inferred = first_bam_reference(bam)

    if not inferred:
        raise BuildError("Could not determine genome name automatically. Please supply -g.")

    info(f"Genome name inferred: {Colours.BOLD}{inferred}{Colours.RESET}")
    return inferred


def process_bam(bam: Path, out_dir: Path, threads: int) -> Path:
    markdup_bam = out_dir / "alignment.markdup.bam"
    tmp_prefix = out_dir / "tmp_markdup"

    step("Processing BAM (collate → fixmate → sort → markdup)")
    info(f"Threads: {threads}")

    # Equivalent to:
    # samtools collate -u -@ THREADS -O BAM TMP_PREFIX \
    #   | samtools fixmate -@ THREADS -m -u - - \
    #   | samtools sort -@ THREADS -u -O bam - \
    #   | samtools markdup -@ THREADS - MARKDUP_BAM
    p1 = subprocess.Popen(
        ["samtools", "collate", "-u", "-@", str(threads), "-O", str(bam), str(tmp_prefix)],
        stdout=subprocess.PIPE,
    )
    p2 = subprocess.Popen(
        ["samtools", "fixmate", "-@", str(threads), "-m", "-u", "-", "-"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
    )
    assert p1.stdout is not None
    p1.stdout.close()
    p3 = subprocess.Popen(
        ["samtools", "sort", "-@", str(threads), "-u", "-O", "bam", "-"],
        stdin=p2.stdout,
        stdout=subprocess.PIPE,
    )
    assert p2.stdout is not None
    p2.stdout.close()
    p4 = subprocess.Popen(
        ["samtools", "markdup", "-@", str(threads), "-", str(markdup_bam)],
        stdin=p3.stdout,
    )
    assert p3.stdout is not None
    p3.stdout.close()

    exit_codes = [p.wait() for p in (p1, p2, p3, p4)]
    if any(code != 0 for code in exit_codes):
        raise BuildError(f"BAM processing pipeline failed; exit codes: {exit_codes}")

    success(f"markdup BAM written: {markdup_bam}")

    step("Indexing BAM")
    run(["samtools", "index", str(markdup_bam)])
    success(f"BAM index written: {markdup_bam}.bai")
    return markdup_bam


def copy_and_index_reference(ref: Path, out_dir: Path) -> Path:
    step("Copying and indexing reference FASTA")
    ref_dest = out_dir / ref.name
    shutil.copy2(ref, ref_dest)
    run(["samtools", "faidx", str(ref_dest)])
    success(f"FASTA index written: {ref_dest}.fai")
    return ref_dest


def copy_and_index_vcf(variants: Optional[Path], out_dir: Path) -> Optional[Path]:
    if variants is None:
        info("No variant file provided — skipping VCF step")
        return None

    step("Copying and indexing variant file")
    vcf_dest = out_dir / variants.name

    if variants.name.endswith(".gz"):
        shutil.copy2(variants, vcf_dest)
    else:
        warn(f"VCF is not gzipped — bgzipping to {vcf_dest}.gz")
        vcf_dest = Path(str(vcf_dest) + ".gz")
        with open(vcf_dest, "wb") as out_handle:
            run(["bgzip", "-c", str(variants)], stdout=out_handle)

    run(["tabix", "-p", "vcf", str(vcf_dest)])
    success(f"Variant index written: {vcf_dest}.tbi")
    return vcf_dest


def run_viewer_generator(
    genome_name: str,
    ref_dest: Path,
    markdup_bam: Path,
    vcf_dest: Optional[Path],
    out_dir: Path,
) -> None:
    step("Running generate_viewer.py")

    # Works both from a source checkout and after pip installation, provided
    # generate_viewer.py is included in the allelescope package.
    generator = Path(__file__).with_name("generate_viewer.py")
    if not generator.is_file():
        raise BuildError(f"generate_viewer.py not found next to this file: {generator}")

    cmd = [
        sys.executable,
        str(generator),
        "--genome-name",
        genome_name,
        "--reference-fasta",
        str(ref_dest),
        "--bam-path",
        str(markdup_bam),
        "--output-dir",
        str(out_dir),
    ]
    if vcf_dest is not None:
        cmd.extend(["--vcf-path", str(vcf_dest)])

    run(cmd)
    success(f"Viewer generated in: {out_dir}")


def build_variant_browser(args: argparse.Namespace) -> None:
    bam = Path(args.bam).expanduser().resolve()
    ref = Path(args.reference).expanduser().resolve()
    variants = Path(args.variants).expanduser().resolve() if args.variants else None
    out_dir = Path(args.out_dir).expanduser().resolve()

    step("Validating inputs")
    require_file(bam, "BAM file")
    require_file(ref, "FASTA file")
    if variants is not None:
        require_file(variants, "Variant file")

    require_command("samtools")
    if variants is not None:
        require_command("tabix")
        require_command("bgzip")
    success("All inputs present")

    step("Resolving genome name")
    genome_name = resolve_genome_name(args.genome_name, variants, bam)

    step(f"Preparing output directory: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    markdup_bam = process_bam(bam, out_dir, args.threads)
    ref_dest = copy_and_index_reference(ref, out_dir)
    vcf_dest = copy_and_index_vcf(variants, out_dir)
    run_viewer_generator(genome_name, ref_dest, markdup_bam, vcf_dest, out_dir)

    print(f"\n{Colours.BOLD}────────────────────────── Summary ──────────────────────────────{Colours.RESET}")
    print(f"  Genome name   : {Colours.CYAN}{genome_name}{Colours.RESET}")
    print(f"  Output dir    : {Colours.CYAN}{out_dir}{Colours.RESET}")
    print(f"  BAM           : {markdup_bam.name}")
    print(f"  FASTA         : {ref_dest.name}")
    print(f"  Variants      : {vcf_dest.name if vcf_dest else '(none)'}")
    print(f"{Colours.BOLD}─────────────────────────────────────────────────────────────────{Colours.RESET}\n")


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="build_variant_browser",
        description="Prepare alignment, reference and optional variant files for an Allelescope browser instance.",
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA file")
    parser.add_argument("-v", "--variants", help="Variant file: VCF or VCF.gz, must be sorted")
    parser.add_argument("-g", "--genome-name", default="", help="Genome name; default: first VCF contig, then first BAM reference")
    parser.add_argument("-o", "--out-dir", default="./out_dir", help="Output directory; default: ./out_dir")
    parser.add_argument("-t", "--threads", type=int, default=10, help="Threads for samtools; default: 10")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        build_variant_browser(args)
    except BuildError as exc:
        print(f"{Colours.RED}[ERROR]{Colours.RESET} {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
