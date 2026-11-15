"""File names and text encodings that must also work on Windows."""
import locale
import re
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from allelescope.process_variant import process_variant, site_stem

WINDOWS_FORBIDDEN = re.compile(r'[<>:"/\\|?*\x00-\x1f]')


@pytest.mark.parametrize("chrom", ["chr1", "1", "chrUn_KI270302v1", "NC_000001.11", "chr1_KI270706v1_random", "MT"])
def test_ordinary_names_unchanged(chrom):
    """Existing per-site files keep their names."""
    assert site_stem(chrom, 1200) == f"variant_{chrom}_1200"


@pytest.mark.parametrize("chrom", [
    "HLA-A*01:01:01:01", "chr1:1-100", "a/b", "..\\x", "..", "", "seq name", 'q"<>|?',
])
def test_unsafe_names_are_valid_file_names(chrom):
    stem = site_stem(chrom, 7)
    assert not WINDOWS_FORBIDDEN.search(stem)
    assert stem not in (".", "..") and not stem.endswith((" ", "."))
    assert Path(stem).name == stem          # no path separators survive


def test_unsafe_names_do_not_collide():
    assert site_stem("HLA-A*01:01", 1) != site_stem("HLA-A_01_01", 1)
    assert site_stem("HLA-A*01:01", 1) != site_stem("HLA-A:01*01", 1)
    assert site_stem("HLA-A*01:01", 1) == site_stem("HLA-A*01:01", 1)


def test_position_must_be_numeric():
    with pytest.raises(ValueError):
        site_stem("chr1", "12/../x")


def test_process_variant_with_hla_style_contig(tmp_path):
    """A site on a contig whose name has '*' and ':' is processed and stored safely."""
    chrom = "HLA-A*01:01:01:01"
    seq = "ACGT" * 250
    (tmp_path / "ref.fa").write_text(f">{chrom}\n{seq}\n", encoding="utf-8")
    pysam.faidx(str(tmp_path / "ref.fa"))

    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": chrom, "LN": len(seq)}]})
    reads = []
    for i, (start, mstart) in enumerate([(100, 300), (120, 320), (140, 340)]):
        for r1, s, ms, rev in ((True, start, mstart, False), (False, mstart, start, True)):
            a = pysam.AlignedSegment(header)
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = s
            a.cigarstring = "100M"
            a.query_sequence = seq[s:s + 100]
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.mapping_quality = 60
            a.flag = 0x1 | 0x2 | (0x40 if r1 else 0x80) | (0x10 if rev else 0) | (0x0 if rev else 0x20)
            a.next_reference_id = 0
            a.next_reference_start = ms
            reads.append(a)
    bam = tmp_path / "alignment.markdup.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in sorted(reads, key=lambda r: r.reference_start):
            out.write(r)
    pysam.index(str(bam))

    assert process_variant(tmp_path, chrom, 150, "ref.fa") == 0

    stem = site_stem(chrom, 150)
    names = {p.name for p in (tmp_path / "bam").iterdir()}
    assert f"{stem}.calmd.bam" in names
    assert not any(WINDOWS_FORBIDDEN.search(n) for n in names)
    with pysam.AlignmentFile(str(tmp_path / "bam" / f"{stem}.calmd.bam")) as f:
        assert f.count() == 6


def test_generate_viewer_writes_utf8(tmp_path):
    """The page contains non-ASCII characters; it must not depend on the locale."""
    seq = "ACGT" * 50
    (tmp_path / "ref.fa").write_text(f">chr1\n{seq}\n", encoding="utf-8")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(seq)}]}
    with pysam.AlignmentFile(str(tmp_path / "a.bam"), "wb", header=header):
        pass
    env = {"PYTHONIOENCODING": "utf-8", "LC_ALL": "C", "LANG": "C", "PYTHONUTF8": "0", "PATH": ""}
    # A non-UTF-8 default encoding, as on Windows (cp1252). ASCII is stricter.
    cmd = [sys.executable, "-X", "warn_default_encoding", "-c",
           "import runpy, sys; sys.argv = sys.argv[1:]; runpy.run_module('allelescope.generate_viewer', run_name='__main__')",
           "generate_viewer", "--genome-name", "chr1", "--reference-fasta", str(tmp_path / "ref.fa"),
           "--bam-path", str(tmp_path / "a.bam"), "--output-dir", str(tmp_path)]
    run = subprocess.run(cmd, capture_output=True, text=True, env={**env, "PYTHONPATH": ":".join(sys.path)})
    assert run.returncode == 0, run.stderr[-600:]
    html = (tmp_path / "interactive_variants.html").read_text(encoding="utf-8")
    assert "◀" in html            # the ◀ pane toggle survived
