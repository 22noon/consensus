"""Detection of BAM preprocessing already done (allelescope/bam_state.py)."""
import os
import time

import pysam
import pytest

from allelescope import build_variant_browser as bvb
from allelescope.bam_state import (PLAN_COPY, PLAN_FULL, PLAN_MARKDUP, PLAN_SORT,
                                   existing_index, inspect_bam)


def write_bam(path, *, so="coordinate", pg=(), mate_tags=False, order="sorted", paired=True, n=40, dup_every=0):
    header = {"HD": {"VN": "1.6", "SO": so}, "SQ": [{"SN": "chr1", "LN": 10_000}], "PG": list(pg)}
    hdr = pysam.AlignmentHeader.from_dict(header)
    reads = []
    for i in range(n):
        start = 100 + 10 * i
        mates = (True, False) if paired else (True,)
        for is_r1 in mates:
            r = pysam.AlignedSegment(hdr)
            r.query_name = f"q{i}"
            r.reference_id = 0
            r.reference_start = start if is_r1 else start + 200
            r.cigarstring = "50M"
            r.query_sequence = "A" * 50
            r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            r.mapping_quality = 60
            if paired:
                r.flag = 0x1 | 0x2 | (0x40 if is_r1 else 0x80 | 0x10) | (0x20 if is_r1 else 0)
                r.next_reference_id = 0
                r.next_reference_start = start + 200 if is_r1 else start
                if mate_tags:
                    r.set_tag("MC", "50M")
                    r.set_tag("ms", 1500)
            if dup_every and i % dup_every == 0:
                r.flag |= 0x400
            reads.append(r)
    if order == "sorted":
        reads.sort(key=lambda r: r.reference_start)
    elif order == "reversed":
        reads.sort(key=lambda r: -r.reference_start)
    with pysam.AlignmentFile(str(path), "wb", header=hdr) as out:
        for r in reads:
            out.write(r)
    return path


MARKDUP_PG = [
    {"ID": "samtools.3", "PN": "samtools", "VN": "1.19.2", "CL": "samtools markdup -@ 4 - out.bam"},
    {"ID": "MarkDuplicates", "PN": "MarkDuplicates", "CL": "MarkDuplicates INPUT=a.bam OUTPUT=b.bam"},
    {"ID": "GATK MarkDuplicatesSpark", "PN": "GATK MarkDuplicatesSpark"},
    {"ID": "sambamba", "PN": "sambamba", "CL": "sambamba markdup in.bam out.bam"},
    {"ID": "bamsormadup", "PN": "bamsormadup"},
    {"ID": "samblaster", "PN": "samblaster", "CL": "samblaster -i in.sam"},
]


@pytest.mark.parametrize("pg", MARKDUP_PG, ids=lambda pg: pg["ID"])
def test_duplicate_marking_tools_are_recognised(tmp_path, pg):
    state = inspect_bam(write_bam(tmp_path / "a.bam", pg=[pg]))
    assert state.duplicates_marked
    assert state.plan() == PLAN_COPY


@pytest.mark.parametrize("pg", [
    {"ID": "bwa", "PN": "bwa", "CL": "bwa mem ref.fa r1.fq r2.fq"},
    {"ID": "samtools", "PN": "samtools", "CL": "samtools sort -o out.bam in.bam"},
    {"ID": "samtools.1", "PN": "samtools", "CL": "samtools view -b -q 20 in.bam"},
], ids=["bwa", "samtools sort", "samtools view"])
def test_other_programs_do_not_count_as_duplicate_marking(tmp_path, pg):
    assert not inspect_bam(write_bam(tmp_path / "a.bam", pg=[pg])).duplicates_marked


def test_duplicate_flags_alone_are_not_trusted(tmp_path):
    """Flags without a known tool are reported, but the BAM is still processed."""
    state = inspect_bam(write_bam(tmp_path / "a.bam", dup_every=5))
    assert state.duplicates_seen > 0 and not state.duplicates_marked
    assert state.plan() == PLAN_FULL
    assert any("flagged records" in line for line in state.describe())


def test_plans(tmp_path):
    md = [MARKDUP_PG[0]]
    assert inspect_bam(write_bam(tmp_path / "1.bam")).plan() == PLAN_FULL
    assert inspect_bam(write_bam(tmp_path / "2.bam", mate_tags=True)).plan() == PLAN_MARKDUP
    assert inspect_bam(write_bam(tmp_path / "3.bam", pg=md)).plan() == PLAN_COPY
    assert inspect_bam(write_bam(tmp_path / "4.bam", pg=md, so="queryname", order="name")).plan() == PLAN_SORT
    # sorted header but out-of-order records: not trusted
    lying = inspect_bam(write_bam(tmp_path / "5.bam", mate_tags=True, order="reversed"))
    assert lying.header_sorted and not lying.is_sorted and lying.plan() == PLAN_FULL
    # unsorted with mate tags still needs the full pipeline (markdup needs sorted input)
    assert inspect_bam(write_bam(tmp_path / "6.bam", so="unsorted", mate_tags=True, order="name")).plan() == PLAN_FULL
    # single-end without duplicate marking: full pipeline, as before
    assert inspect_bam(write_bam(tmp_path / "7.bam", paired=False)).plan() == PLAN_FULL


def test_partial_mate_tags_are_not_enough(tmp_path):
    path = write_bam(tmp_path / "a.bam", mate_tags=True)
    with pysam.AlignmentFile(str(path)) as fin, pysam.AlignmentFile(str(tmp_path / "b.bam"), "wb", template=fin) as fout:
        for i, r in enumerate(fin):
            if i == 7:
                r.set_tag("ms", None)
            fout.write(r)
    state = inspect_bam(tmp_path / "b.bam")
    assert state.with_mate_tags == state.paired_checked - 1
    assert state.plan() == PLAN_FULL


def test_existing_index(tmp_path):
    bam = write_bam(tmp_path / "a.bam")
    assert existing_index(bam) is None
    pysam.index(str(bam))
    assert existing_index(bam) == tmp_path / "a.bam.bai"
    # an index older than the BAM is not reused
    old = time.time() - 3600
    os.utime(tmp_path / "a.bam.bai", (old, old))
    assert existing_index(bam) is None


def _run_build(tmp_path, bam, *extra):
    (tmp_path / "ref.fa").write_text(">chr1\n" + "A" * 10_000 + "\n", encoding="utf-8")
    out = tmp_path / "inst"
    rc = bvb.main(["-b", str(bam), "-r", str(tmp_path / "ref.fa"), "-o", str(out), "-t", "1", *extra])
    return rc, out


def test_build_copies_processed_bam_and_reuses_index(tmp_path):
    bam = write_bam(tmp_path / "done.bam", pg=[MARKDUP_PG[1]])
    pysam.index(str(bam))
    rc, out = _run_build(tmp_path, bam)
    assert rc == 0
    dest = out / "alignment.markdup.bam"
    assert dest.read_bytes() == bam.read_bytes()
    assert (out / "alignment.markdup.bam.bai").read_bytes() == (tmp_path / "done.bam.bai").read_bytes()
    assert os.stat(dest).st_ino != os.stat(bam).st_ino


def test_build_link_option_hard_links(tmp_path):
    bam = write_bam(tmp_path / "done.bam", pg=[MARKDUP_PG[0]])
    rc, out = _run_build(tmp_path, bam, "--link")
    assert rc == 0
    assert os.stat(out / "alignment.markdup.bam").st_ino == os.stat(bam).st_ino
    assert (out / "alignment.markdup.bam.bai").exists()


def test_build_skip_refuses_unsorted(tmp_path):
    bam = write_bam(tmp_path / "raw.bam", so="unsorted", order="name")
    rc, _ = _run_build(tmp_path, bam, "--bam-processing", "skip")
    assert rc == 1
