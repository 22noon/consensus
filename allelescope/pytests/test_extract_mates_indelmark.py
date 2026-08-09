from pathlib import Path

import pysam

from allelescope.extract_mates import extract_mates


def write_source_bam(path: Path) -> None:
    write_read_bam(path, "AACGT", "5M", "1T3")


def write_read_bam(path: Path, sequence: str, cigar: str, md: str) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}

    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        read = pysam.AlignedSegment(bam.header)
        read.query_name = "test-read"
        read.query_sequence = sequence
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 100  # 1-based position 101
        read.mapping_quality = 60
        read.cigarstring = cigar
        read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
        read.set_tag("MD", md)
        bam.write(read)

    pysam.index(str(path))


def test_extract_mates_tags_the_requested_one_based_position(tmp_path: Path) -> None:
    source_bam = tmp_path / "source.bam"
    overlap_bam = tmp_path / "overlap.bam"
    mates_bam = tmp_path / "mates.bam"
    write_source_bam(source_bam)

    allele_counts = extract_mates(
        str(source_bam), str(overlap_bam), str(mates_bam), "chr1", 102, 102
    )

    assert allele_counts == {"MA": 1}
    with pysam.AlignmentFile(overlap_bam, "rb") as bam:
        reads = list(bam)
    assert len(reads) == 1
    assert reads[0].get_tag("XA") == "MA"


def test_extract_mates_tags_an_insertion_and_its_summary(tmp_path: Path) -> None:
    source_bam = tmp_path / "source.bam"
    overlap_bam = tmp_path / "overlap.bam"
    mates_bam = tmp_path / "mates.bam"
    write_read_bam(source_bam, "ACGGTTA", "2M2I3M", "5")

    allele_counts = extract_mates(
        str(source_bam), str(overlap_bam), str(mates_bam), "chr1", 103, 103
    )

    assert allele_counts == {"IGG": 1}
    with pysam.AlignmentFile(overlap_bam, "rb") as bam:
        read = next(bam)
    assert read.get_tag("XA") == "IGG"
    assert read.get_tag("XI") == 2
    assert read.get_tag("XD") == 0
    assert read.get_tag("XN") == 2
    assert read.get_tag("XS") == "I2"


def test_extract_mates_tags_a_deletion_and_its_summary(tmp_path: Path) -> None:
    source_bam = tmp_path / "source.bam"
    overlap_bam = tmp_path / "overlap.bam"
    mates_bam = tmp_path / "mates.bam"
    write_read_bam(source_bam, "ACTTA", "2M3D3M", "2^GTA3")

    allele_counts = extract_mates(
        str(source_bam), str(overlap_bam), str(mates_bam), "chr1", 104, 104
    )

    assert allele_counts == {"D3": 1}
    with pysam.AlignmentFile(overlap_bam, "rb") as bam:
        read = next(bam)
    assert read.get_tag("XA") == "D3"
    assert read.get_tag("XI") == -3
    assert read.get_tag("XD") == 3
    assert read.get_tag("XN") == 0
    assert read.get_tag("XS") == "D3"


def test_extract_mates_tags_a_reference_skip_without_indel_summary(tmp_path: Path) -> None:
    source_bam = tmp_path / "source.bam"
    overlap_bam = tmp_path / "overlap.bam"
    mates_bam = tmp_path / "mates.bam"
    write_read_bam(source_bam, "ACTTA", "2M3N3M", "5")

    allele_counts = extract_mates(
        str(source_bam), str(overlap_bam), str(mates_bam), "chr1", 104, 104
    )

    assert allele_counts == {"D3": 1}
    with pysam.AlignmentFile(overlap_bam, "rb") as bam:
        read = next(bam)
    assert read.get_tag("XA") == "D3"
    assert read.get_tag("XI") == 0
    assert read.get_tag("XD") == 0
    assert read.get_tag("XN") == 0
    assert read.get_tag("XS") == ""


def test_extract_mates_writes_empty_outputs_when_region_has_no_reads(tmp_path: Path) -> None:
    source_bam = tmp_path / "source.bam"
    overlap_bam = tmp_path / "overlap.bam"
    mates_bam = tmp_path / "mates.bam"
    write_source_bam(source_bam)

    allele_counts = extract_mates(
        str(source_bam), str(overlap_bam), str(mates_bam), "chr1", 500, 500
    )

    assert allele_counts == {}
    with pysam.AlignmentFile(overlap_bam, "rb") as bam:
        assert list(bam) == []
    with pysam.AlignmentFile(mates_bam, "rb") as bam:
        assert list(bam) == []
