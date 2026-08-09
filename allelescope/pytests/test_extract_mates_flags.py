import pysam
from allelescope.extract_mates import (
    extract_mates,
    Tag_bits,
)

def make_read(
    qname,
    is_read1,
    start,
    mate_start,
    reverse=False,
    mate_reverse=True,
):
    r = pysam.AlignedSegment()
    r.query_name = qname
    r.query_sequence = "AAAAAAAAAA"

    r.flag = 99 if is_read1 else 147
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = "10M"

    r.is_paired = True
    r.is_proper_pair = True
    r.is_read1 = is_read1
    r.is_read2 = not is_read1

    r.is_reverse = reverse
    r.mate_is_reverse = mate_reverse

    r.next_reference_id = 0
    r.next_reference_start = mate_start
    r.template_length = 100

    return r


def test_extract_mates_adds_tags(tmp_path):
    input_bam = tmp_path / "input.bam"
    overlap_bam = tmp_path / "overlap.bam"
    mates_bam = tmp_path / "mates.bam"

    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }

    # Proper FR pair spanning chr1:100
    r1 = make_read(
        "pair1",
        is_read1=True,
        start=95,
        mate_start=150,
        reverse=False,
        mate_reverse=True,
    )

    r2 = make_read(
        "pair1",
        is_read1=False,
        start=150,
        mate_start=95,
        reverse=True,
        mate_reverse=False,
    )

    with pysam.AlignmentFile(input_bam, "wb", header=header) as bam:
        bam.write(r1)
        bam.write(r2)

    pysam.index(str(input_bam))

    extract_mates(
        str(input_bam),
        str(overlap_bam),
        str(mates_bam),
        "chr1",
        100,
        100,
    )

    with pysam.AlignmentFile(overlap_bam, "rb") as bam:
        reads = list(bam)

    assert len(reads) == 1

    read = reads[0]

    # Orientation tag
    assert read.has_tag("XR")
    assert read.get_tag("XR") == "F1R2"

    # Flag tag
    assert read.has_tag("XO")
    xo = read.get_tag("XO")

    assert not xo & Tag_bits["OVERLAP"]
    assert xo & Tag_bits["PROPER"]

    # Allele/indel tags added
    assert read.has_tag("XA")
    assert read.has_tag("XI")
    assert read.has_tag("XD")
    assert read.has_tag("XN")
    assert read.has_tag("XS")

    # Soft clip tags added
    assert read.has_tag("SL")
    assert read.has_tag("SR")

    assert read.get_tag("SL") == 0
    assert read.get_tag("SR") == 0
