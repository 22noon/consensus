from allelescope.extract_mates import extract_mates, Tag_bits
from bam_helper  import * 


def test_extract_mates_sets_non_spanning_mate(tmp_path):

    in_bam = tmp_path / "input.bam"
    out_bam = tmp_path / "overlap.bam"
    mate_bam = tmp_path / "mates.bam"

    r1 = make_bam_read(
        "pair1",
        is_read1=True,
        start=95,          # overlaps target 100
        mate_start=300,
        reverse=False,
        mate_reverse=True,
    )

    r2 = make_bam_read(
        "pair1",
        is_read1=False,
        start=300,         # non-overlapping mate
        mate_start=95,
        reverse=True,
        mate_reverse=False,
    )

    write_and_index_bam(in_bam, [r1, r2])

    extract_mates(
        str(in_bam),
        str(out_bam),
        str(mate_bam),
        "chr1",
        100,
        100,
    )

    with pysam.AlignmentFile(mate_bam, "rb") as bam:
        mates = list(bam)

    assert len(mates) == 1

    xo = mates[0].get_tag("XO")

    assert xo & Tag_bits["NON_SPANNING_MATE"]

def test_extract_mates_sets_split(tmp_path):

    in_bam = tmp_path / "input.bam"
    out_bam = tmp_path / "overlap.bam"
    mate_bam = tmp_path / "mates.bam"

    r1 = make_bam_read(
        "split1",
        is_read1=True,
        start=95,
        mate_start=300,
        reverse=False,
        mate_reverse=True,
        sa_tag="chr1,300,+,10M,60,0;",
    )

    r2 = make_bam_read(
        "split1",
        is_read1=False,
        start=300,
        mate_start=95,
        reverse=True,
        mate_reverse=False,
    )

    write_and_index_bam(in_bam, [r1, r2])

    extract_mates(
        str(in_bam),
        str(out_bam),
        str(mate_bam),
        "chr1",
        100,
        100,
    )

    with pysam.AlignmentFile(out_bam, "rb") as bam:
        reads = list(bam)

    split_read = next(r for r in reads if r.query_name == "split1")

    xo = split_read.get_tag("XO")

    assert xo & Tag_bits["SPLIT"]

def test_extract_mates_sets_splitx(tmp_path):

    in_bam = tmp_path / "input.bam"
    out_bam = tmp_path / "overlap.bam"
    mate_bam = tmp_path / "mates.bam"

    primary = make_bam_read(
        "splitread",
        is_read1=True,
        start=95,
        mate_start=500,
        reverse=False,
        mate_reverse=True,
        sa_tag="chr1,301,+,10M,60,0;",
    )

    supplementary = make_bam_read(
        "splitread",
        is_read1=True,
        start=300,
        mate_start=500,
        reverse=False,
        mate_reverse=True,
    )
    supplementary.is_supplementary = True

    mate = make_bam_read(
        "splitread",
        is_read1=False,
        start=500,
        mate_start=95,
        reverse=True,
        mate_reverse=False,
    )

    write_and_index_bam(
        in_bam,
        [primary, supplementary, mate],
    )

    extract_mates(
        str(in_bam),
        str(out_bam),
        str(mate_bam),
        "chr1",
        100,
        100,
    )

    with pysam.AlignmentFile(mate_bam, "rb") as bam:
        reads = list(bam)

    splitx_read = next(
        r for r in reads
        if r.reference_start == 300
    )

    xo = splitx_read.get_tag("XO")

    assert xo & Tag_bits["SPLITX"]


