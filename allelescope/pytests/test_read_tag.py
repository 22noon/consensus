import pytest
import pysam

from allelescope.extract_mates import tag_read, Tag_bits
def mk_read(is_read1, read_reverse, mate_reverse):
    r = pysam.AlignedSegment()
    r.is_paired = True
    r.is_read1 = is_read1
    r.is_read2 = not is_read1
    r.is_reverse = read_reverse
    r.mate_is_reverse = mate_reverse
    r.is_unmapped = False
    r.mate_is_unmapped = False
    r.is_proper_pair = True
    return r


def test_tag_read_sets_proper():
    read = mk_read(
        is_read1=True,
        read_reverse=False,
        mate_reverse=True,
    )

    tagged = tag_read(read, 0)

    xo = tagged.get_tag("XO")

    assert xo & Tag_bits["PROPER"]
    assert not (xo & Tag_bits["IMPROPER"])
    assert tagged.get_tag("XR") == "F1R2"


def test_tag_read_sets_improper():
    read = mk_read(
        is_read1=True,
        read_reverse=True,
        mate_reverse=False,
    )

    tagged = tag_read(read, 0)

    xo = tagged.get_tag("XO")

    assert xo & Tag_bits["IMPROPER"]
    assert not (xo & Tag_bits["PROPER"])
    assert tagged.get_tag("XR") == "R1F2"


def test_tag_read_sets_split():
    read = mk_read(
        is_read1=True,
        read_reverse=False,
        mate_reverse=True,
    )

    read.set_tag(
        "SA",
        "chr1,200,+,10M,60,0;"
    )

    tagged = tag_read(read, 0)

    xo = tagged.get_tag("XO")

    assert xo & Tag_bits["SPLIT"]
