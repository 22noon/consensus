"""Pass 2 of extract_mates: merged mate fetches find the same records."""
import pysam
import pytest

from allelescope.extract_mates import _merge_positions, extract_mates, Tag_bits


def test_merge_positions_basic():
    assert _merge_positions([], 1000) == []
    assert _merge_positions([500], 1000) == [(500, 501)]
    # unsorted input, duplicates, and gaps on both sides of the threshold
    assert _merge_positions({3000, 100, 900, 100, 1100}, 1000) == [(100, 1101), (3000, 3001)]
    assert _merge_positions([0, 1000, 2001], 1000) == [(0, 1001), (2001, 2002)]
    assert _merge_positions([0, 1000], 999) == [(0, 1), (1000, 1001)]


HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"},
          "SQ": [{"SN": "chr1", "LN": 100_000}, {"SN": "chr2", "LN": 100_000}]}


def _read(hdr, name, chrom, start, r1, rev, mchrom, mstart, mrev, cigar="50M", extra=0, sa=None):
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.reference_id = 0 if chrom == "chr1" else 1
    r.reference_start = start
    r.cigarstring = cigar
    qlen = r.infer_query_length()
    r.query_sequence = "A" * qlen
    r.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    r.mapping_quality = 60
    r.flag = 0x1 | 0x2 | (0x40 if r1 else 0x80) | (0x10 if rev else 0) | (0x20 if mrev else 0) | extra
    r.next_reference_id = 0 if mchrom == "chr1" else 1
    r.next_reference_start = mstart
    if sa:
        r.set_tag("SA", sa)
    return r


def test_mates_found_near_far_and_on_other_sequence(tmp_path):
    """Many nearby mates, one mate 50 kb away, one on chr2 and a split partner."""
    hdr = pysam.AlignmentHeader.from_dict(HEADER)
    reads, expected = [], set()
    site = 10_000
    for i in range(200):                     # dense pairs around the site
        start = site - 40 + (i % 30)
        mstart = start + 200 + i             # 200 distinct mate positions
        reads.append(_read(hdr, f"p{i}", "chr1", start, True, False, "chr1", mstart, True))
        reads.append(_read(hdr, f"p{i}", "chr1", mstart, False, True, "chr1", start, False))
        expected.add((f"p{i}", "chr1", mstart))
    # a filler read at a mate position with a different name must not be taken
    reads.append(_read(hdr, "other", "chr1", site + 200, True, False, "chr1", site + 900, True))
    # far mate on the same sequence and a mate on another sequence
    reads.append(_read(hdr, "far", "chr1", site - 10, True, False, "chr1", 60_000, True))
    reads.append(_read(hdr, "far", "chr1", 60_000, False, True, "chr1", site - 10, False))
    expected.add(("far", "chr1", 60_000))
    reads.append(_read(hdr, "trans", "chr1", site - 20, True, False, "chr2", 5_000, True))
    reads.append(_read(hdr, "trans", "chr2", 5_000, False, True, "chr1", site - 20, False))
    expected.add(("trans", "chr2", 5_000))
    # split read: primary crosses the site, supplementary on chr2
    reads.append(_read(hdr, "split", "chr1", site - 30, True, False, "chr1", site + 400, True,
                       cigar="40M10S", sa="chr2,20001,+,40S10M,60,0;"))
    reads.append(_read(hdr, "split", "chr2", 20_000, True, False, "chr1", site + 400, True,
                       cigar="40H10M", extra=0x800, sa="chr1,9971,+,40M10S,60,0;"))
    reads.append(_read(hdr, "split", "chr1", site + 400, False, True, "chr1", site - 30, False))
    expected |= {("split", "chr2", 20_000), ("split", "chr1", site + 400)}

    src = tmp_path / "in.bam"
    with pysam.AlignmentFile(str(src), "wb", header=hdr) as out:
        for r in sorted(reads, key=lambda r: (r.reference_id, r.reference_start)):
            out.write(r)
    pysam.index(str(src))

    mates = tmp_path / "mates.bam"
    extract_mates(str(src), str(tmp_path / "overlap.bam"), str(mates), "chr1", site, site)

    with pysam.AlignmentFile(str(mates), "rb", check_sq=False) as f:
        got = {(r.query_name, r.reference_name, r.reference_start): r for r in f}

    assert set(got) == expected
    assert got[("split", "chr2", 20_000)].get_tag("XO") & Tag_bits["SPLITX"]
    assert got[("far", "chr1", 60_000)].get_tag("XO") & Tag_bits["NON_SPANNING_MATE"]
    assert not any(name == "other" for name, _, _ in got)
