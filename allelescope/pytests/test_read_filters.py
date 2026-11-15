"""Tests for the server-side read filters (allelescope/read_filters.py)."""
import pysam
import pytest

from allelescope.baq_filter import baq_at_site
from allelescope.read_filters import ReadFilter, filter_reads, optional_int

SITE = 100          # 1-based site position
SITE0 = SITE - 1    # 0-based

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 10_000}],
})


def mk(name, start=80, cigar="40M", mapq=60, nm=0, xo=0, xa="RA", flag=0x1 | 0x2 | 0x40,
       qual=30, bq=None, sl=0, sr=0):
    """One alignment record carrying the tags Allelescope writes to per-site BAMs."""
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = start
    r.cigarstring = cigar
    qlen = r.infer_query_length()
    r.query_sequence = "A" * qlen
    r.query_qualities = pysam.qualitystring_to_array(chr(33 + qual) * qlen)
    r.mapping_quality = mapq
    r.flag = flag
    r.set_tag("NM", nm, value_type="i")
    r.set_tag("XO", xo, value_type="i")
    r.set_tag("XA", xa, value_type="Z")
    r.set_tag("SL", sl, value_type="i")
    r.set_tag("SR", sr, value_type="i")
    if bq is not None:
        r.set_tag("BQ", bq, value_type="Z")
    return r


def write_bam(tmp_path, reads):
    path = tmp_path / "site.bam"
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as bam:
        for r in sorted(reads, key=lambda x: x.reference_start):
            bam.write(r)
    return path


def names(tmp_path, reads, **params):
    params.setdefault("Chrom", "chr1")
    params.setdefault("Pos", str(SITE))
    _, kept = filter_reads(write_bam(tmp_path, reads), params)
    return sorted(r.query_name for r in kept)


# ---------------------------------------------------------------- defaults

def test_optional_int():
    assert optional_int(None) is None
    assert optional_int("") is None
    assert optional_int("  ") is None
    assert optional_int("undefined") is None
    assert optional_int("abc") is None
    assert optional_int("-1") is None
    assert optional_int("0") == 0
    assert optional_int("12") == 12
    assert optional_int(7) == 7


@pytest.mark.parametrize("blank", ["", "undefined", None])
def test_no_limits_keeps_everything(tmp_path, blank):
    """With every numeric limit left blank nothing is filtered."""
    reads = [
        mk("perfect"),
        mk("lowmapq", mapq=0),
        mk("mismatched", nm=9),
        mk("clipped", cigar="10S30M", start=90, sl=10),
        mk("lowbaq", qual=2),
    ]
    params = {k: blank for k in
              ("minMapQ", "SoftClip", "SoftClipF", "EditDistance", "EditDistanceF", "BAQ", "BAQF")}
    params.update(Flagf="0", FlagF="0", Tagf="0", TagF="0", XAFilter="")
    params = {k: v for k, v in params.items() if v is not None}
    assert names(tmp_path, reads, **params) == ["clipped", "lowbaq", "lowmapq", "mismatched", "perfect"]


def test_zero_minimums_keep_everything(tmp_path):
    reads = [mk("a", mapq=0, qual=0), mk("b", nm=3)]
    assert names(tmp_path, reads, minMapQ="0", SoftClip="0", EditDistance="0", BAQ="0") == ["a", "b"]


def test_zero_maximum_is_a_real_limit(tmp_path):
    """A maximum of 0 is applied, it does not mean 'off'."""
    reads = [mk("perfect"), mk("mismatched", nm=1), mk("clipped", cigar="5S35M", start=85, sl=5)]
    assert names(tmp_path, reads, EditDistanceF="0") == ["clipped", "perfect"]
    assert names(tmp_path, reads, SoftClipF="0") == ["mismatched", "perfect"]


# ---------------------------------------------------------------- MAPQ

def test_min_mapq(tmp_path):
    reads = [mk("q0", mapq=0), mk("q19", mapq=19), mk("q20", mapq=20), mk("q60", mapq=60)]
    assert names(tmp_path, reads, minMapQ="20") == ["q20", "q60"]
    assert names(tmp_path, reads, minMapQ="1") == ["q19", "q20", "q60"]


# ---------------------------------------------------------------- edit distance / soft clip

def test_edit_distance_min_and_max(tmp_path):
    reads = [mk("nm0", nm=0), mk("nm1", nm=1), mk("nm2", nm=2), mk("nm5", nm=5)]
    assert names(tmp_path, reads, EditDistance="1") == ["nm1", "nm2", "nm5"]
    assert names(tmp_path, reads, EditDistanceF="2") == ["nm0", "nm1", "nm2"]
    assert names(tmp_path, reads, EditDistance="1", EditDistanceF="2") == ["nm1", "nm2"]


def test_soft_clip_min_and_max(tmp_path):
    reads = [
        mk("none"),
        mk("left5", cigar="5S35M", start=85, sl=5),
        mk("both", cigar="5S30M5S", start=85, sl=5, sr=5),
    ]
    assert names(tmp_path, reads, SoftClip="1") == ["both", "left5"]
    assert names(tmp_path, reads, SoftClip="6") == ["both"]
    assert names(tmp_path, reads, SoftClipF="5") == ["left5", "none"]


def test_soft_clip_falls_back_to_cigar_without_tags(tmp_path):
    r = mk("clipped", cigar="5S35M", start=85)
    r.set_tag("SL", None)
    r.set_tag("SR", None)
    assert names(tmp_path, [r, mk("plain")], SoftClip="1") == ["clipped"]


# ---------------------------------------------------------------- flags / XO / alleles

def test_flag_and_xo_bits(tmp_path):
    reads = [
        mk("fwd"),
        mk("rev", flag=0x1 | 0x2 | 0x40 | 0x10),
        mk("dup", flag=0x1 | 0x2 | 0x40 | 0x400),
        mk("mate", xo=0x2),
    ]
    assert names(tmp_path, reads, Flagf="16") == ["rev"]
    assert names(tmp_path, reads, FlagF="16") == ["dup", "fwd", "mate"]
    assert names(tmp_path, reads, FlagF="1024") == ["fwd", "mate", "rev"]
    assert names(tmp_path, reads, Tagf="2") == ["mate"]
    assert names(tmp_path, reads, TagF="2") == ["dup", "fwd", "rev"]


def test_allele_filter(tmp_path):
    reads = [mk("ref", xa="RA"), mk("snp", xa="MC"), mk("del", xa="D7"), mk("none", xa=".")]
    assert names(tmp_path, reads, XAFilter="MC|D7") == ["del", "snp"]
    assert names(tmp_path, reads, XAFilter=".") == ["none"]


# ---------------------------------------------------------------- BAQ

def bq_string(length, offsets):
    """BQ tag value: '@' (offset 0) everywhere except {query index: offset}."""
    chars = ["@"] * length
    for i, off in offsets.items():
        chars[i] = chr(64 + off)
    return "".join(chars)


def test_baq_at_site_uses_bq_offset():
    # read starts at 80, so the site (0-based 99) is query index 19
    r = mk("r", qual=30, bq=bq_string(40, {19: 25}))
    assert baq_at_site(r, SITE0) == 5
    assert baq_at_site(r, SITE0 + 1) == 30


def test_baq_at_site_without_bq_tag_is_base_quality():
    assert baq_at_site(mk("r", qual=17), SITE0) == 17


def test_baq_at_site_not_covering():
    assert baq_at_site(mk("r", start=500), SITE0) is None


def test_baq_at_site_deletion_uses_flanks():
    # 19M 3D 21M from 80: the deletion removes 0-based 99..101 (the site)
    r = mk("r", cigar="19M3D21M", qual=30, bq=bq_string(40, {17: 10, 18: 20}))
    # flanks within +-2 bp: query 17, 18 (ref 97, 98); ref 99..101 deleted
    assert baq_at_site(r, SITE0) == 15


def test_baq_limits(tmp_path):
    reads = [
        mk("good", qual=30),
        mk("lowbaq", qual=30, bq=bq_string(40, {19: 25})),   # BAQ 5 at the site
        mk("high", qual=40),
    ]
    assert names(tmp_path, reads, BAQ="20") == ["good", "high"]
    assert names(tmp_path, reads, BAQF="30") == ["good", "lowbaq"]
    assert names(tmp_path, reads, BAQ="20", BAQF="30") == ["good"]


def test_baq_drops_the_whole_pair(tmp_path):
    """A mate that does not cover the site follows the read that does."""
    r1 = mk("pair", qual=30, bq=bq_string(40, {19: 25}))
    r2 = mk("pair", start=400, flag=0x1 | 0x2 | 0x80 | 0x10, xo=0x2)
    other_mate = mk("ok", start=420, flag=0x1 | 0x2 | 0x80 | 0x10, xo=0x2)
    kept = names(tmp_path, [r1, r2, mk("ok"), other_mate], BAQ="20")
    assert kept == ["ok", "ok"]


def test_baq_decided_independently_of_other_filters(tmp_path):
    """The reverse mate is dropped by BAQ even when the strand filter removes the failing read."""
    r1 = mk("pair", qual=30, bq=bq_string(40, {19: 25}))                 # forward, fails BAQ
    r2 = mk("pair", start=400, flag=0x1 | 0x2 | 0x80 | 0x10, xo=0x2)     # reverse mate
    assert names(tmp_path, [r1, r2], BAQ="20", Flagf="16") == []


def test_from_params_site_position():
    f = ReadFilter.from_params({"Pos": "100", "BAQ": "10"})
    assert f.site_pos0 == 99 and f.uses_baq
    assert not ReadFilter.from_params({"Pos": "100"}).uses_baq
    assert not ReadFilter.from_params({"Pos": "0", "BAQ": "10"}).uses_baq
