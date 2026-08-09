"""
pytest suite for extract_mates().

Builds small synthetic BAM files with pysam to exercise:
  - a normal pair where both mates cross the region        -> OVERLAP set
  - a split read where only the primary crosses the region  -> no OVERLAP,
    but the non-crossing supplementary and the distant mate
    both get recovered into mates_bam (SPLITX / NON_SPANNING_MATE)
  - a "duplication-like" read where two of its own segments
    both cross the region                                   -> both get
    OVERLAP, allele_counts is deterministic and NOT double-counted,
    and no redundant mate-fetch happens (both already written directly)
  - an unpaired read crossing the region                     -> never OVERLAP
  - a mate that is expected but not actually present in the BAM
    (simulates a genuinely unmappable/missing mate)          -> no crash,
    not written to mates_bam, reported as missing

extract_mates.py's helper-function dependencies (Tag_bits, tag_read, etc.)
are monkeypatched to the deterministic stand-ins in fake_helpers.py so the
tests don't depend on the real project's implementations.
"""
import sys
import os

import pysam
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import allelescope.extract_mates as em
import fake_helpers as fh


CHROM = "chr1"
CHROM_LEN = 2_000_000


@pytest.fixture(autouse=True)
def patch_helpers(monkeypatch):
    """Wire extract_mates.py's helper dependencies to the deterministic fakes."""
    monkeypatch.setattr(em, "Tag_bits", fh.Tag_bits, raising=False)
    monkeypatch.setattr(em, "tag_read", fh.tag_read, raising=False)
    monkeypatch.setattr(em, "calculate_softclip_lengths", fh.calculate_softclip_lengths, raising=False)
    monkeypatch.setattr(em, "parse_sa_tag", fh.parse_sa_tag, raising=False)
    monkeypatch.setattr(em, "get_allele_at_base", fh.get_allele_at_base, raising=False)
    monkeypatch.setattr(em, "get_indels_at_base", fh.get_indels_at_base, raising=False)
    monkeypatch.setattr(em, "summarise_indels", fh.summarise_indels, raising=False)


def _header():
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
    }


def mk_read(
    header_obj,
    qname,
    pos,                      # 0-based reference_start
    seq,
    *,
    is_paired=True,
    is_read1=False,
    is_read2=False,
    is_reverse=False,
    is_supplementary=False,
    mate_pos=None,            # 0-based; None -> mate_is_unmapped
    mate_reverse=False,
    sa_tag=None,
    cigar=None,
    mapq=60,
):
    r = pysam.AlignedSegment(header_obj)
    r.query_name = qname
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    r.reference_id = 0
    r.reference_start = pos
    r.mapping_quality = mapq
    r.cigarstring = cigar or f"{len(seq)}M"

    r.is_paired = is_paired
    r.is_read1 = is_read1
    r.is_read2 = is_read2
    r.is_reverse = is_reverse
    r.is_supplementary = is_supplementary
    r.is_unmapped = False

    if is_paired and mate_pos is not None:
        r.next_reference_id = 0
        r.next_reference_start = mate_pos
        r.mate_is_reverse = mate_reverse
        r.mate_is_unmapped = False
    else:
        r.next_reference_id = -1
        r.next_reference_start = -1
        r.mate_is_unmapped = True

    if sa_tag is not None:
        r.set_tag("SA", sa_tag, value_type="Z")

    return r


def build_bam(path, reads):
    header_dict = _header()
    with pysam.AlignmentFile(path, "wb", header=header_dict) as fh_out:
        for r in sorted(reads, key=lambda x: x.reference_start):
            fh_out.write(r)
    pysam.index(path)


SEQ60 = "ACGTACGTAC" * 6  # 60bp, base at any offset is deterministic


def run_extract(tmp_path, reads, start=101, end=101):
    src = str(tmp_path / "src.bam")
    out = str(tmp_path / "overlap.bam")
    mates = str(tmp_path / "mates.bam")
    build_bam(src, reads)
    allele_counts = em.extract_mates(src, out, mates, CHROM, start, end)
    return allele_counts, out, mates


def read_names_and_flags(path):
    """Return {(query_name, reference_start): XO_tag_or_None} for every record."""
    result = {}
    with pysam.AlignmentFile(path, "rb") as f:
        for r in f.fetch(until_eof=True):
            xo = r.get_tag("XO") if r.has_tag("XO") else None
            result[(r.query_name, r.reference_start)] = xo
    return result


# ---------------------------------------------------------------------------
# General case: both mates cross the region -> OVERLAP set on both
# ---------------------------------------------------------------------------
def test_both_mates_crossing_get_overlap(tmp_path):
    seq = SEQ60
    r1 = mk_read(pysam.AlignmentHeader.from_dict(_header()), "pairA", 90, seq,
                 is_read1=True, mate_pos=80, mate_reverse=True)
    r2 = mk_read(pysam.AlignmentHeader.from_dict(_header()), "pairA", 80, seq,
                 is_read2=True, is_reverse=True, mate_pos=90)

    allele_counts, out_bam, mates_bam = run_extract(tmp_path, [r1, r2])

    overlap = read_names_and_flags(out_bam)
    assert len(overlap) == 2
    for xo in overlap.values():
        assert xo is not None and (xo & fh.Tag_bits["OVERLAP"])

    # one allele observation per read pair member that covers target_pos --
    # here each mate is its own physical read, so 2 distinct query names -> 2 counts total? No:
    # allele_counts is keyed by allele string, one increment per *query_name*.
    assert sum(allele_counts.values()) == 1  # both r1 and r2 share qname "pairA" -> ONE fragment
    assert set(allele_counts.keys()) == {seq[100 - 90]}  # base at target for the primary (r1)


# ---------------------------------------------------------------------------
# Only one mate crosses the region; the other mate and a split partner are
# fetched separately and must NOT get OVERLAP.
# ---------------------------------------------------------------------------
def test_split_primary_crosses_mate_and_supplementary_recovered(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_header())
    seq = SEQ60

    # Primary R1: crosses target (90-150), points via SA at its own
    # supplementary at 0-based 1000 (SA pos is 1-based -> 1001).
    primary = mk_read(hdr, "splitB", 90, seq, is_read1=True,
                       mate_pos=2000, mate_reverse=True,
                       sa_tag=f"{CHROM},1001,+,{len(seq)}M,60,0;")

    # Supplementary R1: does NOT cross the target region at all.
    supp = mk_read(hdr, "splitB", 1000, seq, is_read1=True, is_supplementary=True,
                    mate_pos=2000, mate_reverse=True,
                    sa_tag=f"{CHROM},91,+,{len(seq)}M,60,0;")

    # R2 mate: far away, also does not cross the target region.
    mate_read = mk_read(hdr, "splitB", 2000, seq, is_read2=True,
                         mate_pos=90)

    allele_counts, out_bam, mates_bam = run_extract(tmp_path, [primary, supp, mate_read])

    overlap = read_names_and_flags(out_bam)
    # Only the primary physically overlaps region [100,101) -> only it is
    # written directly, and since its mate (R2) never crosses the region,
    # crosses_as_pair() is False, so OVERLAP must NOT be set.
    assert list(overlap.keys()) == [("splitB", 90)]
    assert overlap[("splitB", 90)] is not None
    assert not (overlap[("splitB", 90)] & fh.Tag_bits["OVERLAP"])

    mates = read_names_and_flags(mates_bam)
    assert ("splitB", 1000) in mates  # supplementary recovered
    assert mates[("splitB", 1000)] & fh.Tag_bits["SPLITX"]
    assert ("splitB", 2000) in mates  # R2 mate recovered
    assert mates[("splitB", 2000)] & fh.Tag_bits["NON_SPANNING_MATE"]

    assert sum(allele_counts.values()) == 1  # one fragment, counted once


# ---------------------------------------------------------------------------
# "Duplication" case: two of the read's own segments both cross the region.
# Both get OVERLAP (their pair also crosses), allele_counts picks the
# primary's call deterministically and only counts once, and no redundant
# mate-fetch happens since both segments were already written directly.
# ---------------------------------------------------------------------------
def test_duplication_both_segments_cross_deterministic_allele(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_header())
    seq_primary = "A" * 60          # base at target_pos (100) = 'A'
    seq_supp    = "T" * 60          # base at target_pos (100) = 'T' (different allele)

    primary = mk_read(hdr, "dupC", 90, seq_primary, is_read1=True,
                       mate_pos=80, mate_reverse=True,
                       sa_tag=f"{CHROM},96,+,60M,60,0;")   # points at supp (0-based 95 -> 1-based 96)
    supp = mk_read(hdr, "dupC", 95, seq_supp, is_read1=True, is_supplementary=True,
                    mate_pos=80, mate_reverse=True,
                    sa_tag=f"{CHROM},91,+,60M,60,0;")      # points at primary (0-based 90 -> 1-based 91)
    mate_read = mk_read(hdr, "dupC", 80, SEQ60, is_read2=True, is_reverse=True,
                         mate_pos=90)

    allele_counts, out_bam, mates_bam = run_extract(tmp_path, [primary, supp, mate_read])

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("dupC", 90), ("dupC", 95), ("dupC", 80)}
    for xo in overlap.values():
        assert xo & fh.Tag_bits["OVERLAP"]

    # Deterministic: primary's allele ('A') wins over supplementary's ('T'),
    # and it's counted exactly once for the whole fragment.
    assert allele_counts == {"A": 1}

    # Both split partner and mate were already written directly in pass 1 --
    # nothing left to fetch, mates_bam should be empty.
    with pysam.AlignmentFile(mates_bam, "rb") as f:
        assert list(f.fetch(until_eof=True)) == []


# ---------------------------------------------------------------------------
# Unpaired read crossing the region must never get OVERLAP.
# ---------------------------------------------------------------------------
def test_unpaired_read_never_gets_overlap(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_header())
    r = mk_read(hdr, "unpairedD", 90, SEQ60, is_paired=False)

    allele_counts, out_bam, mates_bam = run_extract(tmp_path, [r])

    overlap = read_names_and_flags(out_bam)
    assert list(overlap.keys()) == [("unpairedD", 90)]
    assert not (overlap[("unpairedD", 90)] & fh.Tag_bits["OVERLAP"])
    assert sum(allele_counts.values()) == 1


# ---------------------------------------------------------------------------
# Expected mate position has no actual read there -> should not crash,
# should not appear in mates_bam, and should be reported as missing.
# ---------------------------------------------------------------------------
def test_missing_mate_reported_not_written(tmp_path, capsys):
    hdr = pysam.AlignmentHeader.from_dict(_header())
    # R1 crosses target and claims a mate at 500000 that we never write.
    r1 = mk_read(hdr, "farE", 90, SEQ60, is_read1=True, mate_pos=500_000)

    allele_counts, out_bam, mates_bam = run_extract(tmp_path, [r1])

    overlap = read_names_and_flags(out_bam)
    assert list(overlap.keys()) == [("farE", 90)]
    # only one mate (R1) has a primary crossing -> not both -> no OVERLAP
    assert not (overlap[("farE", 90)] & fh.Tag_bits["OVERLAP"])

    with pysam.AlignmentFile(mates_bam, "rb") as f:
        assert list(f.fetch(until_eof=True)) == []

    captured = capsys.readouterr()
    assert "not found" in captured.out
    assert sum(allele_counts.values()) == 1


# ---------------------------------------------------------------------------
# Edge case: mate lies exactly at the WINDOW boundary and must still be found.
# ---------------------------------------------------------------------------
def test_mate_at_window_boundary_is_found(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_header())
    # WINDOW = 300; put the mate's alignment 0-based start at +290 so its
    # fetch window (pos-300, pos+300) is guaranteed to include it.
    mate_pos = 90 + 290
    r1 = mk_read(hdr, "edgeF", 90, SEQ60, is_read1=True, mate_pos=mate_pos)
    r2 = mk_read(hdr, "edgeF", mate_pos, SEQ60, is_read2=True, mate_pos=90)

    allele_counts, out_bam, mates_bam = run_extract(tmp_path, [r1, r2])

    overlap = read_names_and_flags(out_bam)
    # r2 doesn't cross [100,101) (starts at 380) so only r1 is written directly.
    assert list(overlap.keys()) == [("edgeF", 90)]

    mates = read_names_and_flags(mates_bam)
    assert ("edgeF", mate_pos) in mates
    assert mates[("edgeF", mate_pos)] & fh.Tag_bits["NON_SPANNING_MATE"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
