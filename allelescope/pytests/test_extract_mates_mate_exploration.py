"""
Focused test suite for the MATE-EXPLORATION logic in extract_mates():
i.e., which additional alignments get fetched into mates_bam once the
reads directly overlapping the region of interest are known.

Covers:
  - both mates of a pair overlap the region directly -> no exploration
    needed, mates_bam stays empty
  - only one mate overlaps -> the other is explored and recovered
  - split read: supplementary alignment does NOT overlap -> explored and
    recovered (SPLITX)
  - split read: supplementary alignment DOES overlap directly -> already
    written, no duplicate exploration/fetch
  - a read needing BOTH a split-partner fetch AND a paired-mate fetch at
    once (multiple pending targets for the same query name)
  - mate maps to a different contig entirely
  - mate is flagged unmapped -> never explored, no crash
  - two independent pairs sitting close together -> no cross-contamination
    of mate bookkeeping between different query names
  - a read with two SA entries, one already covered directly and one
    requiring a fetch -> partial dedup

Reuses the BAM-building helpers from test_extract_mates.py.
"""
import sys
import os

import pysam
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import allelescope.extract_mates as em
import fake_helpers as fh

from test_extract_mates import (
    CHROM,
    CHROM_LEN,
    mk_read,
    run_extract,
    read_names_and_flags,
    build_bam,
    SEQ60,
)

CHROM2 = "chr2"


@pytest.fixture(autouse=True)
def patch_helpers(monkeypatch):
    monkeypatch.setattr(em, "Tag_bits", fh.Tag_bits, raising=False)
    monkeypatch.setattr(em, "tag_read", fh.tag_read, raising=False)
    monkeypatch.setattr(em, "calculate_softclip_lengths", fh.calculate_softclip_lengths, raising=False)
    monkeypatch.setattr(em, "parse_sa_tag", fh.parse_sa_tag, raising=False)
    monkeypatch.setattr(em, "get_allele_at_base", fh.get_allele_at_base, raising=False)
    monkeypatch.setattr(em, "get_indels_at_base", fh.get_indels_at_base, raising=False)
    monkeypatch.setattr(em, "summarise_indels", fh.summarise_indels, raising=False)


def _two_contig_header():
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": CHROM, "LN": CHROM_LEN},
            {"SN": CHROM2, "LN": CHROM_LEN},
        ],
    }


def mk_read_multi(
    header_obj,
    qname,
    pos,
    seq,
    *,
    reference_id=0,
    is_paired=True,
    is_read1=False,
    is_read2=False,
    is_reverse=False,
    is_supplementary=False,
    is_secondary=False,
    mate_pos=None,
    mate_reference_id=None,
    mate_reverse=False,
    mate_unmapped=False,
    sa_tag=None,
    cigar=None,
    mapq=60,
):
    """Like test_extract_mates.mk_read, but supports a specific reference_id
    (for multi-contig scenarios) and an explicit mate_unmapped flag."""
    r = pysam.AlignedSegment(header_obj)
    r.query_name = qname
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    r.reference_id = reference_id
    r.reference_start = pos
    r.mapping_quality = mapq
    r.cigarstring = cigar or f"{len(seq)}M"

    r.is_paired = is_paired
    r.is_read1 = is_read1
    r.is_read2 = is_read2
    r.is_reverse = is_reverse
    r.is_supplementary = is_supplementary
    r.is_secondary = is_secondary
    r.is_unmapped = False

    if is_paired and mate_pos is not None and not mate_unmapped:
        r.next_reference_id = mate_reference_id if mate_reference_id is not None else reference_id
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


def run_extract_multi(tmp_path, reads, chrom=CHROM, start=101, end=101):
    src = str(tmp_path / "src.bam")
    out = str(tmp_path / "overlap.bam")
    mates = str(tmp_path / "mates.bam")
    header_dict = _two_contig_header()
    with pysam.AlignmentFile(src, "wb", header=header_dict) as f:
        for r in sorted(reads, key=lambda x: (x.reference_id, x.reference_start)):
            f.write(r)
    pysam.index(src)
    allele_counts = em.extract_mates(src, out, mates, chrom, start, end)
    return allele_counts, out, mates


def mates_bam_records(mates_bam):
    with pysam.AlignmentFile(mates_bam, "rb") as f:
        return {(r.query_name, r.reference_name, r.reference_start): (
            r.get_tag("XO") if r.has_tag("XO") else None
        ) for r in f.fetch(until_eof=True)}


# ---------------------------------------------------------------------------
# 1) Both mates overlap the region directly -> no exploration, mates_bam empty
# ---------------------------------------------------------------------------
def test_both_mates_overlap_directly_no_exploration_needed(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())
    r1 = mk_read_multi(hdr, "pairA", 90, SEQ60, is_read1=True, mate_pos=80, mate_reverse=True)
    r2 = mk_read_multi(hdr, "pairA", 80, SEQ60, is_read2=True, is_reverse=True, mate_pos=90)

    _, out_bam, mates_bam = run_extract_multi(tmp_path, [r1, r2])

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("pairA", 90), ("pairA", 80)}

    assert mates_bam_records(mates_bam) == {}


# ---------------------------------------------------------------------------
# 2) Only one mate overlaps -> the other is explored and recovered
# ---------------------------------------------------------------------------
def test_one_mate_far_away_is_explored_and_recovered(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())
    r1 = mk_read_multi(hdr, "pairB", 90, SEQ60, is_read1=True, mate_pos=5000)
    r2 = mk_read_multi(hdr, "pairB", 5000, SEQ60, is_read2=True, mate_pos=90)

    _, out_bam, mates_bam = run_extract_multi(tmp_path, [r1, r2])

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("pairB", 90)}   # only r1 crosses [100,101)

    mates = mates_bam_records(mates_bam)
    assert (("pairB", CHROM, 5000)) in mates
    assert mates[("pairB", CHROM, 5000)] & fh.Tag_bits["NON_SPANNING_MATE"]


# ---------------------------------------------------------------------------
# 3) Split read: supplementary does NOT overlap -> explored and recovered
# ---------------------------------------------------------------------------
def test_split_supplementary_not_overlapping_is_explored(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())
    primary = mk_read_multi(
        hdr, "splitC", 90, SEQ60, is_read1=True, mate_pos=90,
        sa_tag=f"{CHROM},2001,+,60M,60,0;",
    )
    supp = mk_read_multi(
        hdr, "splitC", 2000, SEQ60, is_read1=True, is_supplementary=True,
        mate_pos=90, sa_tag=f"{CHROM},91,+,60M,60,0;",
    )
    # mate (R2) also overlaps directly, so pair-crossing is satisfied,
    # but that's incidental to this test -- the focus is the supplementary.
    mate_read = mk_read_multi(hdr, "splitC", 85, SEQ60, is_read2=True, is_reverse=True, mate_pos=90)

    _, out_bam, mates_bam = run_extract_multi(tmp_path, [primary, supp, mate_read])

    overlap = read_names_and_flags(out_bam)
    assert ("splitC", 90) in overlap
    assert ("splitC", 2000) not in overlap  # supplementary doesn't cross, not written directly

    mates = mates_bam_records(mates_bam)
    assert ("splitC", CHROM, 2000) in mates
    assert mates[("splitC", CHROM, 2000)] & fh.Tag_bits["SPLITX"]


# ---------------------------------------------------------------------------
# 4) Split read: supplementary DOES overlap directly -> no duplicate fetch
# ---------------------------------------------------------------------------
def test_split_supplementary_overlapping_directly_no_duplicate_fetch(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())
    primary = mk_read_multi(
        hdr, "splitD", 90, SEQ60, is_read1=True, mate_pos=90,
        sa_tag=f"{CHROM},96,+,60M,60,0;",
    )
    supp = mk_read_multi(
        hdr, "splitD", 95, SEQ60, is_read1=True, is_supplementary=True,
        mate_pos=90, sa_tag=f"{CHROM},91,+,60M,60,0;",
    )
    mate_read = mk_read_multi(hdr, "splitD", 85, SEQ60, is_read2=True, is_reverse=True, mate_pos=90)

    _, out_bam, mates_bam = run_extract_multi(tmp_path, [primary, supp, mate_read])

    overlap = read_names_and_flags(out_bam)
    assert {("splitD", 90), ("splitD", 95), ("splitD", 85)} <= set(overlap.keys())

    # Both split partners AND the mate were already written directly ->
    # nothing left to explore for this query name.
    assert mates_bam_records(mates_bam) == {}


# ---------------------------------------------------------------------------
# 5) One query name needs BOTH a split-partner fetch AND a paired-mate fetch
# ---------------------------------------------------------------------------
def test_read_needing_split_and_mate_fetch_simultaneously(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())
    primary = mk_read_multi(
        hdr, "splitE", 90, SEQ60, is_read1=True, mate_pos=8000,
        sa_tag=f"{CHROM},3001,+,60M,60,0;",
    )
    supp = mk_read_multi(
        hdr, "splitE", 3000, SEQ60, is_read1=True, is_supplementary=True,
        mate_pos=8000, sa_tag=f"{CHROM},91,+,60M,60,0;",
    )
    far_mate = mk_read_multi(hdr, "splitE", 8000, SEQ60, is_read2=True, mate_pos=90)

    _, out_bam, mates_bam = run_extract_multi(tmp_path, [primary, supp, far_mate])

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("splitE", 90)}

    mates = mates_bam_records(mates_bam)
    assert ("splitE", CHROM, 3000) in mates
    assert mates[("splitE", CHROM, 3000)] & fh.Tag_bits["SPLITX"]
    assert ("splitE", CHROM, 8000) in mates
    assert mates[("splitE", CHROM, 8000)] & fh.Tag_bits["NON_SPANNING_MATE"]


# ---------------------------------------------------------------------------
# 6) Mate maps to a different contig entirely
# ---------------------------------------------------------------------------
def test_mate_on_different_contig_is_explored(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())
    r1 = mk_read_multi(
        hdr, "pairF", 90, SEQ60, is_read1=True, reference_id=0,
        mate_pos=500, mate_reference_id=1,
    )
    r2 = mk_read_multi(
        hdr, "pairF", 500, SEQ60, is_read2=True, reference_id=1,
        mate_pos=90, mate_reference_id=0,
    )

    _, out_bam, mates_bam = run_extract_multi(tmp_path, [r1, r2])

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("pairF", 90)}

    mates = mates_bam_records(mates_bam)
    assert ("pairF", CHROM2, 500) in mates
    assert mates[("pairF", CHROM2, 500)] & fh.Tag_bits["NON_SPANNING_MATE"]


# ---------------------------------------------------------------------------
# 7) Mate flagged unmapped -> never explored, no crash
# ---------------------------------------------------------------------------
def test_unmapped_mate_never_explored_no_crash(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())
    r1 = mk_read_multi(hdr, "pairG", 90, SEQ60, is_read1=True, mate_unmapped=True)

    _, out_bam, mates_bam = run_extract_multi(tmp_path, [r1])

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("pairG", 90)}
    assert mates_bam_records(mates_bam) == {}


# ---------------------------------------------------------------------------
# 8) Two independent pairs sitting close together -> no cross-contamination
# ---------------------------------------------------------------------------
def test_two_independent_pairs_no_cross_contamination(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())

    # Pair 1: R1 crosses region at 90, its mate is far away at 5000.
    p1_r1 = mk_read_multi(hdr, "pairH1", 90, SEQ60, is_read1=True, mate_pos=5000)
    p1_r2 = mk_read_multi(hdr, "pairH1", 5000, SEQ60, is_read2=True, mate_pos=90)

    # Pair 2: R1 crosses region at 92 (close to pair 1!), its mate is far
    # away at a DIFFERENT position, 6000.
    p2_r1 = mk_read_multi(hdr, "pairH2", 92, SEQ60, is_read1=True, mate_pos=6000)
    p2_r2 = mk_read_multi(hdr, "pairH2", 6000, SEQ60, is_read2=True, mate_pos=92)

    _, out_bam, mates_bam = run_extract_multi(
        tmp_path, [p1_r1, p1_r2, p2_r1, p2_r2]
    )

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("pairH1", 90), ("pairH2", 92)}

    mates = mates_bam_records(mates_bam)
    assert set(mates.keys()) == {
        ("pairH1", CHROM, 5000),
        ("pairH2", CHROM, 6000),
    }
    # Each mate must be tagged for -- and only for -- its own pair.
    assert mates[("pairH1", CHROM, 5000)] & fh.Tag_bits["NON_SPANNING_MATE"]
    assert mates[("pairH2", CHROM, 6000)] & fh.Tag_bits["NON_SPANNING_MATE"]


# ---------------------------------------------------------------------------
# 9) Two SA entries: one already covered directly, one requiring a fetch
# ---------------------------------------------------------------------------
def test_partial_dedup_multiple_sa_entries(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_two_contig_header())

    # Primary at 90 (crosses). SA lists TWO other loci for this fragment:
    # one at 95 (which also crosses and IS written directly -> dedup),
    # and one at 4000 (does not cross -> must be fetched).
    primary = mk_read_multi(
        hdr, "splitI", 90, SEQ60, is_read1=True, mate_pos=90,
        sa_tag=f"{CHROM},96,+,60M,60,0;{CHROM},4001,+,60M,60,0;",
    )
    supp_overlapping = mk_read_multi(
        hdr, "splitI", 95, SEQ60, is_read1=True, is_supplementary=True,
        mate_pos=90, sa_tag=f"{CHROM},91,+,60M,60,0;",
    )
    supp_far = mk_read_multi(
        hdr, "splitI", 4000, SEQ60, is_read1=True, is_supplementary=True,
        mate_pos=90, sa_tag=f"{CHROM},91,+,60M,60,0;",
    )
    mate_read = mk_read_multi(hdr, "splitI", 85, SEQ60, is_read2=True, is_reverse=True, mate_pos=90)

    _, out_bam, mates_bam = run_extract_multi(
        tmp_path, [primary, supp_overlapping, supp_far, mate_read]
    )

    overlap = read_names_and_flags(out_bam)
    assert {("splitI", 90), ("splitI", 95), ("splitI", 85)} <= set(overlap.keys())
    assert ("splitI", 4000) not in overlap

    mates = mates_bam_records(mates_bam)
    # Only the non-overlapping SA target (4000) needed fetching; the
    # overlapping one (95) must NOT be written a second time.
    assert set(mates.keys()) == {("splitI", CHROM, 4000)}
    assert mates[("splitI", CHROM, 4000)] & fh.Tag_bits["SPLITX"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
