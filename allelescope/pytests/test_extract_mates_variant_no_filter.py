"""
Variant test: filter REMOVED ENTIRELY.

    # deleted:
    # if read.is_supplementary or read.is_secondary:
    #     continue

Under this policy, ANY alignment (primary, supplementary, OR secondary)
flagged is_read1/is_read2 counts as "this mate crosses the region" for
OVERLAP-pair determination.

This monkeypatches extract_mates._counts_as_mate_evidence so you don't
have to hand-edit extract_mates.py to try it -- run this file as-is to
see how the "no filter" policy behaves, independent of your source edits.

Includes:
  - the full baseline suite (should behave identically to the default
    policy for every existing scenario -- none of them depend on
    supplementary/secondary exclusion to reach their expected result)
  - new differential cases that specifically exercise supplementary- and
    secondary-only crossings, where this policy's behavior diverges from
    the default (and from the "secondary-only" variant)
"""
import sys
import os

import pysam
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import allelescope.extract_mates as em
import fake_helpers as fh

from test_extract_mates import (  # reuse builders + full baseline suite
    CHROM,
    mk_read,
    run_extract,
    read_names_and_flags,
    SEQ60,
    test_both_mates_crossing_get_overlap,          # noqa: F401  (re-collected by pytest)
    test_split_primary_crosses_mate_and_supplementary_recovered,  # noqa: F401
    test_duplication_both_segments_cross_deterministic_allele,    # noqa: F401
    test_unpaired_read_never_gets_overlap,          # noqa: F401
    test_missing_mate_reported_not_written,         # noqa: F401
    test_mate_at_window_boundary_is_found,          # noqa: F401
)


@pytest.fixture(autouse=True)
def patch_helpers_and_policy(monkeypatch):
    monkeypatch.setattr(em, "Tag_bits", fh.Tag_bits, raising=False)
    monkeypatch.setattr(em, "tag_read", fh.tag_read, raising=False)
    monkeypatch.setattr(em, "calculate_softclip_lengths", fh.calculate_softclip_lengths, raising=False)
    monkeypatch.setattr(em, "parse_sa_tag", fh.parse_sa_tag, raising=False)
    monkeypatch.setattr(em, "get_allele_at_base", fh.get_allele_at_base, raising=False)
    monkeypatch.setattr(em, "get_indels_at_base", fh.get_indels_at_base, raising=False)
    monkeypatch.setattr(em, "summarise_indels", fh.summarise_indels, raising=False)

    # The policy under test: nothing is excluded.
    monkeypatch.setattr(em, "_counts_as_mate_evidence", lambda read: True, raising=False)


def _hdr():
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": 2_000_000}],
    })


# ---------------------------------------------------------------------------
# Differential case A: BOTH mates' primaries miss the region, but each has a
# SUPPLEMENTARY alignment that crosses it.
# Default policy: no OVERLAP (supplementary excluded from mate-evidence).
# This "no filter" policy: OVERLAP IS set on the two crossing supplementaries.
# ---------------------------------------------------------------------------
def test_supplementary_only_crossing_sets_overlap_under_no_filter(tmp_path):
    hdr = _hdr()
    seq = SEQ60

    r1_primary = mk_read(hdr, "suppOnlyG", 5000, seq, is_read1=True,
                          mate_pos=6000,
                          sa_tag=f"{CHROM},91,+,60M,60,0;")
    r1_supp = mk_read(hdr, "suppOnlyG", 90, seq, is_read1=True, is_supplementary=True,
                       mate_pos=6000,
                       sa_tag=f"{CHROM},5001,+,60M,60,0;")

    r2_primary = mk_read(hdr, "suppOnlyG", 6000, seq, is_read2=True,
                          mate_pos=5000,
                          sa_tag=f"{CHROM},81,+,60M,60,0;")
    r2_supp = mk_read(hdr, "suppOnlyG", 80, seq, is_read2=True, is_supplementary=True,
                       mate_pos=5000,
                       sa_tag=f"{CHROM},6001,+,60M,60,0;")

    allele_counts, out_bam, mates_bam = run_extract(
        tmp_path, [r1_primary, r1_supp, r2_primary, r2_supp]
    )

    overlap = read_names_and_flags(out_bam)
    # Only the two supplementaries (90, 80) physically cross [100, 101);
    # the primaries (5000, 6000) don't overlap the fetch region at all and
    # are never even written to out_bam.
    assert set(overlap.keys()) == {("suppOnlyG", 90), ("suppOnlyG", 80)}

    for xo in overlap.values():
        assert xo & fh.Tag_bits["OVERLAP"], (
            "expected OVERLAP under the 'no filter' policy, since a "
            "supplementary from each mate crosses the region"
        )


# ---------------------------------------------------------------------------
# Differential case B: BOTH mates' primaries miss the region, but each has a
# SECONDARY (multi-mapping) alignment that crosses it.
# Under "no filter", secondary counts too -> OVERLAP set. (Contrast with the
# "secondary-only-excluded" variant file, where this must NOT set OVERLAP.)
# ---------------------------------------------------------------------------
def test_secondary_only_crossing_sets_overlap_under_no_filter(tmp_path):
    hdr = _hdr()
    seq = SEQ60

    r1_primary = mk_read(hdr, "secOnlyH", 5000, seq, is_read1=True, mate_pos=6000)
    r1_secondary = mk_read(hdr, "secOnlyH", 90, seq, is_read1=True, mate_pos=6000)
    r1_secondary.is_secondary = True

    r2_primary = mk_read(hdr, "secOnlyH", 6000, seq, is_read2=True, mate_pos=5000)
    r2_secondary = mk_read(hdr, "secOnlyH", 80, seq, is_read2=True, mate_pos=5000)
    r2_secondary.is_secondary = True

    allele_counts, out_bam, mates_bam = run_extract(
        tmp_path, [r1_primary, r1_secondary, r2_primary, r2_secondary]
    )

    overlap = read_names_and_flags(out_bam)
    assert set(overlap.keys()) == {("secOnlyH", 90), ("secOnlyH", 80)}
    for xo in overlap.values():
        assert xo & fh.Tag_bits["OVERLAP"], (
            "expected OVERLAP under the 'no filter' policy, since a "
            "secondary alignment from each mate crosses the region"
        )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
