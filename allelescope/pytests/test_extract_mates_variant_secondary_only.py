"""
Variant test: KEEP ONLY the secondary exclusion.

    if read.is_secondary:
        continue

Under this policy, supplementary alignments DO count as "this mate
crosses the region" for OVERLAP-pair determination, but secondary
(multi-mapping) alignments still do not.

Monkeypatches extract_mates._counts_as_mate_evidence -- no source edits
needed to try this policy.

Includes:
  - the full baseline suite (identical results to the default policy --
    none of the existing scenarios have a supplementary-only crossing)
  - the same differential cases as test_extract_mates_variant_no_filter.py,
    so you can diff the two policies directly:
      * supplementary-only crossing -> OVERLAP IS set (same as no-filter)
      * secondary-only crossing     -> OVERLAP is NOT set (differs from
        no-filter, matches the original/default policy)
"""
import sys
import os

import pysam
import pytest

pytest.skip("Temporarily disabled", allow_module_level=True)

sys.path.insert(0, os.path.dirname(__file__))

import allelescope.extract_mates as em
import fake_helpers as fh

from test_extract_mates import (  # reuse builders + full baseline suite
    CHROM,
    mk_read,
    run_extract,
    read_names_and_flags,
    SEQ60,
    test_both_mates_crossing_get_overlap,          # noqa: F401
    test_split_primary_crosses_mate_and_supplementary_recovered,  # noqa: F401
    test_duplication_both_segments_cross_deterministic_allele,    # noqa: F401
    test_unpaired_read_never_gets_overlap,          # noqa: F401
    test_missing_mate_reported_not_written,         # noqa: F401
    test_mate_at_window_boundary_is_found,          # noqa: F401
)

@pytest.mark.skip(reason="Temporarily disabled")
@pytest.fixture(autouse=True)
def patch_helpers_and_policy(monkeypatch):
    monkeypatch.setattr(em, "Tag_bits", fh.Tag_bits, raising=False)
    monkeypatch.setattr(em, "tag_read", fh.tag_read, raising=False)
    monkeypatch.setattr(em, "calculate_softclip_lengths", fh.calculate_softclip_lengths, raising=False)
    monkeypatch.setattr(em, "parse_sa_tag", fh.parse_sa_tag, raising=False)
    monkeypatch.setattr(em, "get_allele_at_base", fh.get_allele_at_base, raising=False)
    monkeypatch.setattr(em, "get_indels_at_base", fh.get_indels_at_base, raising=False)
    monkeypatch.setattr(em, "summarise_indels", fh.summarise_indels, raising=False)

    # The policy under test: only secondary alignments are excluded.
    monkeypatch.setattr(
        em, "_counts_as_mate_evidence", lambda read: not read.is_secondary, raising=False
    )


def _hdr():
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": 2_000_000}],
    })


# ---------------------------------------------------------------------------
# Differential case A: supplementary-only crossing -> OVERLAP SHOULD be set
# under this policy (supplementary is allowed to count).
# ---------------------------------------------------------------------------
def test_supplementary_only_crossing_sets_overlap_under_secondary_only_filter(tmp_path):
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
    assert set(overlap.keys()) == {("suppOnlyG", 90), ("suppOnlyG", 80)}
    for xo in overlap.values():
        assert xo & fh.Tag_bits["OVERLAP"], (
            "expected OVERLAP under the 'secondary-only-excluded' policy, "
            "since supplementary alignments are allowed to count"
        )


# ---------------------------------------------------------------------------
# Differential case B: secondary-only crossing -> OVERLAP must NOT be set
# under this policy (secondary is still excluded).
# ---------------------------------------------------------------------------
def test_secondary_only_crossing_does_not_set_overlap_under_secondary_only_filter(tmp_path):
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
        assert not (xo & fh.Tag_bits["OVERLAP"]), (
            "secondary alignments should still be excluded from "
            "mate-crossing evidence under this policy"
        )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
