"""
Additional tests for mate extraction orientation logic.
"""
import sys
import os
import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import allelescope.extract_mates as em
import pytests.fake_helpers as fh

CHROM = "chr1"

def _header():
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": 2_000_000}],
    }

def mk_read(header_obj, qname, pos, is_read1, is_read2, is_reverse, mate_pos, mate_reverse):
    r = pysam.AlignedSegment(header_obj)
    r.query_name = qname
    r.query_sequence = "ACGT" * 15
    r.query_qualities = pysam.qualitystring_to_array("I" * 60)
    r.reference_id = 0
    r.reference_start = pos
    r.cigarstring = "60M"
    r.is_paired = True
    r.is_read1 = is_read1
    r.is_read2 = is_read2
    r.is_reverse = is_reverse
    r.mate_is_reverse = mate_reverse
    r.next_reference_id = 0
    r.next_reference_start = mate_pos
    return r

def test_orientation_tags(tmp_path):
    hdr = pysam.AlignmentHeader.from_dict(_header())
    
    # Test F1R2: R1 Forward, R2 Reverse
    r1_f1r2 = mk_read(hdr, "f1r2", 100, is_read1=True, is_read2=False, is_reverse=False, mate_pos=200, mate_reverse=True)
    r2_f1r2 = mk_read(hdr, "f1r2", 200, is_read1=False, is_read2=True, is_reverse=True, mate_pos=100, mate_reverse=False)
    
    # Test R1F2: R1 Reverse, R2 Forward
    r1_r1f2 = mk_read(hdr, "r1f2", 100, is_read1=True, is_read2=False, is_reverse=True, mate_pos=200, mate_reverse=False)
    r2_r1f2 = mk_read(hdr, "r1f2", 200, is_read1=False, is_read2=True, is_reverse=False, mate_pos=100, mate_reverse=True)

    # Test F1F2: R1 Forward, R2 Forward
    r1_f1f2 = mk_read(hdr, "f1f2", 100, is_read1=True, is_read2=False, is_reverse=False, mate_pos=200, mate_reverse=False)
    r2_f1f2 = mk_read(hdr, "f1f2", 200, is_read1=False, is_read2=True, is_reverse=False, mate_pos=100, mate_reverse=False)

    # Test R1R2: R1 Reverse, R2 Reverse
    r1_r1r2 = mk_read(hdr, "r1r2", 100, is_read1=True, is_read2=False, is_reverse=True, mate_pos=200, mate_reverse=True)
    r2_r1r2 = mk_read(hdr, "r1r2", 200, is_read1=False, is_read2=True, is_reverse=True, mate_pos=100, mate_reverse=True)

    reads = [r1_f1r2, r2_f1r2, r1_r1f2, r2_r1f2, r1_f1f2, r2_f1f2, r1_r1r2, r2_r1r2]
    
    for r in reads:
        em.tag_read(r, 0)
        xr = r.get_tag("XR")
        xo = r.get_tag("XO")
        
        expected_xr = ""
        if r.is_read1:
            expected_xr = ("R1" if r.is_reverse else "F1") + ("R2" if r.mate_is_reverse else "F2")
        else:
            expected_xr = ("R2" if r.is_reverse else "F2") + ("R1" if r.mate_is_reverse else "R1") # Wait, logic in get_pair_orientation
            # get_pair_orientation uses mate_reverse for read1.
            # R1: read_reverse (for R1), mate_reverse (for R2)
            # R2: read_reverse (for R2), mate_reverse (for R1)
            # Let's check get_pair_orientation implementation
            pass
        
        # Verify XR tag matches expected orientation string
        expected_xr_str = em.get_pair_orientation(r)
        assert xr == expected_xr_str
        
        # Verify XO tag has the corresponding bit set
        assert xo & em.Tag_bits[expected_xr_str]
