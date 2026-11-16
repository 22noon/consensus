"""Per-allele strand counts and mean qualities in the allele table."""
import pysam

from allelescope.extract_mates import _base_quality_at
from allelescope.process_variant import process_variant, site_stem


def _segment(header, start, cigar, quals):
    a = pysam.AlignedSegment(header)
    a.query_name = "r"
    a.reference_id = 0
    a.reference_start = start
    a.cigarstring = cigar
    n = a.infer_query_length()
    a.query_sequence = "A" * n
    a.query_qualities = quals[:n]
    return a


def test_base_quality_at():
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "c", "LN": 1000}]})
    quals = list(range(10, 40))
    r = _segment(header, 100, "3S5M2I5M4D5M", quals)
    assert _base_quality_at(r, 100) == quals[3]          # after the soft clip
    assert _base_quality_at(r, 105) == quals[3 + 5 + 2]  # after the insertion
    assert _base_quality_at(r, 111) is None              # inside the deletion
    assert _base_quality_at(r, 114) == quals[3 + 5 + 2 + 5]
    assert _base_quality_at(r, 99) is None               # before the read


def test_allele_table_has_strand_and_quality(tmp_path):
    seq = "ACGT" * 250
    (tmp_path / "ref.fa").write_text(f">c\n{seq}\n", encoding="utf-8")
    pysam.faidx(str(tmp_path / "ref.fa"))
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "c", "LN": len(seq)}]})
    site0 = 500
    alt = "T" if seq[site0] != "T" else "G"
    reads = []
    # 12 REF fragments, crossing read on alternating strands, MAPQ 60, BQ 30;
    # 8 ALT fragments, crossing read always reverse, MAPQ 40, BQ 20.
    for i in range(20):
        is_alt = i >= 12
        cross_rev = True if is_alt else bool(i % 2)
        start = site0 - 20 - i
        mstart = start + 300
        for crossing in (True, False):
            a = pysam.AlignedSegment(header)
            a.query_name = f"f{i}"
            a.reference_id = 0
            s = start if crossing else mstart
            a.reference_start = s
            a.cigarstring = "100M"
            bases = list(seq[s:s + 100])
            if crossing and is_alt:
                bases[site0 - s] = alt
            a.query_sequence = "".join(bases)
            q = 20 if is_alt else 30
            a.query_qualities = pysam.qualitystring_to_array(chr(33 + q) * 100)
            a.mapping_quality = 40 if is_alt else 60
            rev = cross_rev if crossing else not cross_rev
            first = crossing
            a.flag = 0x1 | 0x2 | (0x40 if first else 0x80) | (0x10 if rev else 0) | (0x0 if rev else 0x20)
            a.next_reference_id = 0
            a.next_reference_start = mstart if crossing else start
            reads.append(a)
    bam = tmp_path / "alignment.markdup.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in sorted(reads, key=lambda r: r.reference_start):
            out.write(r)
    pysam.index(str(bam))

    assert process_variant(tmp_path, "c", site0 + 1, "ref.fa") == 0
    stats = tmp_path / "bam" / f"{site_stem('c', site0 + 1)}.bam.indelstats.txt"
    rows = {f[0]: f for f in (line.split("\t") for line in stats.read_text(encoding="utf-8").splitlines())}

    ref_row, alt_row = rows["R" + seq[site0]], rows["M" + alt]
    assert ref_row[1] == "12" and alt_row[1] == "8"
    assert (int(ref_row[4]), int(ref_row[5])) == (6, 6)
    assert (int(alt_row[4]), int(alt_row[5])) == (0, 8)
    assert float(ref_row[6]) == 30.0 and float(alt_row[6]) == 20.0
    assert float(ref_row[7]) == 60.0 and float(alt_row[7]) == 40.0
