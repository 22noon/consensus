"""Read evidence per allele: position in read, other edits, mapping, clipping,
insert size, orientation and duplicates (extract_mates._read_evidence and the
<bam>.evidence.json summary)."""
import json

import pysam

from allelescope.extract_mates import (_site_query_pos, _allele_edit_length, _read_evidence,
                                       allele_evidence_summary, NEAR_END_BP)
from allelescope.process_variant import process_variant, site_stem


HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "c", "LN": 1000}]})


def _segment(start, cigar, flag=0, mapq=60, tags=()):
    a = pysam.AlignedSegment(HEADER)
    a.query_name = "r"
    a.reference_id = 0
    a.reference_start = start
    a.cigarstring = cigar
    n = a.infer_query_length()
    a.query_sequence = "A" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.flag = flag
    a.mapping_quality = mapq
    for tag, value in tags:
        a.set_tag(tag, value)
    return a


def test_site_query_pos():
    r = _segment(100, "3S5M2I5M4D5M")
    assert _site_query_pos(r, 100)[0] == 3
    assert _site_query_pos(r, 105)[0] == 3 + 5 + 2
    assert _site_query_pos(r, 111) == (3 + 5 + 2 + 5, None)   # in the deletion
    assert _site_query_pos(r, 99) == (None, None)


def test_allele_edit_length():
    assert _allele_edit_length("RA") == 0
    assert _allele_edit_length("MT") == 1
    assert _allele_edit_length("D7") == 7
    assert _allele_edit_length("IATT") == 3


def test_read_evidence_fields():
    # 5 bp soft clip on the left, site 2 bp into the aligned part.
    r = _segment(100, "5S95M", flag=0x1 | 0x40 | 0x10, mapq=0,
                 tags=[("NM", 4), ("XA", "c,+500,100M,1;")])
    r.next_reference_id = 0
    r.next_reference_start = 300
    r.template_length = -350
    (rev, mq, bq, end_dist, other, ambiguous, chimeric, clipped, isize, f1r2, dup) = \
        _read_evidence(r, 102, "MT")
    assert rev and mq == 0 and bq == 40
    assert end_dist == 2
    assert other == 3                     # NM 4 minus the SNV itself
    assert ambiguous and not chimeric and clipped
    assert isize == 350
    assert f1r2 is False                  # read 1 on the reverse strand: F2R1
    assert not dup

    # Right-hand soft clip only: 100M5S, site well away from it.
    r = _segment(100, "100M5S")
    assert _read_evidence(r, 120, "RA")[7] is False
    r = _segment(100, "100M5S")
    assert _read_evidence(r, 195, "RA")[7] is True
    # bwa's integer XS equal to AS: mapped as well elsewhere.
    r = _segment(100, "100M", tags=[("AS", 100), ("XS", 100)])
    assert _read_evidence(r, 150, "RA")[5] is True
    # This module's own XS is a string and must not count.
    r = _segment(100, "100M", tags=[("AS", 100), ("XS", "")])
    assert _read_evidence(r, 150, "RA")[5] is False


def _write_site_bam(tmp_path, seq, site0):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "c", "LN": len(seq)}]})
    alt = "T" if seq[site0] != "T" else "G"
    reads = []
    # 20 REF fragments: site mid-read, clean, insert 400, orientations balanced.
    # 10 ALT fragments: site 3 bp from the read end, 3 other edits, half with
    # bwa XA hits, insert 250, all F1R2, two marked duplicate.
    for i in range(30):
        is_alt = i >= 20
        if is_alt:
            start = site0 - 96            # site at read offset 96: 3 bp from the end
            insert = 250
            cross_rev = False             # read 1 forward: F1R2
        else:
            start = site0 - 50 + (i % 5)
            insert = 400
            cross_rev = bool(i % 2)
        mstart = start + insert - 100
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
            a.query_qualities = pysam.qualitystring_to_array("?" * 100)
            a.mapping_quality = 60
            rev = cross_rev if crossing else not cross_rev
            a.flag = (0x1 | 0x2 | (0x40 if crossing else 0x80) | (0x10 if rev else 0)
                      | (0x20 if not rev else 0) | (0x400 if is_alt and i >= 28 else 0))
            a.next_reference_id = 0
            a.next_reference_start = mstart if crossing else start
            a.template_length = insert if s == start else -insert
            a.set_tag("NM", 4 if (crossing and is_alt) else 0)
            if crossing and is_alt and i % 2:
                a.set_tag("XA", "c,+10,100M,2;")
            reads.append(a)
    bam = tmp_path / "alignment.markdup.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in sorted(reads, key=lambda r: r.reference_start):
            out.write(r)
    pysam.index(str(bam))
    return alt


def test_evidence_json_written(tmp_path):
    seq = "ACGT" * 250
    (tmp_path / "ref.fa").write_text(f">c\n{seq}\n", encoding="utf-8")
    pysam.faidx(str(tmp_path / "ref.fa"))
    site0 = 500
    alt = _write_site_bam(tmp_path, seq, site0)

    assert process_variant(tmp_path, "c", site0 + 1, "ref.fa") == 0
    path = tmp_path / "bam" / f"{site_stem('c', site0 + 1)}.bam.evidence.json"
    ev = json.loads(path.read_text(encoding="utf-8"))
    ref, a = ev["R" + seq[site0]], ev["M" + alt]

    assert ref["n"] == 20 and a["n"] == 10
    assert a["end_median"] == 3 and a["near_end"] == 10
    assert ref["near_end"] == 0 and ref["end_median"] >= NEAR_END_BP
    assert a["nm_mean"] == 3 and ref["nm_mean"] == 0
    assert a["ambiguous"] == 5 and ref["ambiguous"] == 0
    assert a["isize_median"] == 250 and ref["isize_median"] == 400
    assert (a["f1r2"], a["f2r1"]) == (10, 0)
    assert ref["f1r2"] == 10 and ref["f2r1"] == 10
    assert a["dup"] == 2 and ref["dup"] == 0
    assert a["chimeric"] == 0 and a["clipped"] == 0


def test_summary_handles_missing_values():
    from allelescope.extract_mates import _new_allele_stats
    st = _new_allele_stats()
    st["fwd"] = 1
    s = allele_evidence_summary({"RA": st})["RA"]
    assert s["n"] == 1 and s["end_median"] is None and s["nm_mean"] is None
    assert s["isize_median"] is None


def test_review_vcf_carries_read_flags(tmp_path, monkeypatch):
    from allelescope import file_server
    inst = tmp_path / "inst"
    inst.mkdir()
    vcf = inst / "calls.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=c,length=1000>\n"
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "c\t10\t.\tA\tT\t.\tPASS\t.\nc\t20\t.\tC\tG\t.\tPASS\t.\n")
    pysam.tabix_index(str(vcf), preset="vcf")
    monkeypatch.setattr(file_server, "DATA_DIR", tmp_path)
    client = file_server.app.test_client()
    r = client.post("/inst/api/review_vcf", json={
        "reviews": {"c:10": {"s": "artefact"}},
        "flags": {"c:10": ["read-end", "mismatches", "bad flag;x"], "c:20": []},
    })
    assert r.status_code == 200
    text = r.data.decode()
    assert "##INFO=<ID=READ_FLAGS" in text
    rec10 = next(l for l in text.splitlines() if l.startswith("c\t10\t"))
    rec20 = next(l for l in text.splitlines() if l.startswith("c\t20\t"))
    assert "READ_FLAGS=read-end,mismatches" in rec10 and "REVIEW=artefact" in rec10
    assert "bad" not in rec10                      # malformed names are dropped
    assert "READ_FLAGS" not in rec20
