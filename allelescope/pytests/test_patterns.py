"""Read patterns at a site (allelescope.patterns) and the pattern read filter."""
import random

import pysam
import pytest

import allelescope.patterns as P
from allelescope.patterns import Settings, analyse, cached_analysis, keep_names
from allelescope.process_variant import process_variant, site_stem

SITE = 400
RL = 100


def _other(b, k=0):
    return [x for x in "ACGT" if x != b][k]


@pytest.fixture
def site(tmp_path):
    """
    Pairs crossing SITE: 30 reference, 14 lineage B (SITE, 340, 470),
    8 lineage C (SITE, 430), 6 artefact (SITE + 403, only 2-5 bp from the end
    of reverse reads), 4 noisy (lineage C with scattered mismatches),
    3 chimeras (B left of 420, C right), 1 rare listed variant (at 360), and
    2 pairs whose mate maps to another sequence.
    """
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(1000))
    other = "".join(rng.choice("ACGT") for _ in range(1000))
    (tmp_path / "ref.fa").write_text(f">c\n{seq}\n>d\n{other}\n")
    pysam.faidx(str(tmp_path / "ref.fa"))
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "c", "LN": 1000}, {"SN": "d", "LN": 1000}]})
    alt = {p: _other(seq[p - 1]) for p in (SITE, 340, 470, 430, 403, 360)}
    B = {SITE: alt[SITE], 340: alt[340], 470: alt[470]}
    C = {SITE: alt[SITE], 430: alt[430]}
    reads = []

    def read(name, start0, hap, first, rev, scatter=(), ref_id=0, mate_ref=0, mate_start=0, tlen=0):
        s = seq if ref_id == 0 else other
        bases = [hap.get(start0 + i + 1, s[start0 + i]) if ref_id == 0 else s[start0 + i] for i in range(RL)]
        for p in scatter:
            if start0 < p <= start0 + RL:
                bases[p - start0 - 1] = _other(bases[p - start0 - 1], 1)
        a = pysam.AlignedSegment(header)
        a.query_name, a.reference_id, a.reference_start = name, ref_id, start0
        a.cigarstring = f"{RL}M"
        a.query_sequence = "".join(bases)
        a.query_qualities = pysam.qualitystring_to_array("?" * RL)
        a.mapping_quality = 60
        a.flag = 0x1 | (0x40 if first else 0x80) | (0x10 if rev else 0x20) | (0x2 if ref_id == mate_ref else 0)
        a.next_reference_id, a.next_reference_start, a.template_length = mate_ref, mate_start, tlen
        reads.append(a)

    def pair(name, start, insert, hap, scatter=(), chimera=None, end_at=None):
        if end_at is not None:
            start = end_at - insert
        h1 = h2 = hap
        if chimera:
            h1 = h2 = {p: b for p, b in chimera[0].items() if p < 420} | {p: b for p, b in chimera[1].items() if p >= 420}
        m = start + insert - RL
        read(name, start, h1, True, False, scatter, mate_start=m, tlen=insert)
        read(name, m, h2, False, True, scatter, mate_start=start, tlen=-insert)

    for i in range(30):
        pair(f"ref{i}", SITE - 60 + i % 20, 230, {})
    for i in range(14):
        pair(f"B{i}", 300 + i % 3, 240, B)            # reaches 340 .. 470+
    for i in range(8):
        pair(f"C{i}", 335 + i % 3, 210, C)          # read 1 covers 340, SITE and 430
    for i in range(6):                                # reverse read ends 2-5 bp after 403
        pair(f"art{i}", 0, 200, {SITE: alt[SITE], 403: alt[403]}, end_at=405 + i % 3)
    for i in range(4):
        # scattered mismatches, different in each pair and away from the variants
        pair(f"noisy{i}", 335, 210, C,
             scatter=[p for p in range(336 + i, 545, 12) if p not in alt])
    for i in range(3):
        pair(f"chim{i}", 310, 240, None, chimera=(B, C))
    pair("listed0", 320, 220, {360: alt[360]})
    for i in range(2):                                # mate on sequence d
        read(f"split{i}", SITE - 50, {}, True, False, ref_id=0, mate_ref=1, mate_start=500)
        read(f"split{i}", 500, {}, False, True, ref_id=1, mate_ref=0, mate_start=SITE - 50)
    with pysam.AlignmentFile(str(tmp_path / "alignment.markdup.bam"), "wb", header=header) as out:
        for r in sorted(reads, key=lambda r: (r.reference_id, r.reference_start)):
            out.write(r)
    pysam.index(str(tmp_path / "alignment.markdup.bam"))
    assert process_variant(tmp_path, "c", SITE, "ref.fa") == 0
    bam = tmp_path / "bam" / f"{site_stem('c', SITE)}.calmd.bam"
    return {"dir": tmp_path, "bam": bam, "alt": alt, "seq": seq}


def _by_kind(result):
    out = {}
    for name, pid in result["assignment"].items():
        out.setdefault(name.rstrip("0123456789"), {}).setdefault(pid, 0)
        out[name.rstrip("0123456789")][pid] += 1
    return out


def test_patterns(site):
    alt = site["alt"]
    listed = {SITE: {"M" + alt[SITE]}, 360: {"M" + alt[360]}}
    r = analyse(site["bam"], site["dir"] / "ref.fa", "c", SITE, Settings(listed=listed))
    positions = [p["pos"] for p in r["positions"]]
    assert positions == [340, 360, SITE, 403, 430, 470]
    assert r["total"] == 30 + 14 + 8 + 6 + 4 + 3 + 1 + 2
    pats = {p["id"]: p for p in r["patterns"]}
    kinds = _by_kind(r)

    def only(kind):
        (pid, n), = kinds[kind].items()
        return pats[pid]

    b, c, art, noisy = only("B"), only("C"), only("art"), only("noisy")
    assert b["states"][str(340)] == "M" + alt[340] and b["states"][str(470)] == "M" + alt[470]
    assert c["states"][str(430)] == "M" + alt[430] and c["states"][str(340)] == "R"
    assert art["states"][str(403)] == "M" + alt[403]
    assert art["evidence"]["rev"] == 6 and art["evidence"]["fwd"] == 0 and art["evidence"]["near_end"] == 6
    assert noisy["noisy"] and noisy is not c and noisy["evidence"]["other_edits"] >= P.NOISY_EDITS
    chim = only("chim")
    assert chim["mixed"] == [b["id"], c["id"]] and not chim["common"]
    assert only("listed")["states"][str(360)] == "M" + alt[360]      # rare, but listed
    assert b["common"] and art["common"]
    # Mates on sequence d do not add positions or differences.
    assert set(kinds["split"]) <= {"p1", "uc", "ur"} or only("split")["states"][str(SITE)] == "R"


def test_without_listed_rare_variant_is_not_a_position(site):
    r = analyse(site["bam"], site["dir"] / "ref.fa", "c", SITE, Settings())
    assert 360 not in [p["pos"] for p in r["positions"]]
    assert SITE in [p["pos"] for p in r["positions"]]


def test_md_and_full_comparison_agree(site, monkeypatch):
    a = analyse(site["bam"], site["dir"] / "ref.fa", "c", SITE, Settings())
    monkeypatch.setattr(P, "_md_mismatches", lambda read: None)
    b = analyse(site["bam"], site["dir"] / "ref.fa", "c", SITE, Settings())
    assert a == b


def test_keep_names_and_cache(site):
    r = cached_analysis(site["bam"], site["dir"] / "ref.fa", "c", SITE, Settings())
    again = cached_analysis(site["bam"], site["dir"] / "ref.fa", "c", SITE, Settings())
    assert r == again
    assert len(list(site["bam"].parent.glob("*.patterns.*.json"))) == 1
    common = keep_names(r, "common")
    rare = keep_names(r, "rare")
    assert common.isdisjoint(rare) and len(common | rare) == r["total"]
    assert {n for n in rare if n.startswith("chim")} and not {n for n in common if n.startswith("art")} - common
    art_id = r["assignment"]["art0"]
    assert keep_names(r, art_id) == {n for n in r["assignment"] if n.startswith("art")}
    assert keep_names(r, "") is None


def test_settings_from_params():
    s = Settings.from_params({"PatN": "5", "PatPct": "10", "PatBQ": "30", "PatMin": "20",
                              "PatMinN": "2", "PatSites": "400:MA|MC,360:D2"})
    assert (s.min_reads, s.min_fraction, s.min_bq, s.common_fraction, s.common_reads) == (5, 0.1, 30, 0.2, 2)
    assert s.listed == {400: {"MA", "MC"}, 360: {"D2"}}
    d = Settings.from_params({"PatN": "x"})
    assert d.min_reads == 3 and d.listed == {}


def test_differences_of_a_read():
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "c", "LN": 100}]})
    ref = "ACGTACGTACGTACGTACGT"
    a = pysam.AlignedSegment(header)
    a.reference_id, a.reference_start = 0, 0
    a.cigarstring = "4M2D3M1I5M"
    # ref ACGT|AC deleted|GTA|insert G|CGTAC: mismatch at 4, deletion of 5-6, insertion before 10
    a.query_sequence = "ACGA" + "GTA" + "G" + "CGTAC"
    a.query_qualities = pysam.qualitystring_to_array("IIII" + "III" + "I" + "IIIII")
    diffs, low, inside, span = P.read_differences(a, ref, 1, 20)
    assert diffs == {4: "MA", 5: "D2", 10: "IG"}
    assert inside == {6} and span == (1, 14)


def test_server(site, monkeypatch):
    from allelescope import file_server
    monkeypatch.setattr(file_server, "DATA_DIR", site["dir"].parent)
    monkeypatch.setattr(file_server, "CACHE_ROOT", site["dir"].parent / "cache")
    client = file_server.app.test_client()
    inst = site["dir"].name
    d = client.get(f"/{inst}/api/patterns?Chrom=c&Pos={SITE}&Ref=ref.fa").get_json()
    assert d["success"] and "assignment" not in d and d["total"] == 68
    names = lambda q: {r.query_name for r in _bam(client, f"/{inst}/api/bam?Chrom=c&Pos={SITE}&Ref=ref.fa{q}", site["dir"])}
    everything = names("")
    common, rare = names("&Pat=common"), names("&Pat=rare")
    assert common | rare == everything and not common & rare
    assert any(n.startswith("chim") for n in rare) and not any(n.startswith("chim") for n in common)


def _bam(client, url, tmp):
    data = client.get(url).data
    path = tmp / "served.bam"
    path.write_bytes(data)
    with pysam.AlignmentFile(str(path), "rb") as f:
        return list(f)
