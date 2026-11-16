"""Subsampling of deep sites: exact allele counts, whole fragments, stable choice."""
import json
import random

import pysam
import pytest

from allelescope import file_server
from allelescope.extract_mates import _choose_read_names
from allelescope.process_variant import process_variant, read_site_info, site_sample, site_stem

CHROM = "chr1"
SITE = 1000          # 1-based position of the test SNP
N_PAIRS = 300
N_ALT = 90           # pairs whose site-crossing read carries the ALT base


def _make_site(tmp_path):
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(5000))
    (tmp_path / "ref.fa").write_text(f">{CHROM}\n{seq}\n", encoding="utf-8")
    pysam.faidx(str(tmp_path / "ref.fa"))
    ref_base = seq[SITE - 1]
    alt_base = "C" if ref_base != "C" else "G"

    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": CHROM, "LN": len(seq)}]})
    reads = []
    for i in range(N_PAIRS):
        start = SITE - 1 - 10 - (i % 80)      # crosses the site
        mstart = start + 300 + (i % 50)       # mate does not
        for r1, s, ms, rev in ((True, start, mstart, False), (False, mstart, start, True)):
            a = pysam.AlignedSegment(header)
            a.query_name = f"frag{i:04d}"
            a.reference_id = 0
            a.reference_start = s
            a.cigarstring = "100M"
            bases = list(seq[s:s + 100])
            if r1 and i < N_ALT:
                bases[SITE - 1 - s] = alt_base
            a.query_sequence = "".join(bases)
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.mapping_quality = 60
            a.flag = 0x1 | 0x2 | (0x40 if r1 else 0x80) | (0x10 if rev else 0) | (0x0 if rev else 0x20)
            a.next_reference_id = 0
            a.next_reference_start = ms
            reads.append(a)
    bam = tmp_path / "alignment.markdup.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in sorted(reads, key=lambda r: r.reference_start):
            out.write(r)
    pysam.index(str(bam))


def _served(tmp_path, sample=None):
    path = tmp_path / "bam" / f"{site_stem(CHROM, SITE, sample)}.calmd.bam"
    per_name = {}
    with pysam.AlignmentFile(str(path)) as f:
        for r in f.fetch(until_eof=True):
            per_name.setdefault(r.query_name, []).append(r)
    return per_name


def _counts(tmp_path, sample=None):
    stats = tmp_path / "bam" / f"{site_stem(CHROM, SITE, sample)}.bam.indelstats.txt"
    return [line.split("\t")[:2] for line in stats.read_text(encoding="utf-8").splitlines()]


def test_choose_read_names():
    names = {f"r{i}" for i in range(1000)}
    assert _choose_read_names(names, None) is None
    assert _choose_read_names(names, 0) is None
    assert _choose_read_names(names, 1000) is None          # nothing to drop
    kept = _choose_read_names(names, 100)
    assert len(kept) == 100 and kept <= names
    assert _choose_read_names(set(sorted(names, reverse=True)), 100) == kept   # order-independent
    assert _choose_read_names(names, 200) >= kept                               # nested sizes


def test_site_sample_and_stem():
    assert site_sample("2000") == 2000
    assert site_sample(None) is None and site_sample("") is None
    assert site_sample("0") is None and site_sample("-5") is None and site_sample("x") is None
    assert site_stem("chr1", 5) == "variant_chr1_5"
    assert site_stem("chr1", 5, "2000") == "variant_chr1_5.s2000"
    assert site_stem("chr1", 5, "") == "variant_chr1_5"


def test_subsample_keeps_whole_fragments_and_exact_counts(tmp_path):
    _make_site(tmp_path)
    assert process_variant(tmp_path, CHROM, SITE, "ref.fa") == 0
    assert process_variant(tmp_path, CHROM, SITE, "ref.fa", sample=50) == 0

    full, sub = _served(tmp_path), _served(tmp_path, 50)
    assert len(full) == N_PAIRS
    assert len(sub) == 50
    assert set(sub) <= set(full)
    # every kept fragment comes with its mate
    assert all(len(recs) == 2 for recs in sub.values())

    # allele counts use every read, whatever is displayed
    assert _counts(tmp_path, 50) == _counts(tmp_path)
    assert sum(int(c) for _, c in _counts(tmp_path, 50)) == N_PAIRS
    assert dict(_counts(tmp_path)).get(next(xa for xa, c in _counts(tmp_path) if xa.startswith("M"))) == str(N_ALT)

    assert read_site_info(tmp_path, CHROM, SITE, 50) == {"pairs_total": N_PAIRS, "pairs_kept": 50, "sample": 50}
    assert read_site_info(tmp_path, CHROM, SITE)["pairs_kept"] == N_PAIRS


def test_subsample_is_the_same_every_time(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for d in (a, b):
        d.mkdir()
        _make_site(d)
        process_variant(d, CHROM, SITE, "ref.fa", sample=40)
    assert set(_served(a, 40)) == set(_served(b, 40))


def test_sample_larger_than_depth_keeps_everything(tmp_path):
    _make_site(tmp_path)
    process_variant(tmp_path, CHROM, SITE, "ref.fa", sample=10_000)
    assert len(_served(tmp_path, 10_000)) == N_PAIRS
    assert read_site_info(tmp_path, CHROM, SITE, 10_000)["pairs_kept"] == N_PAIRS


@pytest.fixture
def client(tmp_path, monkeypatch):
    inst = tmp_path / "inst"
    inst.mkdir()
    _make_site(inst)
    monkeypatch.setattr(file_server, "DATA_DIR", tmp_path)
    monkeypatch.setattr(file_server, "CACHE_ROOT", tmp_path / "cache")
    return file_server.app.test_client()


def test_depth_endpoint_and_sampled_extract(client):
    depth = client.get(f"/inst/api/depth?chrom={CHROM}&pos={SITE}").get_json()
    assert depth == {"reads": N_PAIRS, "pairs": N_PAIRS, "prepared": {"all": False, "samples": []}}

    res = client.get(f"/inst/api/extract?chrom={CHROM}&pos={SITE}&ref=ref.fa&sample=25").get_json()
    assert res["success"] and res["sampling"]["pairs_kept"] == 25
    assert sum(a["count"] for a in res["alleles"]) == N_PAIRS

    depth = client.get(f"/inst/api/depth?chrom={CHROM}&pos={SITE}").get_json()
    assert depth["prepared"] == {"all": False, "samples": [25]}

    bam = client.get(f"/inst/api/bam?Chrom={CHROM}&Pos={SITE}&Ref=ref.fa&Flagf=0&FlagF=0"
                     f"&Tagf=0&TagF=0&Token=&XAFilter=&Sample=25")
    assert bam.status_code == 200

    assert client.get(f"/inst/api/depth?chrom={CHROM}&pos=abc").status_code == 400


def test_concurrent_requests_for_one_view(client):
    """IGV asks for a BAM and its index together; a site that is not prepared
    yet must be prepared once and every request must succeed."""
    import threading
    q = (f"Chrom={CHROM}&Pos={SITE}&Ref=ref.fa&Flagf=0&FlagF=0&Tagf=0&TagF=0"
         "&Token=&XAFilter=&minMapQ=5")
    codes = []

    def get(kind):
        codes.append(client.get(f"/inst/api/{kind}?{q}").status_code)

    threads = [threading.Thread(target=get, args=(k,)) for k in ["bam", "bai"] * 6]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert codes == [200] * 12
