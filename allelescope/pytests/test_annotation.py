"""Coding effects from a GFF3/GTF and primer overlap (allelescope.annotation)."""
import pysam
import pytest

from allelescope.annotation import Annotation, revcomp, translate, instance_files

# Gene A, + strand, 11..40:  M  D  F  K  W  P  G  Y  H  *
GENE_A = "ATG" "GAT" "TTT" "AAA" "TGG" "CCC" "GGG" "TAC" "CAT" "TAA"
# Gene B, - strand, 101..130 (coding strand): M K D C F W R T V *
GENE_B = "ATG" "AAA" "GAT" "TGC" "TTT" "TGG" "CGT" "ACC" "GTC" "TAG"
# Gene C, + strand, two parts 201..215 and 215..235 (-1 slippage at 215).
GENE_C1 = "ATG" "AAA" "CCC" "GGG" "TTA"      # 201..215
GENE_C2 = "CGTAAACCCGGGTTTGGTTAA"[1:]       # 216..235 (20 bases)


def _genome():
    seq = ["C"] * 300
    seq[10:40] = GENE_A
    seq[100:130] = revcomp(GENE_B)
    seq[200:215] = GENE_C1
    seq[215:235] = GENE_C2
    return "".join(seq)


GFF = """##gff-version 3
v\t.\tgene\t11\t40\t.\t+\t.\tID=gA;Name=geneA
v\t.\tCDS\t11\t40\t.\t+\t0\tID=cA;Parent=gA
v\t.\tgene\t101\t130\t.\t-\t.\tID=gB;Name=geneB
v\t.\tmRNA\t101\t130\t.\t-\t.\tID=mB;Parent=gB
v\t.\tCDS\t101\t130\t.\t-\t0\tID=cB;Parent=mB
v\t.\tCDS\t201\t215\t.\t+\t0\tID=cC;gene=orf1ab;product=polyprotein
v\t.\tCDS\t215\t235\t.\t+\t0\tID=cC;gene=orf1ab;product=polyprotein
v\t.\tCDS\t201\t215\t.\t+\t0\tID=cC1;gene=orf1ab;product=short form
"""


@pytest.fixture
def ann(tmp_path):
    (tmp_path / "ref.fa").write_text(">v\n" + _genome() + "\n")
    pysam.faidx(str(tmp_path / "ref.fa"))
    (tmp_path / "annotation.gff3").write_text(GFF)
    (tmp_path / "primers.bed").write_text(
        "v\t0\t25\tamp1_LEFT\t1\t+\n"
        "v\t29\t50\tamp1_RIGHT\t1\t-\n"
        "v\t99\t110\tamp2_LEFT\t2\n")
    return Annotation(tmp_path / "ref.fa", tmp_path / "annotation.gff3", tmp_path / "primers.bed")


def one(ann, pos, ref, alt, gene="geneA"):
    effects = [e for e in ann.effects("v", pos, ref, alt) if e["gene"] == gene]
    assert len(effects) == 1, ann.effects("v", pos, ref, alt)
    return effects[0]["change"], effects[0]["type"]


def test_translate():
    assert translate("ATG") == "M" and translate("TAA") == "*" and translate("TGA") == "*"
    assert translate("TGA", 4) == "W" and translate("NNN") == "X"


def test_plus_strand_substitutions(ann):
    assert one(ann, 14, "G", "C") == ("D2H", "missense")
    assert one(ann, 16, "T", "C") == ("D2D", "synonymous")
    assert one(ann, 24, "G", "A") == ("W5*", "stop gained")
    assert one(ann, 12, "T", "C") == ("M1T", "start lost")
    assert one(ann, 38, "T", "C") == ("*10Q", "stop lost")
    # Two bases in one codon, written with a shared leading base.
    assert one(ann, 13, "GGA", "GCC") == ("D2P", "missense")   # GAT -> CCT


def test_plus_strand_indels(ann):
    assert one(ann, 19, "TA", "T") == ("K4fs", "frameshift")
    assert one(ann, 19, "TAAA", "T") == ("K4del1", "in-frame deletion")
    assert one(ann, 22, "A", "AGGG") == ("K4ins1", "in-frame insertion")
    assert one(ann, 9, "CCAT", "C") == ("M1?", "CDS boundary")


def test_minus_strand(ann):
    # Coding offset 3 (codon 2, AAA) is genomic 127, where the + strand has T.
    assert _genome()[126] == "T"
    assert one(ann, 127, "T", "G", "geneB") == ("K2Q", "missense")
    # Gene name found through mRNA -> gene.
    assert ann.cds["v"][1].name == "geneB"


def test_slippage_and_duplicate_names(ann):
    names = sorted(c.name for c in ann.cds["v"])
    assert names == ["geneA", "geneB", "orf1ab [polyprotein]", "orf1ab [short form]"]
    long = next(c for c in ann.cds["v"] if c.name == "orf1ab [polyprotein]")
    assert long.cds_offsets(215) == [14, 15]
    assert len(long.seq) == 15 + 21
    # A change at the slippage base touches a codon in each frame.
    e = [x for x in ann.effects("v", 215, "A", "G") if x["gene"] == "orf1ab [polyprotein]"][0]
    assert e["change"].count(",") == 1


def test_non_coding_and_symbolic(ann):
    assert ann.effects("v", 60, "C", "T") == []
    assert ann.effects("v", 14, "G", "<DEL>") == []
    assert ann.effects("v", 14, "G", "*") == []


def test_allele_codes(ann):
    assert [e["change"] for e in ann.allele_effects("v", 14, "MC") if e["gene"] == "geneA"] == ["D2H"]
    assert [e["type"] for e in ann.allele_effects("v", 20, "D1") if e["gene"] == "geneA"] == ["frameshift"]
    assert [e["type"] for e in ann.allele_effects("v", 23, "IGGG") if e["gene"] == "geneA"] == ["in-frame insertion"]
    assert ann.allele_effects("v", 14, "RG") == []


def test_gtf(tmp_path):
    (tmp_path / "ref.fa").write_text(">v\n" + _genome() + "\n")
    pysam.faidx(str(tmp_path / "ref.fa"))
    (tmp_path / "annotation.gtf").write_text(
        'v\tx\tCDS\t11\t40\t.\t+\t0\tgene_id "gA"; transcript_id "tA"; gene_name "geneA";\n')
    a = Annotation(tmp_path / "ref.fa", tmp_path / "annotation.gtf")
    assert a.effects("v", 14, "G", "C")[0]["change"] == "D2H"
    assert a.effects("v", 14, "G", "C")[0]["gene"] == "geneA"


def test_primers(ann):
    assert [p[2] for p in ann.primers_at("v", 20)] == ["amp1_LEFT"]
    assert [p[2] for p in ann.primers_at("v", 105)] == ["amp2_LEFT"]
    assert ann.primers_at("v", 70) == []
    # Strand from the name when the BED has no strand column.
    assert ann.primers["v"][2][3] == "+"
    # At 20, inside amp1_LEFT (1..25): reads starting at 1..20 carry primer bases.
    assert ann.primer_reads("v", 20, {"1": 5, "15": 2, "26": 3}, {"120": 9}) == 7
    # At 40, inside amp1_RIGHT (30..50): reads ending at 40..53 do.
    assert ann.primer_reads("v", 40, {"1": 5}, {"39": 1, "45": 4, "53": 2, "60": 8}) == 6
    assert ann.primer_reads("v", 70, {"1": 5}, {"80": 1}) == 0


def test_info_and_files(tmp_path, ann):
    info = ann.info()
    assert info["cds"] == 4 and info["primers"] == 3 and info["warnings"] == []
    (tmp_path / "other.fa").write_text(">x\nACGT\n")
    pysam.faidx(str(tmp_path / "other.fa"))
    bad = Annotation(tmp_path / "other.fa", tmp_path / "annotation.gff3")
    assert bad.info()["cds"] == 0
    assert "no sequence names match" in bad.info()["warnings"][0]
    fasta, gff, bed = instance_files(tmp_path)
    assert fasta.name in ("ref.fa", "other.fa") and gff.name == "annotation.gff3" and bed.name == "primers.bed"


def test_server_annotate_extract_and_vcf(tmp_path, monkeypatch):
    from allelescope import file_server
    inst = tmp_path / "inst"
    inst.mkdir()
    (inst / "ref.fa").write_text(">v\n" + _genome() + "\n")
    pysam.faidx(str(inst / "ref.fa"))
    monkeypatch.setattr(file_server, "DATA_DIR", tmp_path)
    client = file_server.app.test_client()
    body = {"variants": [{"chrom": "v", "pos": 14, "ref": "G", "alt": "C,T"},
                         {"chrom": "v", "pos": 60, "ref": "C", "alt": "T"}]}

    # No annotation or primers yet.
    assert client.post("/inst/api/annotate", json=body).get_json()["info"] is None

    # Files added to an existing instance are picked up without a restart.
    (inst / "annotation.gff3").write_text(GFF)
    (inst / "primers.bed").write_text("v\t0\t25\tamp1_LEFT\t1\t+\n")
    d = client.post("/inst/api/annotate", json=body).get_json()
    assert d["info"]["cds"] == 4 and d["info"]["primers"] == 1
    v14 = d["variants"]["v:14"]
    assert [(e["alt"], e["change"]) for e in v14["effects"] if e["gene"] == "geneA"] == [("C", "D2H"), ("T", "D2Y")]
    assert v14["primers"] == ["amp1_LEFT"] and v14["coding"] is True
    assert d["variants"]["v:60"]["effects"] == [] and "primers" not in d["variants"]["v:60"]

    # The reviewed VCF carries EFFECT and PRIMER.
    vcf = inst / "calls.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=v,length=300>\n"
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "v\t14\t.\tG\tC\t.\tPASS\t.\nv\t60\t.\tC\tT\t.\tPASS\t.\n")
    pysam.tabix_index(str(vcf), preset="vcf")
    text = client.post("/inst/api/review_vcf", json={"reviews": {}}).data.decode()
    rec14 = next(l for l in text.splitlines() if l.startswith("v\t14\t"))
    rec60 = next(l for l in text.splitlines() if l.startswith("v\t60\t"))
    assert "EFFECT=geneA:D2H:missense" in rec14 and "PRIMER=amp1_LEFT" in rec14
    assert "EFFECT" not in rec60 and "PRIMER" not in rec60
