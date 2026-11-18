"""
Synthetic amplicon data for gene effects and primer sites (user manual,
section 9.4): a 3 kb virus-like genome with one gene on each strand, a tiled
primer scheme of 400 bp amplicons (primers not trimmed from the reads), and
variants of each kind of coding effect, one hidden under a primer and one
artefact that exists only in primer sequence.

Usage: python simulate_amplicons.py OUT_DIR
  -> ref.fa, reads.bam(.bai), calls.vcf, genes.gff3, primers.bed
Build with --bam-processing skip, --annotation genes.gff3 --primers primers.bed.
"""
import random
import sys
from pathlib import Path

import pysam

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(11)
L, RL, PRIMER, AMP, STEP, PAIRS = 3000, 250, 22, 400, 330, 50
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def orf(codons, stop):
    body = []
    while len(body) < codons - 2:
        c = "".join(rng.choice("ACGT") for _ in range(3))
        if c not in STOPS:
            body.append(c)
    return "ATG" + "".join(body) + stop


seq = [rng.choice("ACGT") for _ in range(L)]
gene_a = orf(400, "TAA")                         # + strand, 101..1300
seq[100:1300] = gene_a
gene_b = orf(300, "TAG")                         # - strand, 1501..2400
seq[1500:2400] = rc(gene_b)
seq = "".join(seq)
(out / "ref.fa").write_text(">virus\n" + "\n".join(seq[i:i + 60] for i in range(0, L, 60)) + "\n")
pysam.faidx(str(out / "ref.fa"))
(out / "genes.gff3").write_text(
    "##gff-version 3\n"
    "virus\tsim\tgene\t101\t1300\t.\t+\t.\tID=gene-polA;Name=polA\n"
    "virus\tsim\tCDS\t101\t1300\t.\t+\t0\tID=cds-polA;Parent=gene-polA;gene=polA;product=polymerase A\n"
    "virus\tsim\tgene\t1501\t2400\t.\t-\t.\tID=gene-capB;Name=capB\n"
    "virus\tsim\tCDS\t1501\t2400\t.\t-\t0\tID=cds-capB;Parent=gene-capB;gene=capB;product=capsid B\n")

# Primer scheme: amplicon k covers [s, s + AMP) (0-based), primers 22 bp at each end.
amplicons = []
with open(out / "primers.bed", "w") as bed:
    k, s = 1, 30
    while s + AMP <= L - 30:
        e = s + AMP
        amplicons.append((k, s, e))
        pool = 1 if k % 2 else 2
        bed.write(f"virus\t{s}\t{s + PRIMER}\tsim_{k}_LEFT\t{pool}\t+\n")
        bed.write(f"virus\t{e - PRIMER}\t{e}\tsim_{k}_RIGHT\t{pool}\t-\n")
        k, s = k + 1, s + STEP


def codon_pos(gene_start1, codon, offset=0):
    """1-based genomic position of a base of a codon on a + strand gene."""
    return gene_start1 + 3 * (codon - 1) + offset


def minus_pos(gene_end1, codon, offset=0):
    """1-based genomic position of a codon base on a - strand gene ending at gene_end1."""
    return gene_end1 - (3 * (codon - 1) + offset)


def other(b, avoid=()):
    return next(x for x in "ACGT" if x != b and x not in avoid)


def mutate_codon(cod, idx, want):
    """A base for position idx of codon cod giving an amino acid of kind want."""
    from itertools import product
    aa = translate_codon(cod)
    for b in "ACGT":
        if b == cod[idx]:
            continue
        new = cod[:idx] + b + cod[idx + 1:]
        a2 = translate_codon(new)
        if want == "syn" and a2 == aa:
            return b
        if want == "mis" and a2 not in (aa, "*"):
            return b
        if want == "stop" and a2 == "*":
            return b
    return None


_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: _AA[16 * i + 4 * j + k] for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}


def translate_codon(c):
    return CODE[c]


variants = []   # (pos1, ref, alt, af, kind, mode)  mode: "sample" or "primer-artefact"


def add_snv(pos1, alt, af, kind, mode="sample"):
    variants.append((pos1, seq[pos1 - 1], alt, af, kind, mode))


# 1. polA missense, fixed.  Pick a codon whose 2nd base can give a missense change.
for c in range(60, 90):
    cod = gene_a[3 * (c - 1): 3 * c]
    b = mutate_codon(cod, 1, "mis")
    if b:
        add_snv(codon_pos(101, c, 1), b, 1.0, "missense")
        break
# 2. polA synonymous, 50 %.
for c in range(120, 160):
    cod = gene_a[3 * (c - 1): 3 * c]
    b = mutate_codon(cod, 2, "syn")
    if b:
        add_snv(codon_pos(101, c, 2), b, 0.5, "synonymous")
        break
# 3. polA stop gained, 20 %.
for c in range(200, 260):
    cod = gene_a[3 * (c - 1): 3 * c]
    for i in range(3):
        b = mutate_codon(cod, i, "stop")
        if b:
            break
    if b:
        add_snv(codon_pos(101, c, i), b, 0.2, "stop gained")
        break
# 4. polA frameshift: 1 bp deletion, 60 %. VCF: anchor base + deleted base.
d = codon_pos(101, 300, 1)
variants.append((d - 1, seq[d - 2:d], seq[d - 2], 0.6, "frameshift", "sample"))
# 5. polA missense inside the right primer of amplicon 3, fixed in the sample;
#    amplicon 3's reads end there with primer sequence, amplicon 4's do not.
k3 = amplicons[2]
site = k3[2] - 8                                   # 0-based, 8 bp from the amplicon end
add_snv(site + 1, other(seq[site]), 1.0, "in primer (masked)")
# 6. Artefact only in primer sequence: left primer of amplicon 5, 70 % of its reads.
k5 = amplicons[4]
site = k5[1] + 10
add_snv(site + 1, other(seq[site]), 0.0, "primer artefact", "primer-artefact")
# 7. capB missense on the minus strand, 80 %.
for c in range(40, 80):
    cod = gene_b[3 * (c - 1): 3 * c]
    b = mutate_codon(cod, 0, "mis")
    if b:
        add_snv(minus_pos(2400, c, 0), rc(b), 0.8, "missense (- strand)")
        break
# 8. capB in-frame insertion of 3 bp, 40 %.
ins_at = minus_pos(2400, 150, 2)                   # 1-based anchor
variants.append((ins_at, seq[ins_at - 1], seq[ins_at - 1] + "GCA", 0.4, "in-frame insertion", "sample"))
# 9. Intergenic SNV, 30 %.
add_snv(2650, other(seq[2649]), 0.3, "non-coding")

variants.sort()

header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                          "SQ": [{"SN": "virus", "LN": L}]})


def hap_edits(primer_artefact_on):
    """Edits of one template: {0-based pos: ("S", base) | ("D", n) | ("I", bases)}."""
    edits = {}
    for pos1, ref, alt, af, kind, mode in variants:
        if mode == "sample" and rng.random() < af:
            if len(ref) == len(alt) == 1:
                edits[pos1 - 1] = ("S", alt)
            elif len(ref) > len(alt):
                edits[pos1] = ("D", len(ref) - len(alt))          # deleted bases start after the anchor
            else:
                edits[pos1] = ("I", alt[1:])                      # inserted before pos1 (0-based pos1)
    return edits


def walk(start0, end0, edits, primer_lo, primer_hi, artefacts):
    """Read bases and CIGAR for reference span [start0, end0) of a template.
    Reference positions in [primer_lo, primer_hi) show the primer (reference) base,
    or an artefact base."""
    bases, ops = [], []

    def op(kind, n=1):
        if ops and ops[-1][0] == kind:
            ops[-1][1] += n
        else:
            ops.append([kind, n])

    p = start0
    while p < end0:
        e = edits.get(p)
        in_primer = primer_lo <= p < primer_hi
        if e and e[0] == "I" and not in_primer and p > start0:
            bases.extend(e[1]); op("I", len(e[1]))
        if e and e[0] == "D" and not in_primer:
            n = min(e[1], end0 - p)
            op("D", n); p += n
            continue
        b = seq[p]
        if in_primer:
            b = artefacts.get(p, b)
        elif e and e[0] == "S":
            b = e[1]
        bases.append(b); op("M")
        p += 1
    return "".join(bases), ops


def nm(read_seq, ops, start0):
    n, q, r = 0, 0, start0
    for kind, length in ops:
        if kind == "M":
            n += sum(1 for i in range(length) if read_seq[q + i] != seq[r + i])
            q += length; r += length
        elif kind == "I":
            n += length; q += length
        else:
            n += length; r += length
    return n


reads = []
frag = 0
art = {pos1 - 1: alt for pos1, ref, alt, af, kind, mode in variants if mode == "primer-artefact"}
for k, s, e in amplicons:
    for i in range(PAIRS):
        frag += 1
        edits = hap_edits(False)
        artefacts = {p: b for p, b in art.items() if rng.random() < 0.7}
        name = f"amp{k}_{i}"
        # read 1: forward from the left primer; read 2: reverse, ending at the right primer.
        r1 = walk(s, s + RL, edits, s, s + PRIMER, artefacts)
        r2 = walk(e - RL, e, edits, e - PRIMER, e, artefacts)
        for which, (rs, ops), start, rev in ((1, r1, s, False), (2, r2, e - RL, True)):
            # trim the read to RL bases (insertions lengthen it)
            a = pysam.AlignedSegment(header)
            a.query_name = name
            a.reference_id = 0
            a.reference_start = start
            a.cigartuples = [({"M": 0, "I": 1, "D": 2}[k_], n) for k_, n in ops]
            a.query_sequence = rs
            a.query_qualities = pysam.qualitystring_to_array("F" * len(rs))
            a.mapping_quality = 60
            a.flag = 0x1 | 0x2 | (0x40 if which == 1 else 0x80) | (0x10 if rev else 0x20)
            a.next_reference_id = 0
            a.next_reference_start = (e - RL) if which == 1 else s
            a.template_length = (e - s) if which == 1 else -(e - s)
            a.set_tag("NM", nm(rs, ops, start))
            reads.append(a)

with pysam.AlignmentFile(str(out / "reads.bam"), "wb", header=header) as bam:
    for r in sorted(reads, key=lambda r: (r.reference_start, r.query_name)):
        bam.write(r)
pysam.index(str(out / "reads.bam"))

with open(out / "calls.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=virus,length=3000>\n"
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
            '##INFO=<ID=KIND,Number=1,Type=String,Description="What the test site imitates">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos1, ref, alt, af, kind, mode in variants:
        vaf = 0.7 if mode == "primer-artefact" else af
        f.write(f"virus\t{pos1}\t.\t{ref}\t{alt}\t60\tPASS\tAF={vaf};KIND={kind.replace(' ', '_')}\n")
print(f"wrote {out}: {len(reads)} reads, {len(amplicons)} amplicons, {len(variants)} variants")
for v in variants:
    print("  ", v[0], v[1], v[2], v[3], v[4])
