"""
Synthetic data for read patterns (user manual, section 9.5): pairs around one
site, from a reference majority, two real minority lineages, an artefact
seen only at read ends on one strand, a contaminant with scattered
mismatches, PCR chimeras of the two lineages, and background sequencing
errors.

Usage: python simulate_patterns.py OUT_DIR   -> ref.fa, reads.bam(.bai), calls.vcf
Build with --bam-processing skip.
"""
import random
import sys
from pathlib import Path

import pysam

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(5)
L, RL = 2000, 150
seq = "".join(rng.choice("ACGT") for _ in range(L))
(out / "ref.fa").write_text(">virus\n" + "\n".join(seq[i:i + 60] for i in range(0, L, 60)) + "\n")
pysam.faidx(str(out / "ref.fa"))
header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                          "SQ": [{"SN": "virus", "LN": L}]})


def alt_of(pos1, k=1):
    b = seq[pos1 - 1]
    return [x for x in "ACGT" if x != b][k % 3]


SITE = 1000
# Lineage B: site + 930 + 1110.  Lineage C: site + 1060.  Artefact: site + 1005.
B = {930: alt_of(930), SITE: alt_of(SITE), 1110: alt_of(1110)}
C = {SITE: alt_of(SITE), 1060: alt_of(1060)}
ART = {SITE: alt_of(SITE), 1005: alt_of(1005, 2)}
reads = []
n = 0


def make_read(name, start0, hap, first, reverse, mapq=60, extra=(), errors=0.002):
    bases, quals = [], []
    for i in range(RL):
        p1 = start0 + i + 1
        b = hap.get(p1, seq[p1 - 1])
        q = 37
        if p1 in extra:
            b = extra[p1]
        elif rng.random() < errors:              # sequencing error, often low quality
            b = rng.choice([x for x in "ACGT" if x != b])
            q = rng.choice((8, 12, 30, 37))
        bases.append(b)
        quals.append(q)
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = start0
    a.cigarstring = f"{RL}M"
    a.query_sequence = "".join(bases)
    a.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))
    a.mapping_quality = mapq
    a.set_tag("NM", sum(1 for i, b in enumerate(bases) if b != seq[start0 + i]))
    return a, first, reverse


def pair(hap, start0=None, insert=None, kind="", mapq=60, scatter=0, chimera=None, end_at=None):
    global n
    n += 1
    insert = insert or rng.randint(320, 380)
    if end_at is not None:
        start0 = end_at - insert
    elif start0 is None:
        start0 = rng.randint(SITE - 1 - RL + 5 if rng.random() < .5 else SITE - insert + 10, SITE - 6)
        start0 = max(0, start0)
    end0 = start0 + insert
    name = f"{kind}{n}"
    extra = {}
    if scatter:
        for p in rng.sample(range(start0 + 1, end0 + 1), scatter):
            extra[p] = alt_of(p, 2)
    h1 = h2 = hap
    if chimera:                                   # left part from one lineage, right from another
        h1, h2 = chimera
        h1 = {p: b for p, b in h1.items() if p < 1030} | {p: b for p, b in h2.items() if p >= 1030}
        h2 = h1
    r1 = make_read(name, start0, h1, True, False, mapq, extra)
    r2 = make_read(name, end0 - RL, h2, False, True, mapq, extra)
    for (a, first, rev), (m, _, _) in ((r1, r2), (r2, r1)):
        a.flag = 0x1 | 0x2 | (0x40 if first else 0x80) | (0x10 if rev else 0x20)
        a.next_reference_id = 0
        a.next_reference_start = m.reference_start
        a.template_length = insert if first else -insert
        reads.append(a)


for _ in range(300):
    pair({}, kind="ref")
for _ in range(130):
    pair(B, kind="B")
for _ in range(60):
    pair(C, kind="C")
for _ in range(30):                               # read 2 ends 7 bp after the site, on the reverse strand
    pair(ART, end_at=SITE + 7, kind="art")
for _ in range(14):
    pair(C, kind="cont", mapq=25, scatter=12)
for _ in range(12):
    pair(None, kind="chim", chimera=(B, C))

with pysam.AlignmentFile(str(out / "reads.bam"), "wb", header=header) as bam:
    for r in sorted(reads, key=lambda r: (r.reference_start, r.query_name)):
        bam.write(r)
pysam.index(str(out / "reads.bam"))

with open(out / "calls.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=virus,length=2000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for p, b in sorted({930: B[930], SITE: B[SITE], 1060: C[1060], 1110: B[1110]}.items()):
        f.write(f"virus\t{p}\t.\t{seq[p - 1]}\t{b}\t60\tPASS\t.\n")
print(f"wrote {out}: {len(reads)} reads; site {SITE}; B={B} C={C} artefact={ART}")
