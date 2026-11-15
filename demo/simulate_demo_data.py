#!/usr/bin/env python3
"""Simulate the small synthetic paired-end dataset in demo/.

Usage:  python simulate_demo_data.py OUT_DIR

Writes demo_ref.fa, demo_reads.bam (coordinate sorted, indexed) and
demo_calls.vcf to OUT_DIR. The random seed is fixed, so the output is the same
every time. Everything is synthetic; the alignments are written directly
(truth alignments), not produced by an aligner.

The committed demo_reads.bam additionally has MD/NM tags and the VCF is
bgzipped and indexed; see demo/README.md for the exact commands.
"""
import random
import sys
from pathlib import Path

import pysam

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
rng = random.Random(20261115)

READ_LEN = 150
CONTIGS = [("demo_seg1", 6000), ("demo_seg2", 3000), ("demo_seg3", 1500)]


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


ref = {name: rand_seq(length) for name, length in CONTIGS}
with open(OUT / "demo_ref.fa", "w") as fh:
    for name, _ in CONTIGS:
        fh.write(f">{name}\n")
        s = ref[name]
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")


def other_base(b, k=0):
    return [x for x in "ACGT" if x != b][k]


S1 = ref["demo_seg1"]
# 0-based positions on demo_seg1
# kind, pos0, payload, allele fraction, strand restriction
VARS = [
    dict(kind="snp", pos=1199, alt=other_base(S1[1199]), af=0.97, strand=None, called=True),
    dict(kind="snp", pos=2449, alt=other_base(S1[2449]), af=0.30, strand=None, called=True),
    dict(kind="del", pos=3300, length=7, af=0.40, strand=None, called=True),       # deletes 3300..3306 (0-based)
    dict(kind="ins", pos=4100, seq="GAT", af=0.25, strand=None, called=True),      # inserted before 0-based 4100
    dict(kind="snp", pos=4999, alt=other_base(S1[4999]), af=0.32, strand="-", called=True),  # reverse-strand-only artefact
    # uncalled, low-frequency mixed site
    dict(kind="snp", pos=2799, alt=other_base(S1[2799], 0), af=0.08, strand=None, called=False),
    dict(kind="snp", pos=2799, alt=other_base(S1[2799], 1), af=0.04, strand=None, called=False),
]


def choose_variants():
    """Pick the variants carried by one fragment (one haplotype)."""
    chosen = []
    taken = set()
    for v in VARS:
        if v["pos"] in taken:
            continue
        if rng.random() < v["af"] / (1.0 if v["strand"] is None else 1.0):
            chosen.append(v)
            taken.add(v["pos"])
    return chosen


def walk(contig, start, nquery, variants, reverse):
    """Walk the reference from `start`, emitting `nquery` read bases.
    Returns (seq, cigartuples, ref_end)."""
    s = ref[contig]
    by_pos = {}
    if contig == "demo_seg1":
        for v in variants:
            if v["strand"] == "-" and not reverse:
                continue
            by_pos[v["pos"]] = v
    seq = []
    cig = []

    def push(op, n=1):
        if cig and cig[-1][0] == op:
            cig[-1][1] += n
        else:
            cig.append([op, n])

    p = start
    while len(seq) < nquery and p < len(s):
        v = by_pos.get(p)
        if v and v["kind"] == "del" and seq:
            push(2, v["length"])
            p += v["length"]
            continue
        if v and v["kind"] == "ins" and seq and len(seq) + len(v["seq"]) < nquery:
            seq.extend(v["seq"])
            push(1, len(v["seq"]))
            by_pos.pop(p)
            continue
        b = s[p]
        if v and v["kind"] == "snp":
            b = v["alt"]
        elif rng.random() < 0.002:           # sequencing error
            b = other_base(b, rng.randrange(3))
        seq.append(b)
        push(0)
        p += 1
    # never end on a deletion/insertion
    while cig and cig[-1][0] in (1, 2):
        op, n = cig.pop()
        if op == 1:
            del seq[-n:]
        else:
            p -= n
    return "".join(seq), [tuple(c) for c in cig], p


def quals(n, reverse):
    q = []
    for i in range(n):
        d = (n - 1 - i) if reverse else i   # distance from the 5' end of the read
        base = 37 if d < 110 else rng.choice([37, 37, 30, 24])
        if rng.random() < 0.03:
            base = rng.choice([12, 18, 22])
        q.append(base)
    return q


header = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "unsorted"},
    "SQ": [{"SN": n, "LN": l} for n, l in CONTIGS],
    "PG": [{"ID": "simulate", "PN": "simulate.py", "CL": "synthetic demo data for the Allelescope manual"}],
})
tid = {n: i for i, (n, _) in enumerate(CONTIGS)}
records = []
counter = [0]


def new_name():
    counter[0] += 1
    return f"DEMO:1:1101:{rng.randrange(1000, 30000)}:{counter[0]:05d}"


def make(name, contig, pos, seq, cig, reverse, read1, proper, qual=None, extra_flag=0, tags=()):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = tid[contig]
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigartuples = cig
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + x) for x in (qual or quals(len(seq), reverse))))
    flag = 0x1 | (0x40 if read1 else 0x80) | (0x10 if reverse else 0) | (0x2 if proper else 0) | extra_flag
    a.flag = flag
    for k, v in tags:
        a.set_tag(k, v)
    return a


def set_mates(a, b, proper_tlen=True):
    for x, y in ((a, b), (b, a)):
        x.next_reference_id = y.reference_id
        x.next_reference_start = y.reference_start
        if y.is_reverse:
            x.flag |= 0x20
    if a.reference_id == b.reference_id:
        left, right = (a, b) if a.reference_start <= b.reference_start else (b, a)
        tlen = right.reference_end - left.reference_start
        left.template_length = tlen
        right.template_length = -tlen


def softclip(seq, cig, pos, side, k):
    """Turn k bases at one end into a soft clip of random sequence (walk again is overkill)."""
    junk = rand_seq(k)
    if side == "L":
        # consume k query bases from the left of the alignment
        newcig = list(map(list, cig))
        need = k
        refshift = 0
        while need > 0 and newcig:
            op, n = newcig[0]
            if op == 0:
                t = min(n, need)
                newcig[0][1] -= t
                need -= t
                refshift += t
                if newcig[0][1] == 0:
                    newcig.pop(0)
            elif op == 1:
                t = min(n, need)
                newcig[0][1] -= t
                need -= t
                if newcig[0][1] == 0:
                    newcig.pop(0)
            else:  # deletion
                refshift += n
                newcig.pop(0)
        while newcig and newcig[0][0] in (1, 2):
            op, n = newcig.pop(0)
            if op == 2:
                refshift += n
            else:
                k += n
                junk += rand_seq(n)
        return junk + seq[k:], [(4, k)] + [tuple(c) for c in newcig], pos + refshift
    else:
        newcig = list(map(list, cig))
        need = k
        while need > 0 and newcig:
            op, n = newcig[-1]
            if op in (0, 1):
                t = min(n, need)
                newcig[-1][1] -= t
                need -= t
                if newcig[-1][1] == 0:
                    newcig.pop()
            else:
                newcig.pop()
        while newcig and newcig[-1][0] in (1, 2):
            op, n = newcig.pop()
            if op == 1:
                k += n
                junk += rand_seq(n)
        return seq[:len(seq) - k] + junk[:k], [tuple(c) for c in newcig] + [(4, k)], pos


def normal_pair(contig, frag_start, frag_len, variants, orientation="FR"):
    name = new_name()
    r1_rev = False
    # read 1 on forward strand half the time
    swap = rng.random() < 0.5
    s1, c1, e1 = walk(contig, frag_start, READ_LEN, variants, reverse=False)
    start2 = max(0, frag_start + frag_len - READ_LEN)
    rev2 = orientation == "FR"
    s2, c2, e2 = walk(contig, start2, READ_LEN, variants, reverse=rev2)
    p1, p2 = frag_start, start2
    # occasional soft clips
    if rng.random() < 0.05:
        s1, c1, p1 = softclip(s1, c1, p1, "L", rng.randrange(6, 45))
    if rng.random() < 0.05:
        s2, c2, p2 = softclip(s2, c2, p2, "R", rng.randrange(6, 45))
    proper = orientation == "FR"
    left = make(name, contig, p1, s1, c1, False, not swap, proper)
    right = make(name, contig, p2, s2, c2, rev2, swap, proper)
    set_mates(left, right)
    return [left, right]


def add_fragments(contig, n, mean=300, sd=60):
    L = len(ref[contig])
    out = []
    for _ in range(n):
        flen = max(READ_LEN + 20, int(rng.gauss(mean, sd)))
        start = rng.randrange(-flen + READ_LEN, L - READ_LEN)
        start = max(0, start)
        flen = min(flen, L - start)
        if flen < READ_LEN:
            continue
        variants = choose_variants() if contig == "demo_seg1" else []
        pair = normal_pair(contig, start, flen, variants)
        out.extend(pair)
        # PCR / optical duplicates
        if rng.random() < 0.08:
            dup_name = new_name()
            for r in pair:
                d = pysam.AlignedSegment.fromstring(r.tostring(), header)
                d.query_name = dup_name
                out.append(d)
    return out


records += add_fragments("demo_seg1", 1700)
records += add_fragments("demo_seg2", 420)
records += add_fragments("demo_seg3", 110)

# Discordant pairs: both mates on the forward strand (F1F2), around demo_seg1:5400
for _ in range(14):
    start = rng.randrange(5150, 5380)
    records += normal_pair("demo_seg1", start, rng.randrange(300, 420), choose_variants(), orientation="FF")

# Split (chimeric) reads: demo_seg1 ...3700 | demo_seg2 800...
JA, JB = 3700, 800     # 0-based: seg1 bases < 3700 joined to seg2 bases >= 800
for _ in range(16):
    name = new_name()
    x = rng.randrange(45, 106)                 # bases of read 1 on demo_seg1
    y = READ_LEN - x
    partA = S1[JA - x:JA]
    partB = ref["demo_seg2"][JB:JB + y]
    full = partA + partB
    q = quals(READ_LEN, False)
    if x >= y:
        prim = make(name, "demo_seg1", JA - x, full, [(0, x), (4, y)], False, True, False, qual=q,
                    tags=[("SA", f"demo_seg2,{JB + 1},+,{x}S{y}M,60,0;")])
        supp = make(name, "demo_seg2", JB, partB, [(5, x), (0, y)], False, True, False, qual=q[x:],
                    extra_flag=0x800, tags=[("SA", f"demo_seg1,{JA - x + 1},+,{x}M{y}S,60,0;")])
    else:
        prim = make(name, "demo_seg2", JB, full, [(4, x), (0, y)], False, True, False, qual=q,
                    tags=[("SA", f"demo_seg1,{JA - x + 1},+,{x}M{y}S,60,0;")])
        supp = make(name, "demo_seg1", JA - x, partA, [(0, x), (5, y)], False, True, False, qual=q[:x],
                    extra_flag=0x800, tags=[("SA", f"demo_seg2,{JB + 1},+,{x}S{y}M,60,0;")])
    # mate: reverse read downstream on demo_seg2
    mstart = JB + rng.randrange(150, 260)
    ms, mc, _ = walk("demo_seg2", mstart, READ_LEN, [], reverse=True)
    mate = make(name, "demo_seg2", mstart, ms, mc, True, False, False)
    set_mates(prim, mate)
    supp.next_reference_id = mate.reference_id
    supp.next_reference_start = mate.reference_start
    supp.flag |= 0x20
    records += [prim, supp, mate]

tmp = OUT / "unsorted.bam"
with pysam.AlignmentFile(str(tmp), "wb", header=header) as bam:
    for r in records:
        bam.write(r)
pysam.sort("-o", str(OUT / "demo_reads.bam"), str(tmp))
tmp.unlink()
pysam.index(str(OUT / "demo_reads.bam"))

# ---- VCF with the called variants; DP taken from the simulated pile-up ----
bam = pysam.AlignmentFile(str(OUT / "demo_reads.bam"))
lines = []
for v in VARS:
    if not v["called"]:
        continue
    if v["kind"] == "snp":
        pos1 = v["pos"] + 1
        r, a = S1[v["pos"]], v["alt"]
    elif v["kind"] == "del":
        pos1 = v["pos"]                       # anchor base (1-based) = 0-based pos-1 +1
        r = S1[v["pos"] - 1:v["pos"] + v["length"]]
        a = S1[v["pos"] - 1]
    else:
        pos1 = v["pos"]
        r = S1[v["pos"] - 1]
        a = r + v["seq"]
    dp = bam.count("demo_seg1", pos1 - 1, pos1)
    lines.append((pos1, f"demo_seg1\t{pos1}\t.\t{r}\t{a}\t60\tPASS\tDP={dp}"))
bam.close()
with open(OUT / "demo_calls.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    fh.write('##source=simulate.py (synthetic demo data)\n')
    for n, l in CONTIGS:
        fh.write(f"##contig=<ID={n},length={l}>\n")
    fh.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for _, l in sorted(lines):
        fh.write(l + "\n")
print("reads:", len(records))
for _, l in sorted(lines):
    print(l)
