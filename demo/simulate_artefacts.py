"""Synthetic test data for the read-evidence warnings (user manual, section 9.3):
one clean variant and one site per kind of artefact.
Usage: python simulate_artefacts.py OUT_DIR   -> ref.fa, reads.bam(.bai), calls.vcf"""
import random
import sys
from pathlib import Path

import pysam

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(7)
L = 6000
seq = "".join(rng.choice("ACGT") for _ in range(L))
(out / "ref.fa").write_text(">virus\n" + "\n".join(seq[i:i + 60] for i in range(0, L, 60)) + "\n")
pysam.faidx(str(out / "ref.fa"))
header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                          "SQ": [{"SN": "virus", "LN": L}],
                                          "RG": [{"ID": "s", "SM": "s"}]})
RL = 100
reads = []
n_frag = 0


def other(b):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


def pair(start, insert, alt_site=None, alt=None, read1_fwd=True, crossing_first=True,
         extra_mm=0, xa=False, mapq=60, dup=False, clip=0):
    """One read pair; the read starting at `start` crosses the site."""
    global n_frag
    n_frag += 1
    name = f"frag{n_frag}"
    mstart = start + insert - RL
    recs = []
    for crossing in (True, False):
        s = start if crossing else mstart
        bases = list(seq[s:s + RL])
        nm = 0
        if crossing and alt_site is not None:
            bases[alt_site - s] = alt
            nm += 1
            for k in range(extra_mm):           # paralog-like extra mismatches
                j = (alt_site - s + 7 * (k + 1)) % RL
                bases[j] = other(bases[j])
                nm += 1
        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.reference_id = 0
        a.reference_start = s
        cigar = f"{RL}M"
        if crossing and clip:
            # the last `clip` bases do not match: soft-clipped
            cigar = f"{RL - clip}M{clip}S"
            for j in range(RL - clip, RL):
                bases[j] = other(seq[s + j])
        a.cigarstring = cigar
        a.query_sequence = "".join(bases)
        a.query_qualities = pysam.qualitystring_to_array("?" * RL)
        a.mapping_quality = mapq if crossing else 60
        is_r1 = crossing == crossing_first
        # read 1 forward / read 2 reverse is F1R2
        fwd = read1_fwd if is_r1 else not read1_fwd
        left = s == min(start, mstart)
        a.flag = (0x1 | 0x2 | (0x40 if is_r1 else 0x80) | (0 if fwd else 0x10)
                  | (0x20 if fwd else 0) | (0x400 if dup else 0))
        a.next_reference_id = 0
        a.next_reference_start = mstart if crossing else start
        a.template_length = insert if left else -insert
        a.set_tag("NM", nm)
        a.set_tag("RG", "s")
        if crossing and xa:
            a.set_tag("XA", "virus,+4000,100M,2;")
        recs.append(a)
    reads.extend(recs)


def site(pos1, depth, alt_frac, **alt_kw):
    """Reads around a 1-based site; alt reads get alt_kw."""
    site0 = pos1 - 1
    alt = other(seq[site0])
    n_alt = round(depth * alt_frac)
    for i in range(depth):
        is_alt = i < n_alt
        kw = alt_kw if is_alt else {}
        start = kw.get("start_fn", lambda i: site0 - rng.randint(10, 89))(i)
        r1f = kw.get("r1f_fn", lambda i: bool(i % 2))(i)
        cf = kw.get("cf_fn", lambda i: bool((i // 2) % 2))(i)
        pass_kw = {k: v for k, v in kw.items() if not k.endswith("_fn")}
        pair(start, pass_kw.pop("insert", 300), site0 if is_alt else None, alt,
             read1_fwd=r1f, crossing_first=cf, **pass_kw)
    return alt


calls = []
# 1. Clean variant, 40 %.
calls.append((500, site(500, 60, 0.4), 0.40, "clean"))
# 2. Read-end artefact: alt only 2-3 bp from the end of the read.
calls.append((1200, site(1200, 60, 0.3, start_fn=lambda i: 1199 - 97 + (i % 2)), 0.30, "read_end"))
# 3. Paralog: 4 extra mismatches, half with other hits, lower MAPQ.
calls.append((1900, site(1900, 60, 0.3, extra_mm=4, xa=True, mapq=3), 0.30, "paralog"))
# 4. Orientation (damage-like): every alt pair F1R2, strands still balanced.
calls.append((2600, site(2600, 60, 0.3, r1f_fn=lambda i: True, cf_fn=lambda i: bool(i % 2)), 0.30, "orientation"))
# 5. VCF says 50 %, reads say 10 %.
calls.append((3300, site(3300, 60, 0.1), 0.50, "af_mismatch"))
# 6. Duplicates: alt from two molecules copied many times.
calls.append((4000, site(4000, 60, 0.3, dup=True, start_fn=lambda i: 3999 - 40 - (i % 2)), 0.30, "duplicates"))
# 7. Clipped near the site, long inserts.
calls.append((4500, site(4500, 60, 0.3, clip=8, insert=650, start_fn=lambda i: 4499 - (100 - 8 - 6) + (i % 3)), 0.30, "clipped_insert"))

with pysam.AlignmentFile(str(out / "reads.bam"), "wb", header=header) as bam:
    for r in sorted(reads, key=lambda r: r.reference_start):
        bam.write(r)
pysam.index(str(out / "reads.bam"))

with open(out / "calls.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=virus,length=6000>\n"
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
            '##INFO=<ID=KIND,Number=1,Type=String,Description="What the test site imitates">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos, alt, af, kind in calls:
        f.write(f"virus\t{pos}\t.\t{seq[pos - 1]}\t{alt}\t60\tPASS\tAF={af};KIND={kind}\n")
print("wrote", out, len(reads), "records")
