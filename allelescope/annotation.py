"""
Gene annotation and amplicon primers for an instance.

Coding effects
    Read CDS features from a GFF3 or GTF file and work out what a variant does
    to the protein: synonymous, missense, stop gained or lost, start lost,
    frameshift, in-frame insertion or deletion. Multi-part CDSs (spliced
    genes, or ribosomal slippage as in coronavirus ORF1ab, written as two
    parts with the same ID) are joined in transcript order, on either strand.

Primers
    Read an amplicon primer scheme (BED: chrom, start, end, name, [pool],
    [strand]) and say which primers cover a site, and how many reads of an
    allele have primer sequence there: reads whose aligned start lies in a
    left (+) primer covering the site, or whose aligned end lies in a right
    (-) primer covering it. Primer-trimmed (soft-clipped) reads do not count,
    because their aligned ends lie past the primer.

All positions here are 1-based and inclusive, as in VCF.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from urllib.parse import unquote

import pysam

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
# NCBI genetic codes that differ from the standard one, by codon.
CODE_CHANGES = {
    2: {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},
    3: {"ATA": "M", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "TGA": "W"},
    4: {"TGA": "W"},
    5: {"AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W"},
    11: {},
}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# A read end this many bases outside a primer still counts as starting in it.
PRIMER_SLACK = 3


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(codon: str, table: int = 1) -> str:
    codon = codon.upper()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "X"
    return CODE_CHANGES.get(table, {}).get(codon, CODON_TABLE[codon])


def _attributes(text: str) -> dict:
    """GFF3 (key=value;...) or GTF (key "value"; ...) attributes."""
    out = {}
    text = text.strip()
    if not text or text == ".":
        return out
    if "=" in text.split(";")[0]:
        for part in text.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                out[k.strip()] = unquote(v.strip())
    else:
        for k, v in re.findall(r'(\S+)\s+"([^"]*)"', text):
            out.setdefault(k, v)
    return out


@dataclass
class CDS:
    key: str
    seqid: str
    strand: str
    parts: list = field(default_factory=list)     # [(start, end, phase)] 1-based inclusive
    name: str = ""
    product: str = ""
    table: int = 1
    # filled by prepare()
    seq: str = ""
    frame: int = 0
    offsets: list = field(default_factory=list)  # [(start, end, cds offset of the 5' end)]

    @property
    def start(self):
        return min(s for s, _, _ in self.parts)

    @property
    def end(self):
        return max(e for _, e, _ in self.parts)

    def prepare(self, fasta: pysam.FastaFile):
        rev = self.strand == "-"
        ordered = sorted(self.parts, key=lambda p: p[0], reverse=rev)
        pieces, off = [], 0
        self.offsets = []
        for s, e, _ in ordered:
            piece = fasta.fetch(self.seqid, s - 1, e).upper()
            pieces.append(revcomp(piece) if rev else piece)
            self.offsets.append((s, e, off))
            off += e - s + 1
        self.seq = "".join(pieces)
        phase = ordered[0][2]
        self.frame = phase if isinstance(phase, int) and 0 <= phase <= 2 else 0

    def cds_offsets(self, pos: int) -> list:
        """CDS offsets (0-based) of genomic position pos; two in a slippage overlap."""
        out = []
        for s, e, off in self.offsets:
            if s <= pos <= e:
                out.append(off + (e - pos if self.strand == "-" else pos - s))
        return out


class Annotation:
    """CDS features and primers of one instance, with a reference to read codons from."""

    def __init__(self, fasta_path, gff_path=None, primer_path=None):
        self.fasta = pysam.FastaFile(str(fasta_path))
        self.contigs = set(self.fasta.references)
        self.cds: dict[str, list[CDS]] = {}
        self.primers: dict[str, list[tuple]] = {}
        self.gff_seqids: set = set()
        self.primer_seqids: set = set()
        self.gff_path = Path(gff_path) if gff_path else None
        self.primer_path = Path(primer_path) if primer_path else None
        if gff_path:
            self._read_gff(gff_path)
        if primer_path:
            self._read_primers(primer_path)

    # ---- reading ----------------------------------------------------------
    def _read_gff(self, path):
        by_key: dict[str, CDS] = {}
        names: dict[str, dict] = {}           # feature ID -> attributes, for name lookup
        cds_lines = []
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("##FASTA"):
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                seqid, _, ftype, start, end, _, strand, phase, attr = cols[:9]
                a = _attributes(attr)
                if "ID" in a:
                    names.setdefault(a["ID"], a)
                if ftype == "CDS":
                    try:
                        s, e = int(start), int(end)
                    except ValueError:
                        continue
                    ph = int(phase) if phase in ("0", "1", "2") else 0
                    cds_lines.append((seqid, s, e, strand, ph, a))
        for seqid, s, e, strand, ph, a in cds_lines:
            self.gff_seqids.add(seqid)
            if seqid not in self.contigs:
                continue
            key = a.get("ID") or a.get("Parent") or a.get("transcript_id") or a.get("protein_id") \
                or f"{seqid}:{s}-{e}"
            c = by_key.get((seqid, key))
            if c is None:
                c = by_key[(seqid, key)] = CDS(key=key, seqid=seqid, strand=strand)
                c.name = self._feature_name(a, names)
                c.product = a.get("product", "")
                try:
                    c.table = int(a.get("transl_table", 1))
                except ValueError:
                    c.table = 1
            c.parts.append((s, e, ph))
        # Two CDSs with the same gene name (ORF1a and ORF1ab): tell them apart.
        seen: dict[str, list[CDS]] = {}
        for c in by_key.values():
            seen.setdefault(c.name, []).append(c)
        for same in seen.values():
            if len(same) > 1:
                for c in same:
                    label = c.product or c.key
                    c.name = f"{c.name} [{label}]" if label and label != c.name else f"{c.name} [{c.key}]"
        for c in by_key.values():
            c.prepare(self.fasta)
            self.cds.setdefault(c.seqid, []).append(c)
        for lst in self.cds.values():
            lst.sort(key=lambda c: c.start)

    @staticmethod
    def _feature_name(a, names):
        for k in ("gene", "Name", "gene_name", "locus_tag"):
            if a.get(k):
                return a[k]
        parent = a.get("Parent", "").split(",")[0]
        for _ in range(3):                      # CDS -> mRNA -> gene
            p = names.get(parent)
            if not p:
                break
            for k in ("gene", "Name", "gene_name", "locus_tag"):
                if p.get(k):
                    return p[k]
            parent = p.get("Parent", "").split(",")[0]
        return a.get("product") or a.get("protein_id") or a.get("ID") or "CDS"

    def _read_primers(self, path):
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3:
                    cols = line.split()
                if len(cols) < 3:
                    continue
                try:
                    s0, e = int(cols[1]), int(cols[2])
                except ValueError:
                    continue
                name = cols[3] if len(cols) > 3 else f"{cols[0]}:{s0 + 1}-{e}"
                strand = cols[5] if len(cols) > 5 and cols[5] in "+-" else ""
                if not strand:
                    u = name.upper()
                    strand = "+" if ("LEFT" in u or u.endswith("_F")) else "-" if ("RIGHT" in u or u.endswith("_R")) else ""
                self.primer_seqids.add(cols[0])
                self.primers.setdefault(cols[0], []).append((s0 + 1, e, name, strand))

    # ---- summary ----------------------------------------------------------
    def info(self) -> dict:
        n_cds = sum(len(v) for v in self.cds.values())
        n_primers = sum(len(v) for k, v in self.primers.items() if k in self.contigs)
        out = {"cds": n_cds, "primers": n_primers, "warnings": []}
        if self.gff_path:
            out["annotation_file"] = self.gff_path.name
            missing = sorted(self.gff_seqids - self.contigs)
            if missing and not n_cds:
                out["warnings"].append(
                    f"{self.gff_path.name}: no sequence names match the reference "
                    f"(annotation has {', '.join(missing[:3])}; reference has {', '.join(sorted(self.contigs)[:3])})")
            elif missing:
                out["warnings"].append(f"{self.gff_path.name}: {len(missing)} sequence name(s) not in the reference, ignored")
        if self.primer_path:
            out["primer_file"] = self.primer_path.name
            missing = sorted(self.primer_seqids - self.contigs)
            if missing and not n_primers:
                out["warnings"].append(
                    f"{self.primer_path.name}: no sequence names match the reference (primers on {', '.join(missing[:3])})")
        return out

    # ---- primers ----------------------------------------------------------
    def primers_at(self, chrom: str, pos: int) -> list:
        return [p for p in self.primers.get(chrom, ()) if p[0] <= pos <= p[1]]

    def primer_reads(self, chrom: str, pos: int, starts: dict, ends: dict) -> int:
        """
        Reads with primer sequence at pos, from the counts of aligned starts and
        ends (1-based, inclusive) of the reads of one allele.
        """
        cover = self.primers_at(chrom, pos)
        if not cover:
            return 0
        n = 0
        counted_starts, counted_ends = set(), set()
        for s, e, _, strand in cover:
            if strand in ("+", ""):
                for st, c in starts.items():
                    st = int(st)
                    if s - PRIMER_SLACK <= st <= pos and st not in counted_starts:
                        counted_starts.add(st)
                        n += c
            if strand in ("-", ""):
                for en, c in ends.items():
                    en = int(en)
                    if pos <= en <= e + PRIMER_SLACK and en not in counted_ends:
                        counted_ends.add(en)
                        n += c
        return n

    # ---- coding effects ---------------------------------------------------
    def has_cds(self, chrom: str) -> bool:
        return bool(self.cds.get(chrom))

    def effects(self, chrom: str, pos: int, ref: str, alt: str) -> list:
        """
        Effects of replacing ref at pos by alt (VCF style; ref and alt may share
        leading or trailing bases, and either may be empty after trimming).
        Returns [{"gene", "change", "type", "codon"}]; [] outside every CDS.
        """
        ref, alt = (ref or "").upper(), (alt or "").upper()
        if alt in ("", ".", "*", "<DEL>") and not ref:
            return []
        if alt.startswith("<") or alt in (".", "*"):
            return []
        # Trim shared bases: suffix first, then prefix.
        while ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        while ref and alt and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
        if not ref and not alt:
            return []
        out = []
        for c in self.cds.get(chrom, ()):
            if len(ref) == len(alt):
                e = self._substitution(c, pos, ref, alt)
            else:
                e = self._indel(c, pos, ref, alt)
            if e:
                out.append(e)
        return out

    def _codon_no(self, c: CDS, off: int) -> int:
        return (off - c.frame) // 3 + 1

    def _substitution(self, c: CDS, pos, ref, alt):
        changes = {}                            # cds offset -> new base (coding strand)
        for i, b in enumerate(alt):
            for off in c.cds_offsets(pos + i):
                changes[off] = revcomp(b) if c.strand == "-" else b
        if not changes:
            return None
        seq = list(c.seq)
        for off, b in changes.items():
            if off < len(seq):
                seq[off] = b
        mutated = "".join(seq)
        codons = sorted({(off - c.frame) // 3 for off in changes if off >= c.frame})
        if not codons:
            return {"gene": c.name, "change": "", "type": "partial codon", "codon": None}
        parts, types = [], set()
        for k in codons:
            a, b = c.frame + 3 * k, c.frame + 3 * k + 3
            if b > len(c.seq):
                continue
            ra, aa = translate(c.seq[a:b], c.table), translate(mutated[a:b], c.table)
            n = k + 1
            parts.append(f"{ra}{n}{aa}")
            if ra == aa:
                types.add("synonymous")
            elif aa == "*":
                types.add("stop gained")
            elif ra == "*":
                types.add("stop lost")
            elif n == 1 and ra == "M":
                types.add("start lost")
            else:
                types.add("missense")
        if not parts:
            return None
        for t in ("stop gained", "start lost", "stop lost", "missense", "synonymous"):
            if t in types:
                kind = t
                break
        return {"gene": c.name, "change": ",".join(parts), "type": kind, "codon": codons[0] + 1}

    def _indel(self, c: CDS, pos, ref, alt):
        if ref:                                  # deletion (or complex): bases pos .. pos+len-1
            span = range(pos, pos + len(ref))
            offs = sorted(o for p in span for o in c.cds_offsets(p))
            if not offs:
                return None
            inside = all(c.cds_offsets(p) for p in span)
        else:                                    # pure insertion before pos
            offs = sorted(set(c.cds_offsets(pos)) | set(c.cds_offsets(pos - 1)))
            if not c.cds_offsets(pos) or not c.cds_offsets(pos - 1):
                return None                      # at the edge of the CDS
            inside = True
        first = offs[0]
        n = self._codon_no(c, max(first, c.frame))
        k = (max(first, c.frame) - c.frame) // 3
        ref_aa = translate(c.seq[c.frame + 3 * k: c.frame + 3 * k + 3], c.table)
        diff = len(alt) - len(ref)
        if not inside:
            return {"gene": c.name, "change": f"{ref_aa}{n}?", "type": "CDS boundary", "codon": n}
        if diff % 3:
            return {"gene": c.name, "change": f"{ref_aa}{n}fs", "type": "frameshift", "codon": n}
        if diff < 0:
            return {"gene": c.name, "change": f"{ref_aa}{n}del{-diff // 3}", "type": "in-frame deletion", "codon": n}
        return {"gene": c.name, "change": f"{ref_aa}{n}ins{diff // 3}", "type": "in-frame insertion", "codon": n}

    def allele_effects(self, chrom: str, pos: int, code: str) -> list:
        """Effects of an allele code from extract_mates (M<base>, D<n>, I<seq>) at site pos."""
        if not code or code[0] not in "MDI":
            return []
        if code[0] == "M":
            ref = self.fasta.fetch(chrom, pos - 1, pos)
            return self.effects(chrom, pos, ref, code[1:])
        if code[0] == "D":
            try:
                n = int(code[1:])
            except ValueError:
                return []
            ref = self.fasta.fetch(chrom, pos - 1, pos - 1 + n)
            return self.effects(chrom, pos, ref, "")
        return self.effects(chrom, pos, "", code[1:])


# ---- finding the files of an instance -----------------------------------
ANNOTATION_NAMES = ("annotation.gff3", "annotation.gff", "annotation.gtf",
                    "annotation.gff3.gz", "annotation.gff.gz")
PRIMER_NAMES = ("primers.bed", "primer.bed", "scheme.primer.bed")


def instance_files(serverpath: Path):
    """(reference FASTA, annotation file or None, primer BED or None) of an instance."""
    serverpath = Path(serverpath)
    fasta = None
    for fai in sorted(serverpath.glob("*.fai")):
        candidate = fai.with_suffix("")
        if candidate.is_file():
            fasta = candidate
            break
    gff = next((serverpath / n for n in ANNOTATION_NAMES if (serverpath / n).is_file()), None)
    bed = next((serverpath / n for n in PRIMER_NAMES if (serverpath / n).is_file()), None)
    return fasta, gff, bed


def open_annotation(path) -> Optional[str]:
    """Plain-text copy of a gzipped annotation (opened once and kept)."""
    p = Path(path)
    if p.suffix != ".gz":
        return str(p)
    import gzip
    plain = p.with_suffix("")
    if not plain.exists() or plain.stat().st_mtime < p.stat().st_mtime:
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as src, \
                open(plain, "w", encoding="utf-8") as dst:
            dst.write(src.read())
    return str(plain)
