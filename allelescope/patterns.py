"""
Read patterns at a site: which combinations of differences from the
reference the read pairs crossing a site carry.

Rather than calling haplotypes, this lists the combinations that occur, how
often, and with what read evidence, so that artefacts (combinations tied to
read ends, one strand, one amplicon, extra mismatches, or mixing two common
patterns) can be told from real linked variants.

Steps, on the per-site read set (every pair crossing the site, with mates):

1. Each read pair's differences from the reference: mismatches at base
   quality >= min_bq, deletions (at their first base) and insertions (at the
   base after them), using the allele codes of extract_mates (M<base>,
   D<len>, I<seq>). Mates that disagree at a position leave it unknown.
2. Informative positions: a difference seen in >= min_reads pairs and
   >= min_fraction of the pairs covering that position, plus the variants of
   the variant list ("listed"), whatever their frequency. At most
   MAX_POSITIONS, keeping the commonest.
3. Each pair's state at each informative position: a code, R (reference),
   x (another difference), * (inside a deletion), or unknown (not covered,
   low quality, mates disagree).
4. Pairs with the same states form a pattern; noisy pairs (NOISY_EDITS)
   form patterns apart from clean ones. A pair that does not reach every
   position joins the one pattern it agrees with; if it agrees with several
   it is "uninformative" (it cannot tell them apart).

All positions are 1-based.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pysam

MAX_POSITIONS = 24
# A pair whose reads carry this many other (non-informative, good-quality)
# differences per read on average is "noisy": scattered mismatches typical of
# a contaminant, a related sequence or a damaged molecule. Noisy pairs form
# patterns of their own.
NOISY_EDITS = 4
NEAR_END_BP = 10
CLIP_MIN_BP = 5
UNKNOWN = "."


@dataclass
class Settings:
    min_reads: int = 3          # informative: seen in at least this many pairs ...
    min_fraction: float = 0.02  # ... and this share of the pairs covering the position
    min_bq: int = 20            # mismatches below this base quality are unknown
    common_fraction: float = 0.05   # a pattern is common from this share of pairs ...
    common_reads: int = 5           # ... and this many pairs
    listed: dict = field(default_factory=dict)   # {pos: set(codes)} from the variant list

    @classmethod
    def from_params(cls, params) -> "Settings":
        def num(key, default, conv):
            try:
                v = params.get(key)
                return conv(v) if v not in (None, "") else default
            except (TypeError, ValueError):
                return default
        listed = {}
        for item in (params.get("PatSites") or "").split(","):
            if ":" in item:
                pos, codes = item.split(":", 1)
                try:
                    listed.setdefault(int(pos), set()).update(c for c in codes.split("|") if c)
                except ValueError:
                    pass
        return cls(
            min_reads=max(1, num("PatN", 3, int)),
            min_fraction=max(0.0, num("PatPct", 2.0, float) / 100),
            min_bq=max(0, num("PatBQ", 20, int)),
            common_fraction=max(0.0, num("PatMin", 5.0, float) / 100),
            common_reads=max(1, num("PatMinN", 5, int)),
            listed=listed,
        )

    def key(self) -> str:
        listed = ",".join(f"{p}:{'|'.join(sorted(c))}" for p, c in sorted(self.listed.items()))
        text = f"{self.min_reads}/{self.min_fraction}/{self.min_bq}/{self.common_fraction}/{self.common_reads}/{listed}"
        return hashlib.sha1(text.encode()).hexdigest()[:12]


# ---- per read ----------------------------------------------------------------
import re as _re
_MD_TOKEN = _re.compile(r"(\d+)|(\^[A-Za-z]+)|([A-Za-z])")


def _md_mismatches(read):
    """
    0-based reference positions of the mismatches of a read, from its MD tag
    (fast: no comparison of every base). None if the read has no MD tag.
    """
    if not read.has_tag("MD"):
        return None
    md = read.get_tag("MD")
    out, k = [], 0                         # k: offset along M and D bases
    for num, dele, mm in _MD_TOKEN.findall(md):
        if num:
            k += int(num)
        elif dele:
            k += len(dele) - 1
        else:
            out.append(k)
            k += 1
    # Offsets count reference bases from the alignment start, skipping
    # nothing but insertions (which MD leaves out): map them to positions.
    if not out:
        return []
    positions, r, off = [], read.reference_start, 0
    targets = iter(out)
    t = next(targets, None)
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8, 2):
            while t is not None and t < off + length:
                if op != 2:
                    positions.append(r + (t - off))
                t = next(targets, None)
            off += length
            r += length
        elif op == 3:
            r += length
        if t is None:
            break
    return positions


def read_differences(read, ref_seq: str, ref_start: int, min_bq: int):
    """
    Differences of one read from the reference.
    ref_seq holds the reference from 1-based position ref_start on.
    Returns (diffs {pos: code}, low {pos}, inside_deletion {pos}, (first, last) covered).
    """
    diffs, low, inside = {}, set(), set()
    seq = read.query_sequence or ""
    quals = read.query_qualities
    md = _md_mismatches(read)
    md_set = set(md) if md is not None else None
    q, r = 0, read.reference_start          # r 0-based
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8) and md_set is not None:
            # Only the mismatches from the MD tag (and Ns) need looking at.
            for r0 in [x for x in md_set if r <= x < r + length]:
                pos1, qi = r0 + 1, q + (r0 - r)
                b = seq[qi].upper()
                if b == "N" or (quals is not None and quals[qi] < min_bq):
                    low.add(pos1)
                else:
                    diffs[pos1] = "M" + b
            q += length
            r += length
        elif op in (0, 7, 8):
            for i in range(length):
                pos1 = r + i + 1
                k = pos1 - ref_start
                if 0 <= k < len(ref_seq):
                    b = seq[q + i].upper()
                    if b != ref_seq[k] and b != "N":
                        if quals is not None and quals[q + i] < min_bq:
                            low.add(pos1)
                        else:
                            diffs[pos1] = "M" + b
                    elif b == "N":
                        low.add(pos1)
            q += length
            r += length
        elif op == 1:
            diffs[r + 1] = "I" + seq[q:q + length].upper()
            q += length
        elif op == 4:
            q += length
        elif op in (2, 3):
            if op == 2:
                diffs[r + 1] = f"D{length}"
                inside.update(range(r + 2, r + length + 1))
            r += length
    return diffs, low, inside, (read.reference_start + 1, read.reference_end)


def _soft_clip(read):
    cig = read.cigartuples or ()
    n = 0
    for ops in (cig, tuple(reversed(cig))):
        for op, length in ops:
            if op == 4:
                n = max(n, length)
            elif op != 5:
                break
    return n


@dataclass
class Fragment:
    name: str
    reads: list = field(default_factory=list)   # [(read, diffs, low, inside, span)]

    def state(self, pos: int) -> str:
        """Code, 'R', '*', or UNKNOWN at pos, mates merged."""
        seen = set()
        for _, diffs, low, inside, (s, e) in self.reads:
            if pos in diffs:
                seen.add(diffs[pos])
            elif pos in inside:
                seen.add("*")
            elif pos in low:
                continue
            elif s <= pos <= e:
                seen.add("R")
        if len(seen) == 1:
            return seen.pop()
        return UNKNOWN                       # not covered, low quality, or mates disagree


# ---- the analysis --------------------------------------------------------------
def analyse(bam_path, fasta_path, chrom: str, site: int, settings: Settings, primers=None) -> dict:
    """
    Patterns of the read pairs in bam_path (a per-site read set) around site.
    primers: optional annotation.Annotation, to name the amplicon of each pair.
    Returns a JSON-ready dict; "assignment" maps read names to pattern ids
    ("u" for uninformative).
    """
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        # Mates on another sequence (or far away on this one) are kept for
        # their own sequence's positions only: compare reads with this one.
        records = [r for r in bam if not (r.is_unmapped or r.is_secondary or r.is_supplementary)
                   and r.reference_name == chrom]
    if not records:
        return {"positions": [], "patterns": [], "total": 0, "assignment": {}}

    lo = min(r.reference_start for r in records) + 1
    hi = max(r.reference_end for r in records)
    with pysam.FastaFile(str(fasta_path)) as fa:
        ref_seq = fa.fetch(chrom, lo - 1, hi).upper()

    frags: dict[str, Fragment] = {}
    for r in records:
        f = frags.setdefault(r.query_name, Fragment(r.query_name))
        f.reads.append((r, *read_differences(r, ref_seq, lo, settings.min_bq)))
    frag_list = list(frags.values())

    # 2. Informative positions.
    counts: Counter = Counter()
    for f in frag_list:
        codes = defaultdict(set)
        for _, diffs, *_ in f.reads:
            for p, c in diffs.items():
                codes[p].add(c)
        for p, cs in codes.items():
            if len(cs) == 1:
                counts[(p, next(iter(cs)))] += 1
    # Pairs covering each position (mates' spans merged), for the share test.
    cover = [0] * (hi - lo + 2)
    for f in frag_list:
        spans = sorted(span for *_, span in f.reads)
        merged = []
        for s_, e_ in spans:
            if merged and s_ <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e_)
            else:
                merged.append([s_, e_])
        for s_, e_ in merged:
            cover[s_ - lo] += 1
            cover[e_ - lo + 1] -= 1
    for i in range(1, len(cover)):
        cover[i] += cover[i - 1]

    def depth(p):
        return cover[p - lo] if lo <= p <= hi else 0

    informative = defaultdict(set)
    ranked = sorted(counts.items(), key=lambda kv: -kv[1])
    for (p, c), n in ranked:
        if n >= settings.min_reads and n >= settings.min_fraction * max(1, depth(p)):
            informative[p].add(c)
    for p, codes in settings.listed.items():
        if lo <= p <= hi:
            informative[p].update(codes)
    informative.setdefault(site, set())          # the site is always shown
    if len(informative) > MAX_POSITIONS:
        keep = {site} | {p for p in settings.listed if p in informative}
        by_count = sorted((p for p in informative if p not in keep),
                          key=lambda p: -max((counts[(p, c)] for c in informative[p]), default=0))
        for p in by_count[:max(0, MAX_POSITIONS - len(keep))]:
            keep.add(p)
        informative = defaultdict(set, {p: informative[p] for p in keep})
    positions = sorted(informative)

    # 3. States, plus a last element: "N" for a noisy pair, "C" for a clean one.
    def other_edits(f):
        n = sum(1 for _, diffs, *_ in f.reads for p, c in diffs.items() if c not in informative.get(p, ()))
        return n / max(1, len(f.reads))

    def states_of(f):
        out = []
        for p in positions:
            s = f.state(p)
            if s not in (UNKNOWN, "R", "*") and s not in informative[p]:
                s = "x"
            out.append(s)
        out.append("N" if other_edits(f) >= NOISY_EDITS else "C")
        return tuple(out)

    keys = defaultdict(list)
    for f in frag_list:
        keys[states_of(f)].append(f)

    # 4. Patterns: full keys first, partial keys join the one pattern they agree with.
    groups = []          # [{"states": list, "frags": list}]

    def compatible(g, key):
        return all(k == UNKNOWN or g["states"][i] in (UNKNOWN, k) for i, k in enumerate(key))

    pending = []
    order = sorted(keys.items(), key=lambda kv: (-sum(k != UNKNOWN for k in kv[0][:-1]), -len(kv[1])))
    for key, fl in order:
        if all(k == UNKNOWN for k in key[:-1]):
            pending.append((key, fl))
            continue
        match = [g for g in groups if compatible(g, key)]
        if len(match) == 1:
            g = match[0]
            g["states"] = [k if k != UNKNOWN else s for k, s in zip(key, g["states"])]
            g["frags"].extend(fl)
        elif not match:
            groups.append({"states": list(key), "frags": list(fl)})
        else:
            pending.append((key, fl))
    # Second chance for the ambiguous ones, now that the patterns are known.
    uninformative = []
    for key, fl in pending:
        match = [g for g in groups if compatible(g, key)]
        if len(match) == 1:
            match[0]["frags"].extend(fl)
        else:
            uninformative.append((key, fl, match))

    groups.sort(key=lambda g: -len(g["frags"]))
    assigned = sum(len(g["frags"]) for g in groups)
    patterns = []
    assignment = {}
    for i, g in enumerate(groups, start=1):
        pid = f"p{i}"
        n = len(g["frags"])
        share = n / assigned if assigned else 0.0
        common = share >= settings.common_fraction and n >= settings.common_reads
        g["id"], g["common"] = pid, common
        for f in g["frags"]:
            assignment[f.name] = pid
        patterns.append({
            "id": pid, "n": n, "share": round(share, 4), "common": common,
            "states": {str(p): s for p, s in zip(positions, g["states"])},
            "noisy": g["states"][-1] == "N",
            "evidence": _evidence(g["frags"], positions, g["states"][:-1], informative, primers, chrom, site),
        })
    # Uninformative pairs follow the patterns they agree with: common if any is.
    unin_common = unin_rare = 0
    for key, fl, match in uninformative:
        is_common = any(g.get("common") for g in match) or not match
        for f in fl:
            assignment[f.name] = "uc" if is_common else "ur"
        if is_common:
            unin_common += len(fl)
        else:
            unin_rare += len(fl)

    _mark_mixed(patterns, positions)

    pos_info = []
    for p in positions:
        ref_base = ref_seq[p - lo] if lo <= p <= hi else "N"
        pos_info.append({
            "pos": p, "ref": ref_base, "codes": sorted(informative[p]),
            "listed": p in settings.listed, "site": p == site,
            "counts": {c: counts[(p, c)] for c in informative[p]},
            "depth": depth(p),
        })
    return {
        "positions": pos_info, "patterns": patterns, "total": len(frag_list),
        "assigned": assigned, "uninformative": {"common": unin_common, "rare": unin_rare},
        "settings": {"min_reads": settings.min_reads, "min_fraction": settings.min_fraction,
                     "min_bq": settings.min_bq, "common_fraction": settings.common_fraction,
                     "common_reads": settings.common_reads},
        "assignment": assignment,
    }


def _amplicon_of(f: Fragment, primers, chrom):
    """Amplicon name of a pair from the primers its ends fall in ("" if none)."""
    if primers is None or not primers.primers.get(chrom):
        return ""
    s = min(span[0] for *_, span in f.reads)
    e = max(span[1] for *_, span in f.reads)
    slack = 3
    for ps, pe, name, strand in primers.primers[chrom]:
        if strand in ("+", "") and ps - slack <= s <= pe:
            return _amplicon_name(name)
        if strand in ("-", "") and ps <= e <= pe + slack:
            return _amplicon_name(name)
    return ""


def _amplicon_name(primer_name):
    n = primer_name
    for suffix in ("_LEFT", "_RIGHT", "_left", "_right"):
        i = n.find(suffix)
        if i > 0:
            return n[:i]
    return n


def _evidence(frags, positions, states, informative, primers, chrom, site) -> dict:
    """
    Read evidence of a pattern, over the reads that show its differences (or,
    for a pattern without differences, the reads covering the site).
    """
    diff_pos = [p for p, s in zip(positions, states) if s not in ("R", UNKNOWN)] or [site]
    fwd = rev = near = shown = 0
    other_edits, mapqs, inserts = [], [], []
    dup = clipped = 0
    amplicons = Counter()
    for f in frags:
        reads = [t for t in f.reads]
        if any(t[0].is_duplicate for t in reads):
            dup += 1
        if any(_soft_clip(t[0]) >= CLIP_MIN_BP for t in reads):
            clipped += 1
        tl = [abs(t[0].template_length) for t in reads if t[0].template_length]
        if tl:
            inserts.append(tl[0])
        amp = _amplicon_of(f, primers, chrom)
        if amp:
            amplicons[amp] += 1
        n_other = 0
        for read, diffs, low, inside, (s, e) in reads:
            mapqs.append(read.mapping_quality)
            n_other += sum(1 for p, c in diffs.items() if c not in informative.get(p, ()))
            carried = [p for p in diff_pos if s <= p <= e and (p in diffs or diff_pos == [site])]
            if carried:
                shown += 1
                if read.is_reverse:
                    rev += 1
                else:
                    fwd += 1
                d = min(min(p - s, e - p) for p in carried)
                if d < NEAR_END_BP:
                    near += 1
        other_edits.append(n_other / max(1, len(reads)))
    n = len(frags)
    inserts.sort()
    return {
        "reads": shown, "fwd": fwd, "rev": rev, "near_end": near,
        "other_edits": round(sum(other_edits) / n, 2) if n else None,
        "mapq": round(sum(mapqs) / len(mapqs), 1) if mapqs else None,
        "dup": dup, "clipped": clipped,
        "insert": inserts[len(inserts) // 2] if inserts else None,
        "amplicons": dict(amplicons.most_common(4)),
    }


def _mark_mixed(patterns, positions):
    """
    A rare pattern whose states match one common pattern on the left of some
    point and another on the right looks like a mix of the two: a PCR chimera
    (or a recombinant). Marks it with "mixed": [left id, right id].
    """
    common = [p for p in patterns if p["common"]]
    if len(common) < 2:
        return
    keys = [str(p) for p in positions]
    for pat in patterns:
        if pat["common"]:
            continue
        st = [pat["states"][k] for k in keys]
        if pat.get("noisy"):
            continue
        known = [i for i, s in enumerate(st) if s != UNKNOWN]
        if len(known) < 2:
            continue
        found = None
        for a in common:
            for b in common:
                if a is b:
                    continue
                sa = [a["states"][k] for k in keys]
                sb = [b["states"][k] for k in keys]
                for cut in range(1, len(known)):
                    left, right = known[:cut], known[cut:]
                    # Same as a on the left and b on the right, on positions both
                    # patterns have reads at, and really different from each.
                    if all(sa[i] == st[i] for i in left) and all(sb[i] == st[i] for i in right):
                        if any(sa[i] not in (UNKNOWN, st[i]) for i in right) and any(sb[i] not in (UNKNOWN, st[i]) for i in left):
                            found = (a["id"], b["id"])
                            break
                if found:
                    break
            if found:
                break
        if found:
            pat["mixed"] = list(found)


# ---- caching -------------------------------------------------------------
def cached_analysis(site_bam: Path, fasta_path, chrom, site, settings: Settings, primers=None) -> dict:
    """analyse(), cached next to the site BAM, keyed by the settings and primers."""
    site_bam = Path(site_bam)
    pkey = ""
    if primers is not None and getattr(primers, "primer_path", None):
        st = primers.primer_path.stat()
        pkey = f"{st.st_mtime}:{st.st_size}"
    key = hashlib.sha1(f"{settings.key()}|{pkey}|{site_bam.stat().st_mtime}".encode()).hexdigest()[:12]
    path = site_bam.with_name(f"{site_bam.name}.patterns.{key}.json")
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            pass
    result = analyse(site_bam, fasta_path, chrom, int(site), settings, primers)
    tmp = path.with_suffix(f".{id(result)}.tmp")
    tmp.write_text(json.dumps(result), encoding="utf-8")
    tmp.replace(path)
    return result


def keep_names(result: dict, mode: str) -> Optional[set]:
    """
    Read names to keep for a pattern filter: "common", "rare", or pattern
    ids ("p1|p3"). None means no filtering.
    """
    if not mode:
        return None
    assignment = result.get("assignment", {})
    common = {p["id"] for p in result.get("patterns", []) if p["common"]}
    if mode == "common":
        return {n for n, pid in assignment.items() if pid in common or pid == "uc"}
    if mode == "rare":
        return {n for n, pid in assignment.items() if (pid not in common and pid != "uc")}
    wanted = set(mode.split("|"))
    return {n for n, pid in assignment.items() if pid in wanted}
