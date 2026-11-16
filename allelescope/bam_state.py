"""
Work out which preprocessing steps a BAM still needs.

build_variant_browser normally runs every BAM through

    samtools collate | samtools fixmate -m | samtools sort | samtools markdup

Most alignment pipelines already deliver a coordinate-sorted,
duplicate-marked BAM, and redoing all four passes on a whole-genome file
takes a long time. Each step leaves a trace that can be checked cheaply:

    collate    nothing; it only groups reads by name so fixmate can see
               both mates. It is needed only when fixmate has to run.
    fixmate -m MC (mate CIGAR) and ms (mate score) tags on paired reads,
               which samtools markdup requires.
    sort       SO:coordinate in the @HD header line (and records in order).
    markdup    an @PG header line from a duplicate-marking program.

Only the header and the first records are read, so this takes well under a
second for any file size.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import pysam

# Programs that mark (or remove) duplicates, matched against the ID, PN and
# CL fields of @PG lines (case-insensitive).
DUPLICATE_MARKERS = [
    (re.compile(r"samtools\b.*\bmarkdup\b|^samtools\.?markdup", re.I), "samtools markdup"),
    (re.compile(r"\bmarkduplicates(spark|withmatecigar)?\b", re.I), "Picard/GATK MarkDuplicates"),
    (re.compile(r"sambamba\b.*\bmarkdup\b|^sambamba[_ .-]?markdup", re.I), "sambamba markdup"),
    (re.compile(r"\bbam(sor)?madup\b|\bbammarkduplicates2?\b", re.I), "biobambam"),
    (re.compile(r"samtools\b.*\brmdup\b", re.I), "samtools rmdup"),
    (re.compile(r"\bsamblaster\b", re.I), "samblaster"),
    (re.compile(r"\bumi_tools\b.*\bdedup\b|^umi_tools", re.I), "umi_tools dedup"),
]

PLAN_FULL = "full"        # collate -> fixmate -> sort -> markdup
PLAN_MARKDUP = "markdup"  # markdup only (sorted, mate tags present)
PLAN_SORT = "sort"        # sort only (duplicates already marked)
PLAN_COPY = "copy"        # copy (or link) and index


@dataclass
class BamState:
    header_sorted: bool = False          # @HD SO:coordinate
    records_in_order: bool = True        # sampled records really are in order
    duplicate_tool: Optional[str] = None  # tool found in @PG, if any
    sampled: int = 0                     # records read
    paired_checked: int = 0              # paired, mapped, mate mapped, primary
    with_mate_tags: int = 0              # ... of those carrying MC and ms
    duplicates_seen: int = 0             # records with the 0x400 flag
    notes: List[str] = field(default_factory=list)

    @property
    def is_sorted(self) -> bool:
        return self.header_sorted and self.records_in_order

    @property
    def duplicates_marked(self) -> bool:
        return self.duplicate_tool is not None

    @property
    def paired(self) -> bool:
        return self.paired_checked > 0

    @property
    def has_mate_tags(self) -> bool:
        return self.paired and self.with_mate_tags == self.paired_checked

    def plan(self) -> str:
        if self.duplicates_marked:
            return PLAN_COPY if self.is_sorted else PLAN_SORT
        if self.is_sorted and self.has_mate_tags:
            return PLAN_MARKDUP
        return PLAN_FULL

    def describe(self) -> List[str]:
        """Human-readable findings, one per line."""
        lines = [
            "sorted by coordinate: " + (
                "yes" if self.is_sorted else
                "header says so, but records are out of order" if self.header_sorted else "no"),
            "duplicates marked: " + (
                f"yes ({self.duplicate_tool} in @PG)" if self.duplicates_marked else
                "no duplicate-marking program in @PG" + (
                    f" ({self.duplicates_seen} flagged records seen, but the tool is unknown)"
                    if self.duplicates_seen else "")),
        ]
        if self.paired:
            lines.append(f"fixmate -m tags (MC, ms): {self.with_mate_tags} of {self.paired_checked} sampled pairs")
        else:
            lines.append("paired reads: none in the first records")
        return lines + self.notes


def inspect_bam(bam_path, sample_size: int = 20000) -> BamState:
    """Read the header and the first ``sample_size`` records of a BAM."""
    state = BamState()
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        header = bam.header.to_dict()
        state.header_sorted = header.get("HD", {}).get("SO") == "coordinate"

        for pg in header.get("PG", []):
            text = " ".join(str(pg.get(k, "")) for k in ("ID", "PN", "CL"))
            for pattern, name in DUPLICATE_MARKERS:
                if pattern.search(text):
                    state.duplicate_tool = name
                    break
            if state.duplicate_tool:
                break

        last = (-1, -1)
        seen_unmapped_tail = False
        for read in bam.fetch(until_eof=True):
            state.sampled += 1
            if read.is_duplicate:
                state.duplicates_seen += 1

            # Coordinate order: (reference id, start) never decreases; reads
            # without a reference come last.
            if read.reference_id < 0:
                seen_unmapped_tail = True
            else:
                key = (read.reference_id, read.reference_start)
                if seen_unmapped_tail or key < last:
                    state.records_in_order = False
                last = key

            if (read.is_paired and not read.is_unmapped and not read.mate_is_unmapped
                    and not read.is_secondary and not read.is_supplementary):
                state.paired_checked += 1
                if read.has_tag("MC") and read.has_tag("ms"):
                    state.with_mate_tags += 1

            if state.sampled >= sample_size:
                break

    if state.sampled == 0:
        state.notes.append("the BAM has no records")
    return state


def existing_index(bam_path) -> Optional[Path]:
    """An index next to the BAM that is at least as new as the BAM, if any."""
    bam_path = Path(bam_path)
    candidates = [Path(str(bam_path) + ".bai"), bam_path.with_suffix(".bai"), Path(str(bam_path) + ".csi")]
    bam_mtime = bam_path.stat().st_mtime
    for idx in candidates:
        if idx.is_file() and idx.stat().st_mtime >= bam_mtime:
            return idx
    return None
