"""
Read filtering for the per-site BAMs served to the viewer.

The viewer sends the filter settings as query parameters; this module turns
them into a predicate and applies it to a per-site BAM. It has no Flask or
filesystem-layout dependencies so it can be tested on its own.

Numeric limits
--------------
Every numeric limit is optional. An empty / missing value means "no limit",
so nothing is filtered by default. A value that is present is applied as
given, including 0 (e.g. a maximum edit distance of 0 keeps only reads that
match the reference exactly).

    minMapQ        keep reads with MAPQ >= value
    SoftClip       keep reads with at least this many soft-clipped bases (SL + SR)
    SoftClipF      keep reads with at most this many soft-clipped bases
    EditDistance   keep reads with NM >= value
    EditDistanceF  keep reads with NM <= value
    BAQ            keep reads whose BAQ at the site is >= value
    BAQF           keep reads whose BAQ at the site is <= value

The BAQ limits are evaluated on the reads that cover the site (see
``baq_filter.baq_at_site``). A read pair is kept or dropped as a whole, so
that mates and split partners, which have no base at the site, follow the
read that does.
"""

from dataclasses import dataclass
from typing import Optional

import pysam

from .baq_filter import baq_at_site


def optional_int(value) -> Optional[int]:
    """Parse an optional non-negative integer limit; None means 'no limit'."""
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in ("undefined", "null", "nan", "none"):
        return None
    try:
        number = int(float(text))
    except ValueError:
        return None
    return number if number >= 0 else None


def _flag_int(value) -> int:
    number = optional_int(value)
    return number or 0


@dataclass
class ReadFilter:
    flag_require: int = 0          # SAM flag bits that must all be set
    flag_exclude: int = 0          # SAM flag bits that must all be unset
    xo_require: int = 0            # XO bits that must all be set
    xo_exclude: int = 0            # XO bits that must all be unset
    min_mapq: Optional[int] = None
    min_soft_clip: Optional[int] = None
    max_soft_clip: Optional[int] = None
    min_nm: Optional[int] = None
    max_nm: Optional[int] = None
    min_baq: Optional[int] = None
    max_baq: Optional[int] = None
    alleles: Optional[frozenset] = None   # allowed XA values
    site_pos0: Optional[int] = None       # 0-based site position, for BAQ

    @classmethod
    def from_params(cls, params) -> "ReadFilter":
        xa = (params.get("XAFilter") or "").strip()
        pos = optional_int(params.get("Pos"))
        return cls(
            flag_require  = _flag_int(params.get("Flagf")),
            flag_exclude  = _flag_int(params.get("FlagF")),
            xo_require    = _flag_int(params.get("Tagf")),
            xo_exclude    = _flag_int(params.get("TagF")),
            min_mapq      = optional_int(params.get("minMapQ")),
            min_soft_clip = optional_int(params.get("SoftClip")),
            max_soft_clip = optional_int(params.get("SoftClipF")),
            min_nm        = optional_int(params.get("EditDistance")),
            max_nm        = optional_int(params.get("EditDistanceF")),
            min_baq       = optional_int(params.get("BAQ")),
            max_baq       = optional_int(params.get("BAQF")),
            alleles       = frozenset(xa.split("|")) if xa else None,
            site_pos0     = pos - 1 if pos else None,
        )

    @property
    def uses_baq(self) -> bool:
        return (self.min_baq is not None or self.max_baq is not None) and self.site_pos0 is not None

    def baq_fails(self, read) -> bool:
        """True if the read covers the site and its BAQ there is outside the limits."""
        baq = baq_at_site(read, self.site_pos0)
        if baq is None:
            return False
        if self.min_baq is not None and baq < self.min_baq:
            return True
        if self.max_baq is not None and baq > self.max_baq:
            return True
        return False

    def keeps(self, read) -> bool:
        """Per-record tests (everything except BAQ, which is decided per read name)."""
        if self.flag_require and (read.flag & self.flag_require) != self.flag_require:
            return False
        if self.flag_exclude and (read.flag & self.flag_exclude):
            return False

        if self.xo_require or self.xo_exclude:
            xo = read.get_tag("XO") if read.has_tag("XO") else 0
            if self.xo_require and (xo & self.xo_require) != self.xo_require:
                return False
            if self.xo_exclude and (xo & self.xo_exclude):
                return False

        if self.min_mapq is not None and read.mapping_quality < self.min_mapq:
            return False

        if self.min_nm is not None or self.max_nm is not None:
            nm = read.get_tag("NM") if read.has_tag("NM") else 0
            if self.min_nm is not None and nm < self.min_nm:
                return False
            if self.max_nm is not None and nm > self.max_nm:
                return False

        if self.alleles is not None:
            # Reads without an allele call cannot match a selected allele.
            if not read.has_tag("XA") or read.get_tag("XA") not in self.alleles:
                return False

        if self.min_soft_clip is not None or self.max_soft_clip is not None:
            clip = soft_clip_length(read)
            if self.min_soft_clip is not None and clip < self.min_soft_clip:
                return False
            if self.max_soft_clip is not None and clip > self.max_soft_clip:
                return False

        return True


def soft_clip_length(read) -> int:
    """Total soft-clipped bases, from the SL/SR tags when present, else the CIGAR."""
    if read.has_tag("SL") and read.has_tag("SR"):
        return read.get_tag("SL") + read.get_tag("SR")
    cigar = read.cigartuples or []
    left = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
    right = cigar[-1][1] if len(cigar) > 1 and cigar[-1][0] == 4 else 0
    return left + right


def filter_reads(bam_path, params):
    """Return (header, reads) for the reads in ``bam_path`` that pass ``params``."""
    read_filter = ReadFilter.from_params(params)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        header = bam.header
        records = list(bam)

    # BAQ is decided per read name over the whole read set, independently of
    # the other filters, so a pair is treated the same whichever of its
    # records the other filters happen to keep.
    failing_baq = set()
    if read_filter.uses_baq:
        failing_baq = {read.query_name for read in records if read_filter.baq_fails(read)}

    reads = [
        read for read in records
        if read.query_name not in failing_baq and read_filter.keeps(read)
    ]
    return header, reads
