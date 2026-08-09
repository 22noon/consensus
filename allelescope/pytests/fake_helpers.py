"""
Minimal, deterministic stand-ins for the project's real helper functions
(Tag_bits, tag_read, calculate_softclip_lengths, parse_sa_tag,
get_allele_at_base, get_indels_at_base, summarise_indels).

These exist ONLY so extract_mates() can be exercised in isolation by the
test suite. Swap extract_mates.py's dependency wiring for the real
project implementations in production -- do not use this module there.
"""

Tag_bits = {
    'NON_SPANNING_MATE': 0x1,
    'F1R2':  0x2,
    'R2F1':  0x4,
    'F2R1':  0x8,
    'R1F2':  0x10,
    'F1F2':  0x20,
    'R1R2':  0x40,
    'SPLIT':   0x80,
    'SPLITX':  0x100,
    'OVERLAP': 0x200,
    'IMPROPER': 0x400,
    'UNPAIRED': 0x800,
    'PROPER':   0x1000,
}


def tag_read(read, XO):
    read.set_tag('XO', XO, value_type='i')
    return read


def calculate_softclip_lengths(read):
    # No-op stand-in; the real implementation adds soft-clip-length tags.
    return read


def parse_sa_tag(sa_string):
    """Parse a standard SAM 'SA:Z:' tag value into a list of dicts."""
    entries = []
    for part in sa_string.strip(';').split(';'):
        if not part:
            continue
        rname, pos, strand, cigar, mapq, nm = part.split(',')
        entries.append({
            'rname': rname,
            'pos': int(pos),
            'strand': strand,
            'cigar': cigar,
            'mapq': int(mapq),
            'nm': int(nm),
        })
    return entries


def get_allele_at_base(read, target_pos):
    """Return the query base aligned to 0-based reference target_pos, or None."""
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == target_pos:
            return read.query_sequence[qpos]
    return None


def get_indels_at_base(read, target_pos):
    # Indel detection isn't exercised by these tests; real implementation
    # inspects the CIGAR around target_pos.
    return []


def summarise_indels(indels):
    return 0, "", 0, 0
