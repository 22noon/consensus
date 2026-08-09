import pysam


def make_bam_read(
    qname,
    is_read1,
    start,
    mate_start,
    reverse=False,
    mate_reverse=True,
    sa_tag=None,
):
    r = pysam.AlignedSegment()

    r.query_name = qname
    r.query_sequence = "AAAAAAAAAA"
    r.query_qualities = pysam.qualitystring_to_array("FFFFFFFFFF")

    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = "10M"

    r.is_paired = True
    r.is_proper_pair = True

    r.is_read1 = is_read1
    r.is_read2 = not is_read1

    r.is_reverse = reverse
    r.mate_is_reverse = mate_reverse

    r.next_reference_id = 0
    r.next_reference_start = mate_start

    if sa_tag:
        r.set_tag("SA", sa_tag)

    return r


def write_and_index_bam(path, reads):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }

    reads = sorted(reads, key=lambda x: x.reference_start)

    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for r in reads:
            bam.write(r)

    pysam.index(str(path))
