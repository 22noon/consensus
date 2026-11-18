#!/usr/bin/env python3
import argparse
import io

from flask import Flask, send_from_directory, Response, request, jsonify,send_file, redirect 
import tempfile
import urllib.parse

import pysam
import mimetypes
import os
import pysam
import sys
import re
import hashlib
import threading
import json
from pathlib import Path

from .process_variant import (process_variant, ensure_current_bam_dir, site_stem, site_sample,
                              read_site_info, PROCESSING_VERSION)
from .read_filters import filter_reads
from .annotation import Annotation, instance_files, open_annotation
import shutil


app = Flask(__name__)

default_dir = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-dir', '-d', default=default_dir)
parser.add_argument('--port', type=int, default=8000)
args, _ = parser.parse_known_args()
DATA_DIR = Path(os.path.abspath(args.data_dir))
CACHE_ROOT = DATA_DIR / "cache"
CACHE_ROOT.mkdir(exist_ok=True)
pending_extractions = {} # stack of read extactions pending

if not DATA_DIR.is_dir():
    print(f"DATA_DIR does not exist: {DATA_DIR}", file=sys.stderr)
    sys.exit(1)

# -----------------------------
# Helpers (from complex server)
# -----------------------------
def get_mimetype(filename):
    mimetype, _ = mimetypes.guess_type(filename)
    overrides = {
        '.bam': 'application/octet-stream',
        '.bai': 'application/octet-stream',
        '.fai': 'application/octet-stream',
        '.fastq': 'text/plain',
        '.fq': 'text/plain',
        '.vcf': 'text/plain',
        '.gff': 'text/plain',
        '.bed': 'text/plain',
        '.js': 'application/javascript',
        '.css': 'text/css',
        '.html': 'text/html',
    }
    suffix = Path(filename).suffix.lower()
    return overrides.get(suffix, mimetype or 'application/octet-stream')


def send_file_with_range(path: Path, mimetype="application/octet-stream"):
    if not path.exists() or not path.is_file():
        return Response(status=404)

    file_size = path.stat().st_size
    range_header = request.headers.get("Range")

    if not range_header:
        return send_from_directory(path.parent, path.name, mimetype=mimetype)

    try:
        unit, rng = range_header.split("=")
        start_str, end_str = rng.split("-")
        if start_str == "":
            length = int(end_str)
            start = max(0, file_size - length)
            end = file_size - 1
        else:
            start = int(start_str)
            end = int(end_str) if end_str else file_size - 1
    except Exception:
        return send_from_directory(path.parent, path.name, mimetype=mimetype)

    if start >= file_size or start < 0 or end < start:
        print(f"Serving {path} with Range: {start}-{end} (size: {file_size}) ", flush=True) 
        return Response(
            status=204,  # No Content
            headers={"Accept-Ranges": "bytes", "X-Empty-File": "true"}
        )
        return send_from_directory(path.parent, path.name, mimetype=mimetype)

    end = min(end, file_size - 1)
    length = end - start + 1

    with open(path, "rb") as f:
        f.seek(start)
        data = f.read(length)


    headers = {
        "Accept-Ranges": "bytes",
        "Content-Type": mimetype,
        "Cache-Control": "no-cache",
        "Content-Length": str(length),
        "Content-Range": f"bytes {start}-{end}/{file_size}",
    }
    return Response(data, status=206, headers=headers)


def build_cache_key(*args):
    # PROCESSING_VERSION is part of the key so filtered BAMs cached by an
    # older version (different tag or filter semantics) are never reused.
    key = "|".join(str(a) for a in (PROCESSING_VERSION, *args)).encode()
    print(f"Cache key: {key}")
    return hashlib.md5(key).hexdigest()

# IGV asks for a BAM and its index at the same time, and quick filter changes
# repeat requests, so one cached file can be requested by several threads at
# once. Each cache entry is built under its own lock, into temporary files
# that are renamed into place when complete.
_cache_locks = {}
_cache_locks_guard = threading.Lock()


def _cache_lock(key):
    with _cache_locks_guard:
        return _cache_locks.setdefault(key, threading.Lock())


def processed_bam(serverpath, Chrom, Pos, Ref, sample=None):
    # One preparation per site at a time; a second request waits and then
    # finds the site prepared.
    key = f"site|{Path(serverpath).resolve()}|{Chrom}|{Pos}|{site_sample(sample)}"
    with _cache_lock(key):
        return process_variant(serverpath, Chrom, Pos, Ref, sample=sample)


def get_filtered_reads_bam(serverpath: Path, filename: str, **params):
    """
    Return (header, reads) for the per-site BAM after applying the filter
    parameters. The filter rules live in read_filters.py.
    """
    token = params.get("Token") or ""
    bam_file = ensure_current_bam_dir(serverpath) / f"{filename}{token}.calmd.bam"
    print(f"Getting filtered reads from BAM for {filename} with params {params} from {bam_file}", flush=True)
    if not bam_file.exists():
        processed_bam(serverpath, params.get("Chrom"), params.get("Pos"), params.get("Ref"),
                      sample=params.get("Sample"))

    return filter_reads(bam_file, params)

def build_filtered_bam(serverpath, cache_dir, filename, cache_key, params):
    bam_path = cache_dir / f"{cache_key}.bam"
    bai_path = cache_dir / f"{cache_key}.bam.bai"

    if bam_path.exists() and bai_path.exists():
        return bam_path

    with _cache_lock(str(bam_path)):
        if bam_path.exists() and bai_path.exists():
            return bam_path
        tmp_bam = cache_dir / f"{cache_key}.{threading.get_ident()}.tmp.bam"
        tmp_bai = Path(f"{tmp_bam}.bai")
        try:
            header, reads = get_filtered_reads_bam(serverpath, filename, **params)
            with pysam.AlignmentFile(str(tmp_bam), "wb", header=header) as out:
                for read in reads:
                    out.write(read)
            pysam.index(str(tmp_bam), str(tmp_bai))
            tmp_bam.replace(bam_path)
            tmp_bai.replace(bai_path)
        finally:
            tmp_bam.unlink(missing_ok=True)
            tmp_bai.unlink(missing_ok=True)
    return bam_path


# -----------------------------
# BAM/BAI routes (from complex server)
# -----------------------------

@app.route("/api/bai")
@app.route("/<path:browser>/api/bai")
def serve_bai(browser=None):
    if browser is None:
        serverpath = DATA_DIR          # default location
        cache_dir = CACHE_ROOT 
    else:
        serverpath = DATA_DIR / browser
        cache_dir = CACHE_ROOT / browser
    cache_dir.mkdir(parents=True, exist_ok=True)
    print("SERVE_BAI HIT:", browser,serverpath,cache_dir, flush=True)

    Chrom = request.args.get("Chrom")
    Pos = request.args.get("Pos")
    filename = site_stem(Chrom, Pos, request.args.get("Sample"))
    params = request.args.to_dict()
    cache_key = build_cache_key(filename, *params.values())
    print(f"Cache key for BAI: {cache_key}", flush=True)

    bam_path = cache_dir / f"{cache_key}.bam"
    bai_path = cache_dir / f"{cache_key}.bam.bai"

    if not bam_path.exists() or not bai_path.exists():
        build_filtered_bam(serverpath, cache_dir, filename, cache_key, params)

    return send_file_with_range(bai_path)


@app.route("/api/bam")
@app.route("/<path:browser>/api/bam")
def serve_bam(browser=None):
    if browser is None:
        serverpath = DATA_DIR          # default location
        cache_dir = CACHE_ROOT 
    else:
        serverpath = DATA_DIR / browser
        cache_dir = CACHE_ROOT / browser

    print("SERVE_BAM HIT:", browser,serverpath,cache_dir, flush=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    Chrom = request.args.get("Chrom")
    Pos = request.args.get("Pos")
    Ref = request.args.get("Ref")
    Sample = request.args.get("Sample")
    filename = site_stem(Chrom, Pos, Sample)

    bam_file = ensure_current_bam_dir(serverpath) / f"{filename}.calmd.bam"
    if not bam_file.exists():
        processed_bam(serverpath, Chrom, Pos, Ref, sample=Sample)

    params = request.args.to_dict()
    cache_key = build_cache_key(filename, *params.values())
    print(f"Cache key for BAM: {cache_key}", flush=True)
    bam_path = build_filtered_bam(serverpath, cache_dir, filename, cache_key, params)

    return send_file_with_range(bam_path)

@app.route("/api/isolate_reads", methods=["POST"])
@app.route("/<path:browser>/api/isolate_reads", methods=["POST"])
def isolate_reads(browser=None):
    if browser is None:
        serverpath = DATA_DIR          # default location
    else:
        serverpath = DATA_DIR / browser

    data = request.get_json()

    target_reads = set(data.get("reads"))
    filter_string = data.get("filterString")
    token = hashlib.sha256( "\n".join(sorted(target_reads)).encode()).hexdigest()[:16]

    params = dict(urllib.parse.parse_qsl(filter_string.lstrip("?"), keep_blank_values=True))
    chrom = params.get("Chrom")
    pos   = params.get("Pos")
    print(f"ISOLATE_READS HIT: {browser}, {chrom}, {pos}, {filter_string}, token={token}", flush=True)

    stem = site_stem(chrom, pos, params.get("Sample"))
    out_bam = serverpath / "bam" / f"{stem}{token}.calmd.bam"
    in_bam = serverpath / "bam" / f"{stem}.calmd.bam"
    print("EXTRACT_READS HIT:", browser, chrom, pos, filter_string, in_bam, out_bam, token, flush=True)

    if not in_bam or not os.path.exists(in_bam):
        return jsonify({"error": "No BAM found to extract reads from"}), 404

    with pysam.AlignmentFile(in_bam, "rb") as inbam:
        with pysam.AlignmentFile(out_bam, "wb", header=inbam.header) as outbam:
            for read in inbam.fetch(chrom):
                if read.query_name in target_reads:
                    outbam.write(read)

    pysam.index(str(out_bam))
    return(jsonify({"success": True, "token": token}))

@app.route("/api/remove_reads", methods=["POST"])
@app.route("/<path:browser>/api/remove_reads", methods=["POST"])
def remove_reads(browser=None):
    if browser is None:
        serverpath = DATA_DIR          # default location
    else:
        serverpath = DATA_DIR / browser

    data = request.get_json()

    target_reads = set(data.get("reads"))
    filter_string = data.get("filterString")
    token = "X" + hashlib.sha256( "\n".join(sorted(target_reads)).encode()).hexdigest()[:16]

    params = dict(urllib.parse.parse_qsl(filter_string.lstrip("?"), keep_blank_values=True))
    chrom = params.get("Chrom")
    pos   = params.get("Pos")
    print(f"REMOVE_READS HIT: {browser}, {chrom}, {pos}, {filter_string}, token={token}", flush=True)

    stem = site_stem(chrom, pos, params.get("Sample"))
    out_bam = serverpath / "bam" / f"{stem}{token}.calmd.bam"
    in_bam = serverpath / "bam" / f"{stem}.calmd.bam"
    print("REMOVE_READS HIT:", browser, chrom, pos, filter_string, in_bam, out_bam, token, flush=True)

    if not in_bam or not os.path.exists(in_bam):
        return jsonify({"error": "No BAM found to extract reads from"}), 404

    with pysam.AlignmentFile(in_bam, "rb") as inbam:
        with pysam.AlignmentFile(out_bam, "wb", header=inbam.header) as outbam:
            for read in inbam.fetch(chrom):
                if read.query_name not in target_reads:
                    outbam.write(read)

    pysam.index(str(out_bam))
    return(jsonify({"success": True, "token": token}))

@app.route("/api/extract_reads", methods=["POST"])
@app.route("/<path:browser>/api/extract_reads", methods=["POST"])
def stage_extract_reads(browser=None):
    data = request.get_json()
    token = hashlib.sha256( "\n".join(sorted(data.get("reads"))).encode()).hexdigest()[:16]
    pending_extractions[token] = {
        "reads": data.get("reads"),
        "filterString": data.get("filterString"),
        "browser": browser
    }
    return jsonify({"token": token})

@app.route("/api/extract_reads/<token>", methods=["GET"])
@app.route("/<path:browser>/api/extract_reads/<token>", methods=["GET"])
def extract_reads(browser=None, token=None):
    pending = pending_extractions.pop(token, None)
    if not pending:
        return jsonify({"error": "Invalid or expired token"}), 404

    target_reads = set(pending["reads"])
    filter_string = pending["filterString"]
    SIZE_THRESHOLD = 50 * 1024 * 1024  # 50MB

    if not target_reads:
        return jsonify({"error": "No reads specified"}), 400

    # Parse filter string into an ordered dict, same as request.args.to_dict()
    params = dict(urllib.parse.parse_qsl(filter_string.lstrip("?"), keep_blank_values=True))

    chrom = params.get("Chrom")
    pos   = params.get("Pos")

    # Mirror serve_bam exactly
    filename  = site_stem(chrom, pos, params.get("Sample"))
    cache_key = build_cache_key(filename, *params.values())
    cache_dir = CACHE_ROOT / browser
    bam_path = cache_dir / f"{cache_key}.bam"
    print("EXTRACT_READS HIT:", browser, chrom, pos, filter_string, bam_path, token, flush=True)

    if not bam_path or not os.path.exists(bam_path):
        return jsonify({"error": "No cached BAM found — navigate to a variant first"}), 404

    with tempfile.TemporaryDirectory() as tmpdir:
        extracted_path = os.path.join(tmpdir, "extracted.bam")

        with pysam.AlignmentFile(bam_path, "rb") as inbam:
            with pysam.AlignmentFile(extracted_path, "wb", header=inbam.header) as outbam:
                for read in inbam.fetch(chrom):
                    if read.query_name in target_reads:
                        outbam.write(read)

        if os.path.getsize(extracted_path) > SIZE_THRESHOLD:
            # Small file — safe in memory before tmpdir cleanup
            with open(extracted_path, "rb") as f:
                bam_data = io.BytesIO(f.read())

            return send_file(
                bam_data,
                mimetype="application/octet-stream",
                as_attachment=True,
                download_name=f"selected_reads_{chrom}_{pos}.bam"
            )
        else:
            # Large file — copy outside tmpdir before it's deleted
            final_path = CACHE_ROOT / browser / f"extract_{chrom}_{pos}_{token}.bam"
            shutil.copy(extracted_path, final_path)

    # tmpdir safely deleted here, final_path survives

    # Large file — serve from persistent location, cron will clean up
    return send_file(
        final_path,
        mimetype="application/octet-stream",
        as_attachment=True,
        download_name=f"selected_reads_{chrom}_{pos}.bam"
    )

# ---- Gene annotation and primers -----------------------------------------
_annotations = {}                  # instance folder -> (signature, Annotation or None)
_annotations_guard = threading.Lock()


def instance_annotation(serverpath):
    """
    The Annotation of an instance (annotation.gff3 / .gtf and primers.bed in
    its folder), or None when it has neither. Re-read when a file changes, so
    files can be added to an instance without rebuilding it.
    """
    serverpath = Path(serverpath)
    fasta, gff, bed = instance_files(serverpath)
    if fasta is None or (gff is None and bed is None):
        return None
    sig = tuple((str(f), f.stat().st_mtime, f.stat().st_size) for f in (fasta, gff, bed) if f)
    with _annotations_guard:
        cached = _annotations.get(serverpath)
        if cached and cached[0] == sig:
            return cached[1]
        try:
            ann = Annotation(fasta, open_annotation(gff) if gff else None, bed)
        except Exception as exc:                      # a broken file must not break the viewer
            print(f"Could not read the annotation of {serverpath}: {exc}", flush=True)
            ann = None
        _annotations[serverpath] = (sig, ann)
        return ann


def _effect_text(e):
    """geneA:D2H:missense, as in INFO/EFFECT."""
    return f"{e['gene']}:{e['change']}:{e['type']}".replace(" ", "_").replace(",", "+").replace(";", "_")


@app.route("/api/annotate", methods=["POST"])
@app.route("/<path:browser>/api/annotate", methods=["POST"])
def annotate(browser=None):
    """
    Coding effects and primers for a list of variants.
    Body: {"variants": [{"chrom", "pos", "ref", "alt"}, ...]}   (alt may be "A,T")
    Returns {"info": {...}, "variants": {"chrom:pos": {"effects": [...], "primers": [...]}}}
    or {"info": null} when the instance has no annotation or primers.
    """
    serverpath = DATA_DIR if browser is None else DATA_DIR / browser
    ann = instance_annotation(serverpath)
    if ann is None:
        return jsonify({"info": None, "variants": {}})
    body = request.get_json(silent=True) or {}
    out = {}
    for v in body.get("variants") or []:
        try:
            chrom, pos = str(v["chrom"]), int(v["pos"])
        except (KeyError, TypeError, ValueError):
            continue
        ref = str(v.get("ref") or "")
        effects = []
        if ann.has_cds(chrom):
            for alt in str(v.get("alt") or "").split(","):
                for e in ann.effects(chrom, pos, ref, alt.strip()):
                    effects.append({**e, "alt": alt.strip()})
        entry = {"effects": effects, "coding": ann.has_cds(chrom)}
        # A variant overlaps a primer if any of its reference bases do.
        names = {p[2] for i in range(max(1, len(ref))) for p in ann.primers_at(chrom, pos + i)}
        if names:
            entry["primers"] = sorted(names)
        out[f"{chrom}:{pos}"] = entry
    return jsonify({"info": ann.info(), "variants": out})


@app.route("/api/extract")
@app.route("/<path:browser>/api/extract")
def extract_bam(browser=None):

    if browser is None:
        serverpath = DATA_DIR          # default location
    else:
        serverpath = DATA_DIR / browser
    print("EXTRACT HIT:", browser,serverpath, flush=True)
    Chrom = request.args.get("chrom")
    Pos = request.args.get("pos")
    Ref = request.args.get("ref")
    Sample = site_sample(request.args.get("sample"))

    mock_alleles = []

    result = processed_bam(serverpath, Chrom, Pos, Ref, sample=Sample)
    if result == 0:
        stat_file = serverpath / "bam" / f"{site_stem(Chrom, Pos, Sample)}.bam.indelstats.txt"
        evidence = {}
        evidence_file = Path(f"{stat_file}".replace(".indelstats.txt", ".evidence.json"))
        if evidence_file.exists():
            try:
                evidence = json.loads(evidence_file.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                evidence = {}
        with open(stat_file, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                fields = line.split("\t")
                xa, count, pct, allele_type = fields[:4]
                count = int(count)
                pct = round(float(pct)*100, 2)

                allele = {
                    "xa": xa,
                    "xd": 0,
                    "xn": 0,
                    "count": count,
                    "pct": pct,
                    "type": allele_type,
                }
                if len(fields) >= 8:
                    allele["fwd"] = int(fields[4])
                    allele["rev"] = int(fields[5])
                    allele["bq"] = float(fields[6]) if fields[6] else None
                    allele["mq"] = float(fields[7]) if fields[7] else None
                if xa in evidence:
                    allele["ev"] = evidence[xa]
                mock_alleles.append(allele)

        result = {"success": True, "alleles": mock_alleles,
                  "sampling": read_site_info(serverpath, Chrom, Pos, Sample)}
        ann = instance_annotation(serverpath)
        if ann is not None:
            pos = int(Pos)
            result["primers"] = [p[2] for p in ann.primers_at(Chrom, pos)]
            result["coding"] = ann.has_cds(Chrom)
            for allele in mock_alleles:
                if ann.has_cds(Chrom):
                    allele["effects"] = ann.allele_effects(Chrom, pos, allele["xa"])
                ev = allele.get("ev")
                if result["primers"] and ev and "starts" in ev:
                    ev["primer"] = ann.primer_reads(Chrom, pos, ev["starts"], ev["ends"])
                if ev:                                  # not needed in the browser
                    ev.pop("starts", None)
                    ev.pop("ends", None)
        else:
            for allele in mock_alleles:
                if allele.get("ev"):
                    allele["ev"].pop("starts", None)
                    allele["ev"].pop("ends", None)
        return jsonify(result)
    return jsonify({"success": False}), 500


def _vcf_text(value: str) -> str:
    """Percent-encode the characters a VCF INFO value cannot hold."""
    out = []
    for ch in value:
        if ch in "%;=,:" or ch in "\t\r\n" or ch == " ":
            out.append(f"%{ord(ch):02X}")
        else:
            out.append(ch)
    return "".join(out)


def instance_vcf(serverpath: Path):
    """The VCF an instance was built from (the bgzipped copy in its folder)."""
    candidates = sorted(p for p in Path(serverpath).glob("*.vcf.gz")
                        if not p.name.endswith(".reviewed.vcf.gz"))
    return candidates[0] if candidates else None


@app.route("/api/review_vcf", methods=["POST"])
@app.route("/<path:browser>/api/review_vcf", methods=["POST"])
def review_vcf(browser=None):
    """
    The instance's VCF with the browser's review added to each record:
        INFO/REVIEW       confirmed | doubtful | artefact
        INFO/REVIEW_NOTE  the note, percent-encoded (%3B for ";", %20 for space ...)
        INFO/READ_FLAGS   read-evidence warnings for sites opened in the browser
        INFO/EFFECT       coding effects, if the instance has an annotation
        INFO/PRIMER       primers overlapping the variant, if it has a primer scheme
    Body: {"reviews": {"chrom:pos": {"s": status, "n": note}, ...},
           "flags":   {"chrom:pos": ["read-end", "strand", ...], ...}}
    """
    serverpath = DATA_DIR if browser is None else DATA_DIR / browser
    src = instance_vcf(serverpath)
    if src is None:
        return jsonify({"error": "This instance has no VCF"}), 404
    body = request.get_json(silent=True) or {}
    reviews = body.get("reviews") or {}
    flags = body.get("flags") or {}

    with pysam.VariantFile(str(src)) as vin:
        header = vin.header.copy()
        if "REVIEW" not in header.info:
            header.add_meta("INFO", items=[("ID", "REVIEW"), ("Number", "1"), ("Type", "String"),
                                           ("Description", "Allelescope review: confirmed, doubtful or artefact")])
        if "REVIEW_NOTE" not in header.info:
            header.add_meta("INFO", items=[("ID", "REVIEW_NOTE"), ("Number", "1"), ("Type", "String"),
                                           ("Description", "Allelescope review note (percent-encoded)")])
        ann = instance_annotation(serverpath)
        if ann is not None and ann.cds and "EFFECT" not in header.info:
            header.add_meta("INFO", items=[("ID", "EFFECT"), ("Number", "."), ("Type", "String"),
                                           ("Description", "Allelescope coding effect per ALT: gene:change:type, "
                                            "for example S:D614G:missense (from the instance's annotation)")])
        if ann is not None and ann.primers and "PRIMER" not in header.info:
            header.add_meta("INFO", items=[("ID", "PRIMER"), ("Number", "."), ("Type", "String"),
                                           ("Description", "Allelescope: amplicon primers overlapping the variant")])
        if flags and "READ_FLAGS" not in header.info:
            header.add_meta("INFO", items=[("ID", "READ_FLAGS"), ("Number", "."), ("Type", "String"),
                                           ("Description", "Allelescope read-evidence warnings: strand, orientation, "
                                            "read-end, mismatches, ambiguous, clipped, split, insert-size, "
                                            "duplicates, low-bq, af-mismatch, primer, primer-masked")])
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "reviewed.vcf"
            with pysam.VariantFile(str(out_path), "w", header=header) as vout:
                for rec in vin:
                    rec.translate(vout.header)
                    r = reviews.get(f"{rec.chrom}:{rec.pos}") or {}
                    status = r.get("s")
                    note = (r.get("n") or "").strip()
                    if status in ("confirmed", "doubtful", "artefact"):
                        rec.info["REVIEW"] = status
                    if note:
                        rec.info["REVIEW_NOTE"] = _vcf_text(note)
                    rf = [f for f in (flags.get(f"{rec.chrom}:{rec.pos}") or [])
                          if isinstance(f, str) and re.fullmatch(r"[a-z-]+", f)]
                    if rf:
                        rec.info["READ_FLAGS"] = rf
                    if ann is not None:
                        effects = [_effect_text(e) for alt in (rec.alts or ())
                                   for e in ann.effects(rec.chrom, rec.pos, rec.ref, alt)]
                        if effects and ann.cds:
                            rec.info["EFFECT"] = effects
                        names = sorted({p[2] for i in range(len(rec.ref))
                                        for p in ann.primers_at(rec.chrom, rec.pos + i)})
                        if names and ann.primers:
                            rec.info["PRIMER"] = [re.sub(r"[;=,\s]", "_", n) for n in names]
                    vout.write(rec)
            data = out_path.read_bytes()

    name = src.name[:-len(".vcf.gz")] + ".reviewed.vcf"
    return send_file(io.BytesIO(data), mimetype="text/plain", as_attachment=True, download_name=name)


@app.route("/api/depth")
@app.route("/<path:browser>/api/depth")
def site_depth(browser=None):
    """
    Quick look at a site before it is prepared: how many reads and read pairs
    (fragments) cross it, and which versions of it are already prepared.

        {"reads": 171234, "pairs": 120345,
         "prepared": {"all": false, "samples": [2000]}}
    """
    serverpath = DATA_DIR if browser is None else DATA_DIR / browser
    chrom = request.args.get("chrom")
    pos = request.args.get("pos")
    try:
        pos0 = int(pos) - 1
    except (TypeError, ValueError):
        return jsonify({"error": "pos must be a number"}), 400

    bam_dir = ensure_current_bam_dir(serverpath)
    stem = site_stem(chrom, pos)
    sample_re = re.compile(re.escape(stem) + r"\.s(\d+)\.calmd\.bam$")
    samples = sorted(int(m.group(1)) for p in bam_dir.glob(f"{stem}.s*.calmd.bam")
                     if (m := sample_re.match(p.name)))
    prepared = {"all": (bam_dir / f"{stem}.calmd.bam").exists(), "samples": samples}

    try:
        with pysam.AlignmentFile(str(serverpath / "alignment.markdup.bam"), "rb") as bam:
            names = set()
            reads = 0
            for read in bam.fetch(chrom, pos0, pos0 + 1):
                if read.is_unmapped:
                    continue
                reads += 1
                names.add(read.query_name)
    except (OSError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify({"reads": reads, "pairs": len(names), "prepared": prepared})

@app.route("/api/cache_info")
@app.route("/<path:browser>/api/cache_info")
def cache_info(browser=None):
    cache_dir = CACHE_ROOT if browser is None else CACHE_ROOT / browser
    print("CACHE_INFO HIT:", browser, cache_dir, flush=True)
    if not cache_dir.exists():
        return jsonify({"size": 0})
    
    total_size = sum(f.stat().st_size for f in cache_dir.rglob('*') if f.is_file())
    return jsonify({"size": total_size})

@app.route("/api/clear_cache", methods=["POST"])
@app.route("/<path:browser>/api/clear_cache", methods=["POST"])
def clear_cache(browser=None):
    cache_dir = CACHE_ROOT if browser is None else CACHE_ROOT / browser
    print("CACHE_CLEAR HIT:", browser, cache_dir, flush=True)
    skipped = []
    if cache_dir.exists():
        # On Windows (including a Windows folder mounted into Docker) a file
        # that is being served cannot be deleted. Skip such files rather than
        # failing; they are removed by the next clear.
        def _skip(func, path, exc):
            skipped.append(str(path))

        shutil.rmtree(cache_dir, onexc=_skip)
        cache_dir.mkdir(parents=True, exist_ok=True)
        if skipped:
            print(f"Cache clear skipped {len(skipped)} file(s) in use", flush=True)
    return jsonify({"success": True, "skipped": len(skipped)})

@app.route("/api/load_default_presets")
def load_default_presets():
    preset_path = DATA_DIR / "filters" / "filters.json"
    if not preset_path.exists():
        return jsonify({"error": "Presets file not found"}), 404
    
    try:
        with open(preset_path, "r", encoding="utf-8") as f:
            return jsonify(json.load(f))
    except Exception as e:
        return jsonify({"error": str(e)}), 500



# -----------------------------
# Static routes (original server)
# -----------------------------

# @app.route('/')
# @app.route("/<path:browser>")
# def index(browser=None):
#     if browser is None:
#         serverpath = DATA_DIR
#     else:
#         serverpath = DATA_DIR / browser

#     # If it's a directory, serve the HTML file from within it
#     if serverpath.is_dir():
#         print("INDEX HIT:", browser, serverpath, flush=True)
#         return send_from_directory(serverpath, 'interactive_variants.html')

#     # Otherwise treat it as a file request
#     print("FILE HIT:", DATA_DIR, browser, flush=True)
#     return send_from_directory(DATA_DIR, browser, mimetype=get_mimetype(browser))

@app.route("/")
@app.route("/<path:browser>")
def index(browser=None):
    # A single catch-all rule: registering "/<path:browser>/" as well makes
    # Werkzeug redirect every file URL (e.g. ".../interactive_variants.html")
    # to a trailing-slash version, which breaks the page's relative URLs.
    serverpath = DATA_DIR if browser is None else DATA_DIR / browser

    if serverpath.is_dir():
        if browser is not None and not request.path.endswith("/"):
            # Directory URLs need the trailing slash so that config.js and the
            # /api/... calls resolve relative to the instance.
            print("redirecting: ", DATA_DIR, request.path, flush=True)
            return redirect(request.path + "/", code=308)

        print("Serve interactive variants : ", DATA_DIR, serverpath, flush=True)
        return send_from_directory(serverpath, "interactive_variants.html")

    print("SERVERPATH: ", DATA_DIR, serverpath, flush=True)
    return send_from_directory(DATA_DIR, browser, mimetype=get_mimetype(browser))

def main():
    parser = argparse.ArgumentParser(description="Start Allelescope file server")
    parser.add_argument("--data-dir", "-d", default=".")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--debug", action="store_true")

    args = parser.parse_args()

    global DATA_DIR, CACHE_ROOT

    DATA_DIR = Path(args.data_dir).resolve()
    CACHE_ROOT = DATA_DIR / "cache"
    CACHE_ROOT.mkdir(exist_ok=True)

    if not DATA_DIR.is_dir():
        print(f"DATA_DIR does not exist: {DATA_DIR}", file=sys.stderr)
        sys.exit(1)

    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
