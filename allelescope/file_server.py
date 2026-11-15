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
import json
from pathlib import Path

from .process_variant import process_variant, ensure_current_bam_dir, PROCESSING_VERSION
from .read_filters import filter_reads
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

def processed_bam(serverpath, Chrom, Pos, Ref):
    return process_variant(serverpath, Chrom, Pos, Ref)


def get_filtered_reads_bam(serverpath: Path, filename: str, **params):
    """
    Return (header, reads) for the per-site BAM after applying the filter
    parameters. The filter rules live in read_filters.py.
    """
    token = params.get("Token") or ""
    bam_file = ensure_current_bam_dir(serverpath) / f"{filename}{token}.calmd.bam"
    print(f"Getting filtered reads from BAM for {filename} with params {params} from {bam_file}", flush=True)
    if not bam_file.exists():
        processed_bam(serverpath, params.get("Chrom"), params.get("Pos"), params.get("Ref"))

    return filter_reads(bam_file, params)

def build_filtered_bam(serverpath, cache_dir, filename, cache_key, params):
    bam_path = cache_dir / f"{cache_key}.bam"
    bai_path = cache_dir / f"{cache_key}.bam.bai"

    if bam_path.exists() and bai_path.exists():
        return bam_path

    header, reads = get_filtered_reads_bam(serverpath, filename, **params)
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.index(str(bam_path))
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
    filename = f"variant_{Chrom}_{Pos}"
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
    filename = f"variant_{Chrom}_{Pos}"

    bam_file = ensure_current_bam_dir(serverpath) / f"{filename}.calmd.bam"
    if not bam_file.exists():
        processed_bam(serverpath, Chrom, Pos, Ref)

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

    out_bam = serverpath / f"bam/variant_{chrom}_{pos}{token}.calmd.bam"
    in_bam = serverpath / f"bam/variant_{chrom}_{pos}.calmd.bam"
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

    out_bam = serverpath / f"bam/variant_{chrom}_{pos}{token}.calmd.bam"
    in_bam = serverpath / f"bam/variant_{chrom}_{pos}.calmd.bam"
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
    filename  = f"variant_{chrom}_{pos}"
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

    mock_alleles = []

    result = processed_bam(serverpath, Chrom, Pos, Ref)
    if result == 0:
        stat_file=f"bam/variant_{Chrom}_{Pos}.bam.indelstats.txt"
        with open(serverpath / stat_file) as f:
            for line in f:
                line = line.strip()
                fields = line.split("\t")
                xa, count, pct, allele_type = fields[:4]
                count = int(count)
                pct = round(float(pct)*100, 2)

                mock_alleles.append({
                    "xa": xa,
                    "xd": 0,
                    "xn": 0,
                    "count": count,
                    "pct": pct,
                    "type": allele_type,
                })

        return jsonify({"success": True, "alleles": mock_alleles})
    return jsonify({"success": False}), 500

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
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
    return jsonify({"success": True})

@app.route("/api/load_default_presets")
def load_default_presets():
    preset_path = DATA_DIR / "filters" / "filters.json"
    if not preset_path.exists():
        return jsonify({"error": "Presets file not found"}), 404
    
    try:
        with open(preset_path, "r") as f:
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
