#!/usr/bin/env python3
from flask import Flask, send_from_directory, Response, request, jsonify
import mimetypes
import os
import argparse
import sys
import hashlib
import subprocess
from pathlib import Path

app = Flask(__name__)

default_dir = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-dir', '-d', default=default_dir)
args, _ = parser.parse_known_args()
DATA_DIR = Path(os.path.abspath(args.data_dir))
CACHE_ROOT = DATA_DIR / "cache"
CACHE_ROOT.mkdir(exist_ok=True)

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
    key = "|".join(str(a) for a in args).encode()
    return hashlib.md5(key).hexdigest()


def processed_bam(serverpath, Chrom, Pos, Ref):
    print(f"Processing BAM for {Chrom}:{Pos} ref={Ref}", flush=True)
    return subprocess.run(
        ["bash", "extract_bam.sh", Chrom, Pos, Ref, str(serverpath)],
        check=False
    )


def build_samtools_view_cmd(serverpath, filename, **params):
    bam_file = serverpath / "bam" / f"{filename}.bam"
    if not bam_file.exists():
        processed_bam(serverpath, params.get("Chrom"), params.get("Pos"), params.get("Ref"))

    cmd = ["samtools", "view", "-b"]

    if int(params.get("Flagf", 0)):
        cmd += ["-f", params["Flagf"]]
    if int(params.get("FlagF", 0)):
        cmd += ["-F", params["FlagF"]]

    conditions = []
    if int(params.get("Tagf", 0)):
        conditions.append(f"([XO] & {params['Tagf']}) == {params['Tagf']}")
    if int(params.get("TagF", 0)):
        conditions.append(f"([XO] & {params['TagF']}) == 0")
    if int(params.get("EditDistance", 0)):
        conditions.append(f"([NM] <= {params['EditDistance']})")

    xa_filter = params.get("XAFilter", "").strip()
    if xa_filter:
        # Reads with no XA tag should pass through — samtools -e will
        # error on tag absence, so we guard with a type check first.
        # typeof() returns the SAM type char: 'Z' for string tags.
        conditions.append(f'[XA] =~ "{xa_filter}"')
    
    if conditions:
        cmd += ["-e", " && ".join(conditions)]

    cmd.append(str(bam_file))
    return cmd


def build_filtered_bam(serverpath, cache_dir, filename, cache_key, params):
    bam_path = cache_dir / f"{cache_key}.bam"
    bai_path = cache_dir / f"{cache_key}.bam.bai"

    if bam_path.exists() and bai_path.exists():
        return bam_path

    cmd = build_samtools_view_cmd(serverpath, filename, **params)
    print ("Running samtools view with command:", " ".join(cmd), flush=True)
    with open(bam_path, "wb") as out:
        subprocess.run(cmd, check=True, stdout=out)
    subprocess.run(["samtools", "index", bam_path], check=True)
    return bam_path


# -----------------------------
# BAM/BAI routes (from complex server)
# -----------------------------

@app.route("/<path:browser>/api/bai")
def serve_bai(browser):
    serverpath = DATA_DIR / browser
    cache_dir = CACHE_ROOT / browser
    cache_dir.mkdir(parents=True, exist_ok=True)
    print("SERVE_BAI HIT:", browser,serverpath,cache_dir, flush=True)

    Chrom = request.args.get("Chrom")
    Pos = request.args.get("Pos")
    filename = f"variant_{Chrom}_{Pos}"
    params = request.args.to_dict()
    cache_key = build_cache_key(filename, *params.values())

    bam_path = cache_dir / f"{cache_key}.bam"
    bai_path = cache_dir / f"{cache_key}.bam.bai"

    if not bam_path.exists() or not bai_path.exists():
        build_filtered_bam(serverpath, cache_dir, filename, cache_key, params)

    return send_file_with_range(bai_path)


@app.route("/<path:browser>/api/bam")
def serve_bam(browser):
    serverpath = DATA_DIR / browser
    cache_dir = CACHE_ROOT / browser

    print("SERVE_BAM HIT:", browser,serverpath,cache_dir, flush=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    Chrom = request.args.get("Chrom")
    Pos = request.args.get("Pos")
    Ref = request.args.get("Ref")
    filename = f"variant_{Chrom}_{Pos}"

    bam_file = serverpath / "bam" / f"{filename}.bam"
    if not bam_file.exists():
        processed_bam(serverpath, Chrom, Pos, Ref)

    params = request.args.to_dict()
    cache_key = build_cache_key(filename, *params.values())
    bam_path = build_filtered_bam(serverpath, cache_dir, filename, cache_key, params)

    return send_file_with_range(bam_path)


@app.route("/<path:browser>/api/extract")
def extract_bam(browser):
    serverpath = DATA_DIR / browser
    print("EXTRACT HIT:", browser,serverpath, flush=True)
    Chrom = request.args.get("chrom")
    Pos = request.args.get("pos")
    Ref = request.args.get("ref")

    mock_alleles = []

    result = processed_bam(serverpath, Chrom, Pos, Ref)
    if result.returncode == 0:
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


# -----------------------------
# Static routes (original server)
# -----------------------------

@app.route('/')
def index():
    return send_from_directory(DATA_DIR, 'interactive_viewer.html')


@app.route('/<path:filename>')
def serve_file(filename):
    return send_from_directory(DATA_DIR, filename, mimetype=get_mimetype(filename))


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, debug=True)#