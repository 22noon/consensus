import hashlib
import os
import subprocess
from flask import Flask, send_from_directory, request, send_file, abort
import mimetypes
from pathlib import Path


app = Flask(__name__)

CACHE_DIR = "cache"
os.makedirs(CACHE_DIR, exist_ok=True)
CACHE_DIR = Path("cache")


def build_cache_key(filename,rg, Flagf, FlagF, custom_tag):
    key = f"{filename}|{rg}|{Flagf}|{FlagF}|{custom_tag}".encode()
    return hashlib.md5(key).hexdigest()

def build_filtered_bam(filename,cache_key, rg, Flagf, FlagF, custom_tag):
    filtered_bam = os.path.join(CACHE_DIR, f"{cache_key}.bam")
    filtered_bai = filtered_bam + ".bai"

    if os.path.exists(filtered_bam) and os.path.exists(filtered_bai):
        return filtered_bam

    tag_filters = []
    Extra_Params = ["-f",Flagf,"-F",FlagF]
    if rg:
        Extra_Params.append(f"-r {rg}")

    if custom_tag:
        Extra_Params.append(custom_tag)  # e.g. "XX:i:1"


    cmd_view = ["samtools", "view", "-b", f"bam/{filename}.bam"] #+ tag_filters
    if rg:
        cmd_view[2:2] = Extra_Params # insert before INPUT_BAM
    with open(filtered_bam, "wb") as out:
        subprocess.run(cmd_view, check=True, stdout=out)
        print(f"Command: {' '.join(cmd_view)} > {filtered_bam} ")

    subprocess.run(["samtools", "index", filtered_bam], check=True)

    return filtered_bam

@app.route("/")
def root():
    return app.send_static_file("index.html")

@app.route("/bam/<filename>")
def serve_bam(filename):
    rg = request.args.get("rg")
    Flagf = request.args.get("ff")
    FlagF = request.args.get("FF")
    custom_tag = request.args.get("tag")

    cache_key = build_cache_key(filename,rg, Flagf, FlagF, custom_tag)
    bam_path = build_filtered_bam(filename, cache_key, rg, Flagf, FlagF, custom_tag)

    return send_file(bam_path, mimetype="application/octet-stream")

@app.route("/bai/<filename>")
def serve_bai(filename):
    rg = request.args.get("rg")
    Flagf = request.args.get("ff")
    FlagF = request.args.get("FF")
    custom_tag = request.args.get("tag")

# Compute cache key
    key = build_cache_key(filename,rg, Flagf, FlagF, custom_tag)
    bam_path = CACHE_DIR / f"{key}.bam"
    bai_path = CACHE_DIR / f"{key}.bam.bai"

    # Ensure BAM and BAI exist — generate on demand!
    if not bam_path.exists() or not bai_path.exists():
        build_filtered_bam(filename, key, rg, Flagf, FlagF, custom_tag)

    return send_file(bai_path, mimetype="application/octet-stream")

@app.route('/<path:filename>')
def serve_file(filename):
    # Serve a file by name/path relative to the current working directory
    requested = Path(filename)  # convert filename (str) to a Path

    if not requested.is_absolute():
        requested = Path.cwd() / requested

    print(f"Requested file: {requested}")
    # Optional: prevent directory traversal outside cwd

    if not requested.exists() or not requested.is_file():
        return abort(404)

    # Guess MIME type and override for genomics/text files
    mimetype, _ = mimetypes.guess_type(str(requested))
    suffix = requested.suffix.lower()
    if suffix in ('.bam', '.bai', '.fai'):
        mimetype = 'application/octet-stream'
    elif suffix in ('.fastq', '.fq'):
        mimetype = 'text/plain'
    elif suffix == '.html':
        mimetype = 'text/html'

    return send_file(str(requested), mimetype=mimetype or 'application/octet-stream')

if __name__ == "__main__":
    app.run(debug=True)
