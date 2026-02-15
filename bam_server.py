import hashlib
import os
import subprocess
from flask import Flask, send_from_directory, request, send_file, abort, Response,request
import mimetypes
from pathlib import Path


app = Flask(__name__)

CACHE_DIR = "cache"
os.makedirs(CACHE_DIR, exist_ok=True)
CACHE_DIR = Path("cache")

def send_file_with_range(path: Path, mimetype="application/octet-stream"):
    path = Path(path)
    if not path.exists() or not path.is_file():
        return Response(status=404)

    file_size = path.stat().st_size
    range_header = request.headers.get("Range")

    # Always advertise range support
    def base_headers(length=None):
        h = {
            "Accept-Ranges": "bytes",
            "Content-Type": mimetype,
            "Cache-Control": "no-cache",
        }
        if length is not None:
            h["Content-Length"] = str(length)
        return h

    if not range_header:
        with open(path, "rb") as f:
            data = f.read()
        return Response(data, status=200, headers=base_headers(len(data)))

    # Parse Range
    try:
        unit, rng = range_header.split("=")
        if unit.strip().lower() != "bytes":
            raise ValueError
        start_str, end_str = rng.split("-")
        if start_str == "":
            # suffix range: bytes=-N
            length = int(end_str)
            if length <= 0:
                raise ValueError
            start = max(0, file_size - length)
            end = file_size - 1
        else:
            start = int(start_str)
            end = int(end_str) if end_str else file_size - 1
    except Exception:
        # Malformed or unsupported → just return full file
        with open(path, "rb") as f:
            data = f.read()
        return Response(data, status=200, headers=base_headers(len(data)))

    # If out of range → return full file (workaround for tiny empty BAMs)
    if start >= file_size or start < 0 or end < start:
        with open(path, "rb") as f:
            data = f.read()
        return Response(data, status=200, headers=base_headers(len(data)))

    end = min(end, file_size - 1)
    length = end - start + 1

    with open(path, "rb") as f:
        f.seek(start)
        data = f.read(length)

    headers = base_headers(length)
    headers["Content-Range"] = f"bytes {start}-{end}/{file_size}"
    return Response(data, status=206, headers=headers)


def build_samtools_view_cmd(filename, Flagf=0, FlagF=0, Tagf=0, TagF=0, tag_name="XO"):
    cmd = ["samtools", "view","-b","-h"]
    if Flagf:
        cmd += ["-f", str(Flagf)]
    if FlagF:
        cmd += ["-F", str(FlagF)]
    # Custom tag bitmask filtering
    conditions = []
    if Tagf:
        conditions.append(f"([{tag_name}] & {Tagf}) == {Tagf}")
    if TagF:
        conditions.append(f"([{tag_name}] & {TagF}) == 0")

    if conditions:
        expr = " && ".join(conditions)
        cmd += ["-e", expr]
    cmd.append(f"bam/{filename}.bam")
    return cmd



def build_cache_key(filename,rg, Flagf, FlagF, Tagf, TagF):
    print(f"Building cache key for filename: {filename}, rg: {rg}, Flagf: {Flagf}, FlagF: {FlagF}, Tagf: {Tagf}, TagF: {TagF}")
    key = f"{filename}|{rg}|{Flagf}|{FlagF}|{Tagf}|{TagF}".encode()
    return hashlib.md5(key).hexdigest()

def build_filtered_bam(filename,cache_key, rg, Flagf, FlagF, Tagf, TagF):
    filtered_bam = os.path.join(CACHE_DIR, f"{cache_key}.bam")
    filtered_bai = filtered_bam + ".bai"

    if os.path.exists(filtered_bam) and os.path.exists(filtered_bai):
        return filtered_bam

    cmd_view = build_samtools_view_cmd(filename, Flagf, FlagF, Tagf, TagF)
    print(f"Command: {' '.join(cmd_view)} > {filtered_bam} ")
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
    Flagf = request.args.get("Flagf")
    FlagF = request.args.get("FlagF")
    Tagf = request.args.get("Tagf")
    TagF = request.args.get("TagF")

    cache_key = build_cache_key(filename,rg, Flagf, FlagF, Tagf, TagF)
    print(f"BAM: Checking cache for BAM: {filename} with key: {cache_key}")
    bam_path = build_filtered_bam(filename, cache_key, rg, Flagf, FlagF, Tagf, TagF)

    return send_file_with_range(bam_path)

@app.route("/bai/<filename>")
def serve_bai(filename):
    rg = request.args.get("rg")
    Flagf = request.args.get("Flagf")
    FlagF = request.args.get("FlagF")
    Tagf = request.args.get("Tagf")
    TagF = request.args.get("TagF")

# Compute cache key
    key = build_cache_key(filename,rg, Flagf, FlagF, Tagf, TagF)
    bam_path = CACHE_DIR / f"{key}.bam"
    bai_path = CACHE_DIR / f"{key}.bam.bai"

    # Ensure BAM and BAI exist — generate on demand!
    print(f"BAI: Checking cache for BAM: {bam_path} and BAI: {bai_path}")
    if not bam_path.exists() or not bai_path.exists():
        build_filtered_bam(filename, key, rg, Flagf, FlagF, Tagf, TagF)

    return send_file_with_range(bai_path)

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
