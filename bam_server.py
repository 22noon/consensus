import hashlib
import os
import subprocess
from flask import Flask, send_from_directory, request, send_file, abort, Response,request, jsonify
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
    print(f"Serving file: {path} with MIME type: {mimetype} range: {range_header}")


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


def build_samtools_view_cmd(filename, Flagf=0, FlagF=0, Tagf=0, TagF=0, soft_clip_threshold=0,soft_clip_thresholdF=0, edit_distance_threshold=0, edit_distance_thresholdF=0 ):
    tag_name="XO"
    FlagF = "0" if FlagF is None else FlagF
    Flagf = "0" if Flagf is None else Flagf
    TagF = "0" if TagF is None else TagF
    Tagf = "0" if Tagf is None else Tagf
    soft_clip_threshold = "0" if soft_clip_threshold is None else soft_clip_threshold
    soft_clip_thresholdF = "0" if soft_clip_thresholdF is None else soft_clip_thresholdF
    edit_distance_threshold = "0" if edit_distance_threshold is None else edit_distance_threshold
    edit_distance_thresholdF = "0" if edit_distance_thresholdF is None else edit_distance_thresholdF

    print(f"Building samtools view command for filename: {filename}, Flagf: {Flagf}, FlagF: {FlagF}, Tagf: {Tagf}, TagF: {TagF}, soft_clip_threshold: {soft_clip_threshold} tag_name: {tag_name}    ")
    cmd = ["samtools", "view","-b"]
    if int(Flagf):
        cmd += ["-f", Flagf]
    if int(FlagF):
        cmd += ["-F", FlagF]
    # Custom tag bitmask filtering
    conditions = []
    if int(Tagf):
        conditions.append(f"([{tag_name}] & {Tagf}) == {Tagf}")
    if int(TagF):
        conditions.append(f"([{tag_name}] & {TagF}) == 0")
    if int(soft_clip_threshold) > 0:
        conditions.append(f"([SL] + [SR]) <= {soft_clip_threshold}")
    if int(soft_clip_thresholdF) > 0:
        conditions.append(f"([SL] + [SR]) >= {soft_clip_thresholdF}")
    if int(edit_distance_threshold) > 0:
        conditions.append(f"([NM] <= {edit_distance_threshold})")
    if int(edit_distance_thresholdF) > 0:
        conditions.append(f"([NM] >= {edit_distance_thresholdF})")

    if conditions:
        expr = " && ".join(conditions)
        cmd += ["-e", expr]
    cmd.append(f"bam/{filename}.bam")
    return cmd


def build_cache_key(filename,Flagf, FlagF, Tagf, TagF, soft_clip_threshold,soft_clip_thresholdF, edit_distance_threshold, edit_distance_thresholdF):
    print(f"Building cache key for filename: {filename}, Flagf: {Flagf}, FlagF: {FlagF}, Tagf: {Tagf}, TagF: {TagF}, soft_clip_threshold: {soft_clip_threshold}, soft_clip_thresholdF: {soft_clip_thresholdF}, edit_distance_threshold: {edit_distance_threshold}, edit_distance_thresholdF: {edit_distance_thresholdF}")
    key = f"{filename}|{Flagf}|{FlagF}|{Tagf}|{TagF}|{soft_clip_threshold}|{soft_clip_thresholdF}|{edit_distance_threshold}|{edit_distance_thresholdF}".encode()
    return hashlib.md5(key).hexdigest()

def build_filtered_bam(filename,cache_key, Flagf, FlagF, Tagf, TagF, soft_clip_threshold, soft_clip_thresholdF, edit_distance_threshold, edit_distance_thresholdF):
    filtered_bam = os.path.join(CACHE_DIR, f"{cache_key}.bam")
    filtered_bai = filtered_bam + ".bai"

    if os.path.exists(filtered_bam) and os.path.exists(filtered_bai):
        return filtered_bam

    cmd_view = build_samtools_view_cmd(filename, Flagf, FlagF, Tagf, TagF,soft_clip_threshold,soft_clip_thresholdF, edit_distance_threshold, edit_distance_thresholdF)
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
    Flagf = request.args.get("Flagf")
    FlagF = request.args.get("FlagF")
    Tagf = request.args.get("Tagf")
    TagF = request.args.get("TagF")
    soft_clip_threshold = request.args.get("SoftClip")
    soft_clip_thresholdF = request.args.get("SoftClipF")
    edit_distance_threshold = request.args.get("EditDistance")
    edit_distance_thresholdF = request.args.get("EditDistanceF")
    print(f"Received request for BAM: {filename} with parameters Flagf: {Flagf}, FlagF: {FlagF}, Tagf: {Tagf}, TagF: {TagF}, soft_clip_threshold: {soft_clip_threshold}, soft_clip_thresholdF: {soft_clip_thresholdF}, edit_distance_threshold: {edit_distance_threshold}, edit_distance_thresholdF: {edit_distance_thresholdF}")

    cache_key = build_cache_key(filename,Flagf, FlagF, Tagf, TagF, soft_clip_threshold, soft_clip_thresholdF, edit_distance_threshold, edit_distance_thresholdF)
    print(f"BAM: Checking cache for BAM: {filename} with key: {cache_key}")
    bam_path = build_filtered_bam(filename, cache_key,  Flagf, FlagF, Tagf, TagF, soft_clip_threshold, soft_clip_thresholdF, edit_distance_threshold, edit_distance_thresholdF)
    print(f"Serving BAM from path: {bam_path}")

    return send_file_with_range(bam_path)

@app.route("/bai/<filename>")
def serve_bai(filename):
    rg = request.args.get("rg")
    Flagf = request.args.get("Flagf")
    FlagF = request.args.get("FlagF")
    Tagf = request.args.get("Tagf")
    TagF = request.args.get("TagF")
    soft_clip_threshold = request.args.get("SoftClip")
    soft_clip_thresholdF = request.args.get("SoftClipF")
    edit_distance_threshold = request.args.get("EditDistance")
    edit_distance_thresholdF = request.args.get("EditDistanceF")
    print(f"Received request for BAi: {filename} with parameters Flagf: {Flagf}, FlagF: {FlagF}, Tagf: {Tagf}, TagF: {TagF}, soft_clip_threshold: {soft_clip_threshold}, soft_clip_thresholdF: {soft_clip_thresholdF}, edit_distance_threshold: {edit_distance_threshold}, edit_distance_thresholdF: {edit_distance_thresholdF}")

# Compute cache key
    key = build_cache_key(filename,Flagf, FlagF, Tagf, TagF, soft_clip_threshold, soft_clip_thresholdF, edit_distance_threshold, edit_distance_thresholdF)
    bam_path = CACHE_DIR / f"{key}.bam"
    bai_path = CACHE_DIR / f"{key}.bam.bai"

    # Ensure BAM and BAI exist — generate on demand!
    print(f"BAI: Checking cache for BAM: {bam_path} and BAI: {bai_path}")
    if not bam_path.exists() or not bai_path.exists():
        build_filtered_bam(filename, key,  Flagf, FlagF, Tagf, TagF, soft_clip_threshold, soft_clip_thresholdF, edit_distance_threshold, edit_distance_thresholdF)
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

@app.route("/extract")
def extract_bam():
    Chrom = request.args.get("chrom")
    Pos = request.args.get("pos")
    print(f"Received request for BAM: with parameters Chrom: {Chrom}, Pos: {Pos}")
    if subprocess.run(["bash", "extract_bam.sh", Chrom, Pos], check=False):
        return jsonify({"success": True,"message": "Task completed",}), 200
    else:
        return jsonify({"success": False,"message": "Task failed",}), 500


if __name__ == "__main__":
    app.run(debug=True)
