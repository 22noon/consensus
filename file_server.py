#!/usr/bin/env python3
from flask import Flask, send_from_directory, Response
import mimetypes
import os
import argparse
import sys

app = Flask(__name__)

# Set the directory where your files are located
# Default to the directory of this file
default_dir = os.path.dirname(os.path.abspath(__file__))

# Parse command-line argument (use parse_known_args so Flask/other tooling args are ignored)
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-dir', '-d', default=default_dir, help='Path to the data directory to serve')
args, _ = parser.parse_known_args()

DATA_DIR = os.path.abspath(args.data_dir)
if not os.path.isdir(DATA_DIR):
    print(f"DATA_DIR does not exist: {DATA_DIR}", file=sys.stderr)
    sys.exit(1)

@app.route('/')
def index():
    return send_from_directory(DATA_DIR, 'interactive_viewer.html')

@app.route('/<path:filename>')
def serve_file(filename):
    # Ensure proper MIME types
    mimetype, _ = mimetypes.guess_type(filename)

    # Set specific MIME types for genomics files
    if filename.endswith('.bam'):
        mimetype = 'application/octet-stream'
    elif filename.endswith('.bai'):
        mimetype = 'application/octet-stream'
    elif filename.endswith('.fastq') or filename.endswith('.fq'):
        mimetype = 'text/plain'
    elif filename.endswith('.html'):
        mimetype = 'text/html'

    return send_from_directory(DATA_DIR, filename, mimetype=mimetype)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, debug=True)
