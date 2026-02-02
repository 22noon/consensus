#!/usr/bin/env python3
from flask import Flask, send_from_directory, Response
import mimetypes
import os

app = Flask(__name__)

# Set the directory where your files are located
DATA_DIR = '/home/chandana.tennakoon@pirbright.ac.uk/git/consensus/variant_browser'  # Change this to your actual directory

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
