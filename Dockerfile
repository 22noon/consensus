# Allelescope in a container: Python 3.12, pysam and samtools/htslib 1.19.2
# on Ubuntu 24.04. Runs on Linux, macOS and Windows (Docker Desktop / WSL2).
#
#   docker build -t allelescope .
#   docker run --rm -p 127.0.0.1:8000:8000 allelescope demo
#
# See docker/README.md for building and serving your own data.

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH

# samtools and tabix/bgzip are needed by build_variant_browser.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates \
        python3 \
        python3-venv \
        samtools \
        tabix \
 && update-ca-certificates \
 && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv "$VIRTUAL_ENV"

# Install the package (only what the build needs, so code edits elsewhere
# in the repository do not invalidate this layer).
WORKDIR /opt/allelescope/src
COPY pyproject.toml poetry.lock ./
COPY allelescope ./allelescope
RUN pip install . \
 && python -c "import pysam, allelescope; print('pysam', pysam.__version__)"

# Demo dataset, used by the 'demo' command.
COPY demo /opt/allelescope/demo
COPY docker/entrypoint.sh /usr/local/bin/allelescope-entrypoint
RUN chmod 0755 /usr/local/bin/allelescope-entrypoint

# Run as an unprivileged user. On Linux hosts add
# --user "$(id -u):$(id -g)" so files written to /data belong to you.
RUN useradd --create-home --uid 1000 allelescope \
 && mkdir -p /data \
 && chown allelescope:allelescope /data
USER allelescope

# /data holds the browser instances (default: /data/Browsers) and any input
# files you want to build from. Mount a host folder here.
WORKDIR /data
VOLUME ["/data"]
EXPOSE 8000

ENTRYPOINT ["allelescope-entrypoint"]
CMD ["serve"]
