#!/bin/sh
# Entry point for the Allelescope container.
#
#   serve                     serve every instance under $ALLELESCOPE_DATA_DIR (default)
#   demo                      build the bundled demo instance if needed, then serve
#   build <args...>           run build_variant_browser with <args>
#   <any other command>       run it as given (e.g. bash, samtools, generate_viewer)
#
# Environment:
#   ALLELESCOPE_DATA_DIR  folder holding the instances   (default /data/Browsers)
#   ALLELESCOPE_PORT      port the server listens on      (default 8000)
#   ALLELESCOPE_THREADS   samtools threads for 'demo'     (default 2)
set -eu

DATA_DIR="${ALLELESCOPE_DATA_DIR:-/data/Browsers}"
PORT="${ALLELESCOPE_PORT:-8000}"

serve() {
    mkdir -p "$DATA_DIR"
    echo "Serving Allelescope instances from $DATA_DIR"
    echo "Open http://localhost:<published port>/<instance>/ in your browser."
    # Inside the container the server must listen on all interfaces so the
    # published port reaches it; limit access with -p 127.0.0.1:PORT:PORT.
    exec start_server --data-dir "$DATA_DIR" --host 0.0.0.0 --port "$PORT"
}

cmd="${1:-serve}"
[ "$#" -gt 0 ] && shift

case "$cmd" in
    serve)
        serve
        ;;
    demo)
        if [ ! -f "$DATA_DIR/demo/interactive_variants.html" ]; then
            echo "Building the demo instance in $DATA_DIR/demo ..."
            mkdir -p "$DATA_DIR"
            build_variant_browser \
                -b /opt/allelescope/demo/demo_reads.bam \
                -r /opt/allelescope/demo/demo_ref.fa \
                -v /opt/allelescope/demo/demo_calls.vcf.gz \
                -o "$DATA_DIR/demo" \
                -t "${ALLELESCOPE_THREADS:-2}"
        fi
        echo "Demo ready: http://localhost:<published port>/demo/"
        serve
        ;;
    build)
        exec build_variant_browser "$@"
        ;;
    *)
        exec "$cmd" "$@"
        ;;
esac
