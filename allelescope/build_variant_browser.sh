#!/bin/bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# build_variant_browser.sh
# Prepares alignment and reference files for the variant browser.
#
# Usage:
#   ./build_variant_browser.sh -b BAM -r REF [-v VARIANTS] [-g GENOME_NAME] [-o OUT_DIR] [-t THREADS]
#
# Required:
#   -b  Input BAM file
#   -r  Reference FASTA file
#
# Optional:
#   -v  Variant file (VCF/VCF.gz)
#   -g  Genome name (default: first contig in VCF, or first reference in BAM)
#   -o  Output directory (default: ./out_dir)
#   -t  Threads (default: 10)
# ─────────────────────────────────────────────────────────────────────────────

# ── Defaults ─────────────────────────────────────────────────────────────────
THREADS=10
OUT_DIR="./out_dir"
GENOME_NAME=""
VARIANTS=""

# ── Colours ──────────────────────────────────────────────────────────────────
RED=$'\033[0;31m'; GREEN=$'\033[0;32m'; YELLOW=$'\033[1;33m'
CYAN=$'\033[0;36m'; BOLD=$'\033[1m'; RESET=$'\033[0m'

info()    { printf "${CYAN}[INFO]${RESET}  %s\n" "$*"; }
success() { printf "${GREEN}[OK]${RESET}    %s\n" "$*"; }
warn()    { printf "${YELLOW}[WARN]${RESET}  %s\n" "$*"; }
error()   { printf "${RED}[ERROR]${RESET} %s\n" "$*" >&2; exit 1; }
step()    { printf "\n${BOLD}▶ %s${RESET}\n" "$*"; }

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF

${BOLD}Usage:${RESET}
  $(basename "$0") -b BAM -r FASTA [-v VCF] [-g GENOME_NAME] [-o OUT_DIR] [-t THREADS]

${BOLD}Required:${RESET}
  -b  Input BAM file
  -r  Reference FASTA file

${BOLD}Optional:${RESET}
  -v  Variant file (VCF or VCF.gz, must be sorted)
  -g  Genome name  (default: first contig in VCF → first reference in BAM)
  -o  Output directory  (default: ./out_dir)
  -t  Threads  (default: 10)
  -h  Show this help

EOF
    exit 1
}

# ── Parse arguments ───────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

while getopts ":b:r:v:g:o:t:h" opt; do
    case $opt in
        b) BAM="$OPTARG"         ;;
        r) REF="$OPTARG"         ;;
        v) VARIANTS="$OPTARG"    ;;
        g) GENOME_NAME="$OPTARG" ;;
        o) OUT_DIR="$OPTARG"     ;;
        t) THREADS="$OPTARG"     ;;
        h) usage                 ;;
        :) error "Option -$OPTARG requires an argument." ;;
        \?) error "Unknown option: -$OPTARG" ;;
    esac
done

# ── Validate required inputs ──────────────────────────────────────────────────
step "Validating inputs"

[[ -z "${BAM:-}" ]] && error "BAM file is required (-b)"
[[ -z "${REF:-}" ]] && error "Reference FASTA is required (-r)"

[[ -f "$BAM" ]] || error "BAM file not found: $BAM"
[[ -f "$REF" ]] || error "FASTA file not found: $REF"
[[ -n "$VARIANTS" && ! -f "$VARIANTS" ]] && error "Variant file not found: $VARIANTS"

NEED_TABIX=false
[[ -n "$VARIANTS" ]] && NEED_TABIX=true

for cmd in samtools; do
    command -v "$cmd" &>/dev/null || error "'$cmd' not found in PATH"
done
if $NEED_TABIX; then
    for cmd in tabix bgzip; do
        command -v "$cmd" &>/dev/null || error "'$cmd' not found in PATH"
    done
fi

success "All inputs present"

# ── Resolve genome name ───────────────────────────────────────────────────────
step "Resolving genome name"

if [[ -z "$GENOME_NAME" ]]; then
    # Try first CHROM in the VCF (skip header lines)
    if [[ -n "$VARIANTS" ]]; then
       GENOME_NAME=$(zcat "$VARIANTS" | grep -v "^#" | awk 'NR==1{print $1}' || true)
    fi

    # Fallback: first reference sequence name in BAM header
    if [[ -z "$GENOME_NAME" ]]; then
        [[ -n "$VARIANTS" ]] && warn "No variants found in VCF; falling back to first reference in BAM"
        GENOME_NAME=$(samtools view -H "$BAM" | awk '/^@SQ/ { split($2, a, ":"); print a[2]; exit }')
    fi

    [[ -z "$GENOME_NAME" ]] && error "Could not determine genome name automatically. Please supply -g."
    info "Genome name inferred: ${BOLD}$GENOME_NAME${RESET}"
else
    info "Genome name provided: ${BOLD}$GENOME_NAME${RESET}"
fi

# ── Prepare output directory ──────────────────────────────────────────────────
step "Preparing output directory: $OUT_DIR"
mkdir -p "$OUT_DIR"

MARKDUP_BAM="$OUT_DIR/alignment.markdup.bam"
TMP_PREFIX="$OUT_DIR/tmp_markdup"

# ── Process BAM ───────────────────────────────────────────────────────────────
step "Processing BAM (collate → fixmate → sort → markdup)"
info "Threads: $THREADS"

samtools collate -u -@ "$THREADS" -O "$BAM" "$TMP_PREFIX" \
    | samtools fixmate -@ "$THREADS" -m -u - - \
    | samtools sort   -@ "$THREADS" -u -O bam - \
    | samtools markdup -@ "$THREADS" - "$MARKDUP_BAM"

success "markdup BAM written: $MARKDUP_BAM"

step "Indexing BAM"
samtools index "$MARKDUP_BAM"
success "BAM index written: ${MARKDUP_BAM}.bai"

# ── Reference FASTA ───────────────────────────────────────────────────────────
step "Copying and indexing reference FASTA"

REF_DEST="$OUT_DIR/$(basename "$REF")"
cp "$REF" "$REF_DEST"
samtools faidx "$REF_DEST"
success "FASTA index written: ${REF_DEST}.fai"

# ── Variant file (optional) ───────────────────────────────────────────────────
VCF_DEST=""
if [[ -n "$VARIANTS" ]]; then
    step "Copying and indexing variant file"

    VCF_DEST="$OUT_DIR/$(basename "$VARIANTS")"

    # Ensure the VCF is bgzipped before tabix indexing
    if [[ "$VARIANTS" == *.gz ]]; then
        cp "$VARIANTS" "$VCF_DEST"
    else
        warn "VCF is not gzipped — bgzipping to ${VCF_DEST}.gz"
        VCF_DEST="${VCF_DEST}.gz"
        bgzip -c "$VARIANTS" > "$VCF_DEST"
    fi

    tabix -p vcf "$VCF_DEST"
    success "Variant index written: ${VCF_DEST}.tbi"
else
    info "No variant file provided — skipping VCF step"
fi

# ── Run viewer generator ──────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENERATOR="$SCRIPT_DIR/generate_viewer.py"

# Build VCF arg only if we have a file
VCF_ARG=""
[[ -n "$VCF_DEST" ]] && VCF_ARG="--vcf-path \"$VCF_DEST\""

if [[ -f "$GENERATOR" ]]; then
    step "Running generate_viewer.py"
    python "$GENERATOR" \
        --genome-name     "$GENOME_NAME" \
        --reference-fasta "$REF_DEST" \
        --bam-path        "$MARKDUP_BAM" \
        ${VCF_DEST:+--vcf-path "$VCF_DEST"} \
        --output-dir      "$OUT_DIR"
    success "Viewer generated in: $OUT_DIR"
else
    warn "generate_viewer.py not found at $GENERATOR — skipping viewer generation."
    warn "Run manually with:"
    echo
    echo "  python generate_viewer.py \\"
    echo "      --genome-name    \"$GENOME_NAME\" \\"
    echo "      --reference-fasta \"$REF_DEST\" \\"
    echo "      --bam-path       \"$MARKDUP_BAM\" \\"
    [[ -n "$VCF_DEST" ]] && echo "      --vcf-path       \"$VCF_DEST\" \\"
    echo "      --output-dir     \"$OUT_DIR\""
    echo
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo -e "\n${BOLD}────────────────────────── Summary ──────────────────────────────${RESET}"
echo -e "  Genome name   : ${CYAN}$GENOME_NAME${RESET}"
echo -e "  Output dir    : ${CYAN}$OUT_DIR${RESET}"
echo -e "  BAM           : $(basename "$MARKDUP_BAM")"
echo -e "  FASTA         : $(basename "$REF_DEST")"
if [[ -n "$VCF_DEST" ]]; then
    echo -e "  Variants      : $(basename "$VCF_DEST")"
else
    echo -e "  Variants      : (none)"
fi
echo -e "${BOLD}─────────────────────────────────────────────────────────────────${RESET}\n"
