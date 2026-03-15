#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------
# PacBio HiFi Variant Calling Script (bcftools)
# --------------------------------------------------------

LOG_DIR="${LOG_DIR:-"./pipeline_logs"}"
mkdir -p "$LOG_DIR"

exec > >(tee -a "$LOG_DIR/bcftools_general.log") \
     2> >(tee -a "$LOG_DIR/bcftools_error.log" >&2)

usage() {
    echo
    echo "Usage:"
    echo "  $0 -r reference.fa -b input.bam -o output_prefix [-t threads]"
    echo
    echo "Options:"
    echo "  -r  Reference FASTA"
    echo "  -b  Input BAM file (sorted and indexed)"
    echo "  -o  Output prefix"
    echo "  -t  Threads (default: 4)"
    echo
    exit 1
}

THREADS=4

while getopts ":r:b:o:t:" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        b) BAM="$OPTARG" ;;
        o) OUTPREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG"; usage ;;
    esac
done

if [[ -z "${REF:-}" || -z "${BAM:-}" || -z "${OUTPREFIX:-}" ]]; then
    usage
fi

# --------------------------------------------------------
# Tool checks
# --------------------------------------------------------

if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools not found in PATH"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found in PATH"
    exit 1
fi

# --------------------------------------------------------
# Input validation
# --------------------------------------------------------

if [[ ! -f "$BAM" ]]; then
    echo "Error: BAM file not found: $BAM"
    exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
    echo "Error: BAM index not found: ${BAM}.bai"
    echo "Run: samtools index $BAM"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "Error: Reference FASTA not found: $REF"
    exit 1
fi

# Ensure reference index exists
if [[ ! -f "${REF}.fai" ]]; then
    echo "Reference index not found. Creating with samtools faidx..."
    samtools faidx "$REF"
fi

# --------------------------------------------------------
# File naming
# --------------------------------------------------------

SAMPLE_NAME=$(basename "$BAM" | sed 's/\.bam//')

RAW_VCF="${OUTPREFIX}_raw.vcf.gz"
FILTERED_VCF="${OUTPREFIX}_filtered.vcf.gz"
FINAL_VCF="${OUTPREFIX}_norm.vcf.gz"

echo
echo "--------------------------------------------"
echo "PacBio HiFi Variant Calling (bcftools)"
echo "Reference : $REF"
echo "BAM       : $BAM"
echo "Sample    : $SAMPLE_NAME"
echo "Threads   : $THREADS"
echo "Output    : $FINAL_VCF"
echo "--------------------------------------------"
echo

# --------------------------------------------------------
# Skip if output exists
# --------------------------------------------------------

if [[ -f "$FINAL_VCF" ]]; then
    echo "Output VCF already exists. Skipping variant calling."
    exit 0
fi

# --------------------------------------------------------
# Variant calling
# --------------------------------------------------------

echo "Running bcftools mpileup + call..."

set +e

bcftools mpileup \
    --fasta-ref "$REF" \
    --output-type u \
    --max-depth 200 \
    --min-MQ 20 \
    --min-BQ 20 \
    --annotate FORMAT/AD,FORMAT/DP \
    --threads "$THREADS" \
    "$BAM" \
| bcftools call \
    --multiallelic-caller \
    --variants-only \
    --output-type z \
    --threads "$THREADS" \
    -o "$RAW_VCF"

PIPE_STATUS=(${PIPESTATUS[@]})

set -e

if [[ ${PIPE_STATUS[0]} -ne 0 || ${PIPE_STATUS[1]} -ne 0 ]]; then
    echo "Error: bcftools mpileup/call failed"
    exit 1
fi

# --------------------------------------------------------
# Index raw VCF
# --------------------------------------------------------

echo
echo "Indexing raw VCF..."

bcftools index "$RAW_VCF"

# --------------------------------------------------------
# Variant filtering
# --------------------------------------------------------

echo
echo "Filtering variants..."

bcftools filter \
    --include 'QUAL >= 20 && INFO/DP >= 5' \
    --soft-filter LowQual \
    --output-type z \
    -o "$FILTERED_VCF" \
    "$RAW_VCF"

bcftools index "$FILTERED_VCF"

# --------------------------------------------------------
# Normalize variants
# --------------------------------------------------------

echo
echo "Normalizing variants..."

bcftools norm \
    --fasta-ref "$REF" \
    --multiallelics -any \
    --output-type z \
    -o "$FINAL_VCF" \
    "$FILTERED_VCF"

bcftools index "$FINAL_VCF"

# --------------------------------------------------------
# Variant statistics
# --------------------------------------------------------

echo
echo "Variant summary:"

bcftools stats "$FINAL_VCF" | grep "^SN"

STATS_OUT="$LOG_DIR/${SAMPLE_NAME}_bcftools_stats.txt"

bcftools stats "$FINAL_VCF" > "$STATS_OUT"

echo
echo "Full stats saved to: $STATS_OUT"

echo
echo "--------------------------------------------"
echo "Variant calling completed successfully"
echo "Raw VCF      : $RAW_VCF"
echo "Filtered VCF : $FILTERED_VCF"
echo "Final VCF    : $FINAL_VCF"
echo "--------------------------------------------"
echo
