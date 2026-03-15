#!/usr/bin/env bash

set -euo pipefail

# --------------------------------------------------------
# PacBio HiFi Alignment Script (minimap2)
# --------------------------------------------------------

LOG_DIR="${LOG_DIR:-"./pipeline_logs"}"
mkdir -p "$LOG_DIR"

exec > >(tee -a "$LOG_DIR/alignment_general.log") \
     2> >(tee -a "$LOG_DIR/alignment_error.log" >&2)

usage() {
    echo
    echo "Usage:"
    echo "  $0 -r reference.fa -i reads.fastq -o output_prefix [-t threads]"
    echo
    echo "Options:"
    echo "  -r  Reference FASTA"
    echo "  -i  Input FASTQ (HiFi reads)"
    echo "  -o  Output prefix"
    echo "  -t  Threads (default: 4)"
    echo
    exit 1
}

THREADS=4

while getopts ":r:i:o:t:" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        i) READS="$OPTARG" ;;
        o) OUTPREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG"; usage ;;
    esac
done

if [[ -z "${REF:-}" || -z "${READS:-}" || -z "${OUTPREFIX:-}" ]]; then
    usage
fi

if ! command -v minimap2 &> /dev/null; then
    echo "Error: minimap2 not found in PATH"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found in PATH"
    exit 1
fi

SAMPLE_NAME=$(basename "$READS" .fastq)
SAMPLE_NAME=$(basename "$SAMPLE_NAME" .fq)
SAMPLE_NAME=$(basename "$SAMPLE_NAME" .fq.gz)

echo
echo "--------------------------------------------"
echo "PacBio HiFi Alignment (minimap2)"
echo "Sample    : $SAMPLE_NAME"
echo "Reference : $REF"
echo "Reads     : $READS"
echo "Threads   : $THREADS"
echo "Output    : ${OUTPREFIX}.bam"
echo "--------------------------------------------"
echo

# Skip if output exists
if [[ -f "${OUTPREFIX}.bam" ]]; then
    echo "Output BAM already exists. Skipping alignment."
    exit 0
fi

echo "Running minimap2..."

set +e
minimap2 \
    -t "$THREADS" \
    -ax map-hifi \
    --secondary=no \
    -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:PACBIO" \
    "$REF" \
    "$READS" \
| samtools view -b - \
| samtools sort -@ "$THREADS" -o "${OUTPREFIX}.bam"

PIPE_STATUS=${PIPESTATUS[0]}
set -e

if [[ $PIPE_STATUS -ne 0 ]]; then
    echo "Error: minimap2 alignment failed"
    exit 1
fi

echo
echo "Indexing BAM..."

samtools index "${OUTPREFIX}.bam"

echo
echo "Alignment statistics"

FLAGSTAT_OUT="$LOG_DIR/$(basename ${OUTPREFIX})_flagstat.txt"
samtools flagstat "${OUTPREFIX}.bam" | tee "$FLAGSTAT_OUT"

echo
echo "Flagstat saved to: $FLAGSTAT_OUT"

echo
echo "Alignment completed successfully"
echo "Output: ${OUTPREFIX}.bam"
echo
