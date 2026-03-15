#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------
# PacBio HiFi Variant Pipeline — Local Runner
# --------------------------------------------------------
# Runs the full pipeline in a single command:
#   1. Alignment        (minimap2)
#   2. Variant calling  (bcftools)
#   3. QC report        (samtools + bcftools + seqkit)
#
# Usage:
#   bash scripts/run_pipeline.sh \
#       -r reference/chr20.fa \
#       -i test_data/synthetic_reads.fastq \
#       -s sample_name \
#       -o results \
#       [-t threads] [-p reports_dir]
# --------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/pipeline"

usage() {
    echo
    echo "Usage:"
    echo "  $0 -r reference.fa -i reads.fastq -s sample_name -o outdir [-t threads] [-p reports_dir]"
    echo
    echo "Options:"
    echo "  -r  Reference FASTA"
    echo "  -i  Input FASTQ (HiFi reads)"
    echo "  -s  Sample name (used for output file naming)"
    echo "  -o  Output directory for BAM and VCF files (default: ./results)"
    echo "  -t  Threads (default: 4)"
    echo "  -p  Reports directory (default: ./reports)"
    echo
    echo "Example:"
    echo "  bash scripts/run_pipeline.sh \\"
    echo "      -r reference/chr20.fa \\"
    echo "      -i test_data/synthetic_reads.fastq \\"
    echo "      -s synthetic_reads \\"
    echo "      -o results \\"
    echo "      -t 4"
    echo
    exit 1
}

# Defaults
THREADS=4
OUTDIR="./results"
REPORT_DIR="./reports"

while getopts ":r:i:s:o:t:p:" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        i) READS="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) REPORT_DIR="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG"; usage ;;
    esac
done

# --------------------------------------------------------
# Argument validation
# --------------------------------------------------------

if [[ -z "${REF:-}" || -z "${READS:-}" || -z "${SAMPLE:-}" ]]; then
    usage
fi

if [[ ! -f "$REF" ]]; then
    echo "Error: Reference FASTA not found: $REF"
    exit 1
fi

if [[ ! -f "$READS" ]]; then
    echo "Error: Input FASTQ not found: $READS"
    exit 1
fi

# --------------------------------------------------------
# Pipeline script validation
# --------------------------------------------------------

for script in align_minimap2.sh call_bcftools.sh qc_summary.sh; do
    if [[ ! -f "${PIPELINE_DIR}/${script}" ]]; then
        echo "Error: Missing pipeline script: ${PIPELINE_DIR}/${script}"
        exit 1
    fi
done

# --------------------------------------------------------
# Reference index
# --------------------------------------------------------

if [[ ! -f "${REF}.fai" ]]; then
    echo "Reference index not found. Creating with samtools faidx..."
    samtools faidx "$REF"
fi

# --------------------------------------------------------
# Setup
# --------------------------------------------------------

mkdir -p "$OUTDIR"
mkdir -p "$REPORT_DIR"

# Centralise logs — exported so all subscripts inherit the same location
export LOG_DIR="${OUTDIR}/logs"
mkdir -p "$LOG_DIR"

OUTPREFIX="${OUTDIR}/${SAMPLE}"
START_TIME=$(date +%s)

echo
echo "============================================"
echo "PacBio HiFi Variant Pipeline"
echo "============================================"
echo "Reference   : $REF"
echo "Reads       : $READS"
echo "Sample      : $SAMPLE"
echo "Output dir  : $OUTDIR"
echo "Report dir  : $REPORT_DIR"
echo "Log dir     : $LOG_DIR"
echo "Threads     : $THREADS"
echo "Started     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo

# --------------------------------------------------------
# Step 1 — Alignment
# --------------------------------------------------------

echo ">>> STEP 1/3 — Alignment (minimap2)"
echo

bash "${PIPELINE_DIR}/align_minimap2.sh" \
    -r "$REF" \
    -i "$READS" \
    -o "$OUTPREFIX" \
    -t "$THREADS"

echo
echo ">>> STEP 1/3 complete"
echo

# --------------------------------------------------------
# Step 2 — Variant calling
# --------------------------------------------------------

echo ">>> STEP 2/3 — Variant calling (bcftools)"
echo

bash "${PIPELINE_DIR}/call_bcftools.sh" \
    -r "$REF" \
    -b "${OUTPREFIX}.bam" \
    -o "$OUTPREFIX" \
    -t "$THREADS"

# Verify final VCF exists before proceeding to QC
if [[ ! -f "${OUTPREFIX}_norm.vcf.gz" ]]; then
    echo "Error: Final VCF not found after variant calling: ${OUTPREFIX}_norm.vcf.gz"
    exit 1
fi

echo
echo ">>> STEP 2/3 complete"
echo

# --------------------------------------------------------
# Step 3 — QC report
# --------------------------------------------------------

echo ">>> STEP 3/3 — QC report"
echo

bash "${PIPELINE_DIR}/qc_summary.sh" \
    -b "${OUTPREFIX}.bam" \
    -v "${OUTPREFIX}_norm.vcf.gz" \
    -r "$REF" \
    -o "$OUTPREFIX" \
    -p "$REPORT_DIR"

echo
echo ">>> STEP 3/3 complete"
echo

# --------------------------------------------------------
# Summary
# --------------------------------------------------------

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
MINUTES=$(( ELAPSED / 60 ))
SECONDS=$(( ELAPSED % 60 ))

echo "============================================"
echo "Pipeline completed successfully"
echo "--------------------------------------------"
echo "Sample      : $SAMPLE"
echo "BAM         : ${OUTPREFIX}.bam"
echo "VCF         : ${OUTPREFIX}_norm.vcf.gz"
echo "Logs        : ${LOG_DIR}/"
echo "Text report : ${REPORT_DIR}/${SAMPLE}_qc_report.txt"
echo "HTML report : ${REPORT_DIR}/${SAMPLE}_qc_report.html"
echo "Elapsed     : ${MINUTES}m ${SECONDS}s"
echo "============================================"
echo
