#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------
# PacBio HiFi QC Summary Script
# --------------------------------------------------------
# Generates QC metrics and reports from pipeline outputs.
#
# Coverage scaling strategy:
#   - Genomes < 200Mb  : full per-base depth (samtools depth)
#   - Genomes >= 200Mb : downsampled every 200 bases (1 in 200 positions)
#
# This keeps the coverage plot visually identical while reducing
# memory usage ~200x for WGS datasets.
# --------------------------------------------------------

LOG_DIR="${LOG_DIR:-"./pipeline_logs"}"
mkdir -p "$LOG_DIR"

exec > >(tee -a "$LOG_DIR/qc_general.log") \
     2> >(tee -a "$LOG_DIR/qc_error.log" >&2)

usage() {
    echo
    echo "Usage:"
    echo "  $0 -b aligned.bam -v variants.vcf.gz -r reference.fa -o output_prefix [-p reports_dir]"
    echo
    echo "Options:"
    echo "  -b  Sorted and indexed BAM file"
    echo "  -v  Final VCF file (compressed)"
    echo "  -r  Reference FASTA"
    echo "  -o  Output prefix for BAM/VCF files"
    echo "  -p  Reports directory (default: ./reports)"
    echo
    exit 1
}

while getopts ":b:v:r:o:p:" opt; do
    case $opt in
        b) BAM="$OPTARG" ;;
        v) VCF="$OPTARG" ;;
        r) REF="$OPTARG" ;;
        o) OUTPREFIX="$OPTARG" ;;
        p) REPORT_DIR="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG"; usage ;;
    esac
done

if [[ -z "${BAM:-}" || -z "${VCF:-}" || -z "${REF:-}" || -z "${OUTPREFIX:-}" ]]; then
    usage
fi

# Default reports directory
REPORT_DIR="${REPORT_DIR:-"./reports"}"
mkdir -p "$REPORT_DIR"

# --------------------------------------------------------
# Tool checks
# --------------------------------------------------------

for tool in samtools bcftools seqkit python3; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool not found in PATH"
        exit 1
    fi
done

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

if [[ ! -f "$VCF" ]]; then
    echo "Error: VCF file not found: $VCF"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "Error: Reference FASTA not found: $REF"
    exit 1
fi

# Auto-generate VCF index if missing
if [[ ! -f "${VCF}.csi" && ! -f "${VCF}.tbi" ]]; then
    echo "VCF index not found. Creating with bcftools index..."
    bcftools index "$VCF"
fi

# --------------------------------------------------------
# File naming
# --------------------------------------------------------

SAMPLE=$(basename "$BAM" .bam)
QC_TXT="${REPORT_DIR}/${SAMPLE}_qc_report.txt"
QC_HTML="${REPORT_DIR}/${SAMPLE}_qc_report.html"
METRICS_JSON="${REPORT_DIR}/${SAMPLE}_metrics.json"
COVERAGE_TSV="${REPORT_DIR}/${SAMPLE}_coverage.tsv"

echo
echo "--------------------------------------------"
echo "PacBio HiFi QC Summary"
echo "BAM        : $BAM"
echo "VCF        : $VCF"
echo "Reference  : $REF"
echo "Sample     : $SAMPLE"
echo "Report dir : $REPORT_DIR"
echo "--------------------------------------------"
echo

# --------------------------------------------------------
# Collect metrics
# --------------------------------------------------------

echo "Collecting alignment metrics..."
ALIGN_STATS=$(samtools flagstat "$BAM")

# Extract mapping rate for JSON
MAPPING_RATE=$(echo "$ALIGN_STATS" | grep "primary mapped" | grep -oE '[0-9]+\.[0-9]+%' || echo "NA")

echo "Collecting coverage metrics..."
COV_STATS=$(samtools coverage "$BAM")

echo "Collecting read metrics..."
READ_STATS=$(samtools stats "$BAM" | grep "^SN")

echo "Collecting variant metrics..."
VARIANT_COUNT=$(bcftools view -H "$VCF" | wc -l | tr -d ' ')

# --------------------------------------------------------
# Ti/Tv with fallback
# TSTV line may be absent if no SNPs were called
# (e.g. indel-only VCF from some callers)
# --------------------------------------------------------
TITV=$(bcftools stats "$VCF" | grep "^TSTV" | awk '{print $5}' || echo "NA")
if [[ -z "$TITV" ]]; then
    TITV="NA"
fi

VARIANT_SUMMARY=$(bcftools stats "$VCF" | grep "^SN")

# --------------------------------------------------------
# Coverage depth — hybrid downsampling strategy
#
# Genome size is determined from the BAM index stats.
# For genomes >= 200Mb, coverage is downsampled every 200
# bases to keep file size and Python memory usage manageable
# while preserving the visual shape of the coverage plot.
#
# Downsampling factors:
#   Small region / chr  : every base    (~70KB for chr20 region)
#   Full chromosome     : every 200bp   (~700KB for chr20)
#   Human WGS           : every 200bp   (~200MB → ~1MB)
# --------------------------------------------------------

echo "Determining genome size..."
GENOME_SIZE=$(samtools idxstats "$BAM" | awk '{sum+=$2} END {print sum}')
echo "Genome size: ${GENOME_SIZE}bp"

LARGE_GENOME_THRESHOLD=200000000  # 200Mb

if (( GENOME_SIZE > LARGE_GENOME_THRESHOLD )); then
    echo "Large genome detected (${GENOME_SIZE}bp > ${LARGE_GENOME_THRESHOLD}bp)"
    echo "Generating downsampled coverage depth (every 200 bases)..."
    samtools depth "$BAM" | awk 'NR % 200 == 0' > "$COVERAGE_TSV"
    COVERAGE_NOTE="downsampled (1 in 200 bases) — genome size: ${GENOME_SIZE}bp"
else
    echo "Generating per-base coverage depth..."
    samtools depth "$BAM" > "$COVERAGE_TSV"
    COVERAGE_NOTE="full per-base coverage — genome size: ${GENOME_SIZE}bp"
fi

echo "Coverage TSV: $COVERAGE_TSV"

# --------------------------------------------------------
# Write metrics JSON for HTML report
# --------------------------------------------------------

cat > "$METRICS_JSON" <<EOF
{
    "sample": "$SAMPLE",
    "mapping_rate": "$MAPPING_RATE",
    "variant_count": $VARIANT_COUNT,
    "titv": "$TITV",
    "genome_size": $GENOME_SIZE,
    "coverage_note": "$COVERAGE_NOTE"
}
EOF

# --------------------------------------------------------
# Write text QC report
# --------------------------------------------------------

{
echo "============================================"
echo "PacBio HiFi Pipeline QC Report"
echo "Sample    : $SAMPLE"
echo "Generated : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo

echo "--------------------------------------------"
echo "ALIGNMENT METRICS (samtools flagstat)"
echo "--------------------------------------------"
echo "$ALIGN_STATS"
echo

echo "--------------------------------------------"
echo "COVERAGE METRICS (samtools coverage)"
echo "--------------------------------------------"
echo "$COV_STATS"
echo

echo "--------------------------------------------"
echo "READ METRICS (samtools stats)"
echo "--------------------------------------------"
echo "$READ_STATS"
echo

echo "--------------------------------------------"
echo "VARIANT SUMMARY (bcftools stats)"
echo "--------------------------------------------"
echo "Total variants : $VARIANT_COUNT"
echo "Ti/Tv ratio    : $TITV"
echo
echo "$VARIANT_SUMMARY"
echo

echo "--------------------------------------------"
echo "VARIANT PREVIEW (first 10)"
echo "--------------------------------------------"
printf "%-10s %-10s %-6s %-6s %-8s %-10s\n" "CHROM" "POS" "REF" "ALT" "QUAL" "FILTER"
bcftools view -H "$VCF" | head -10 | awk '{printf "%-10s %-10s %-6s %-6s %-8s %-10s\n", $1,$2,$4,$5,$6,$7}'
echo

echo "--------------------------------------------"
echo "COVERAGE NOTE"
echo "--------------------------------------------"
echo "$COVERAGE_NOTE"

} > "$QC_TXT"

echo
echo "Text report saved: $QC_TXT"

# --------------------------------------------------------
# Generate HTML report with plots
# --------------------------------------------------------

echo
echo "Generating HTML report with plots..."

python3 scripts/report/generate_html_report.py \
    --sample      "$SAMPLE" \
    --metrics     "$METRICS_JSON" \
    --qc          "$QC_TXT" \
    --coverage    "$COVERAGE_TSV" \
    --vcf         "$VCF" \
    --output      "$QC_HTML"

echo "HTML report saved: $QC_HTML"

echo
echo "--------------------------------------------"
echo "QC summary completed successfully"
echo "Text report   : $QC_TXT"
echo "HTML report   : $QC_HTML"
echo "Coverage TSV  : $COVERAGE_TSV"
echo "Coverage      : $COVERAGE_NOTE"
echo "--------------------------------------------"
echo
