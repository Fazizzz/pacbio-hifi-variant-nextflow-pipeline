#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------
# render_dag.sh
# --------------------------------------------------------
# Renders a clean, readable SVG from a Nextflow DAG dot file.
#
# Nextflow's default DAG output has overlapping text and
# cluttered edge labels. This script:
#   - Strips the workflow prefix from node labels
#   - Applies left-to-right layout for better readability
#   - Increases node spacing and font size
#   - Outputs a clean SVG suitable for documentation and README
#
# Requirements:
#   graphviz (conda install -c conda-forge graphviz=12.2.*)
#
# Usage:
#   bash scripts/report/render_dag.sh
#   bash scripts/report/render_dag.sh --dot path/to/dag.dot --out path/to/dag_clean.svg
#
# Default paths:
#   input  : dag.dot      (repo root)
#   output : dag_clean.svg (repo root)
#
# Generate the input dot file with:
#   nextflow run main.nf -profile test -stub-run -with-dag dag.dot
# --------------------------------------------------------

# --------------------------------------------------------
# Defaults
# --------------------------------------------------------
DOT_FILE="dag.dot"
OUT_FILE="dag_clean.svg"

# --------------------------------------------------------
# Argument parsing
# --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dot)
            DOT_FILE="$2"
            shift 2
            ;;
        --out)
            OUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            echo
            echo "Usage: $0 [--dot dag.dot] [--out dag_clean.svg]"
            echo
            echo "Options:"
            echo "  --dot  Path to input Nextflow DAG dot file (default: dag.dot)"
            echo "  --out  Path to output SVG file (default: dag_clean.svg)"
            echo
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# --------------------------------------------------------
# Tool check
# --------------------------------------------------------
if ! command -v dot &> /dev/null; then
    echo "Error: graphviz not found in PATH"
    echo "Install with: conda install -c conda-forge graphviz=12.2.*"
    exit 1
fi

# --------------------------------------------------------
# Input validation
# --------------------------------------------------------
if [[ ! -f "$DOT_FILE" ]]; then
    echo "Error: DOT file not found: $DOT_FILE"
    echo
    echo "Generate it first with:"
    echo "  nextflow run main.nf -profile test -stub-run -with-dag dag.dot"
    exit 1
fi

echo
echo "--------------------------------------------"
echo "Rendering pipeline DAG"
echo "Input  : $DOT_FILE"
echo "Output : $OUT_FILE"
echo "--------------------------------------------"
echo

# --------------------------------------------------------
# Render
# Strip workflow prefix from node labels then apply
# layout parameters for a clean readable diagram.
#
# -Grankdir=LR       left-to-right layout
# -Gnodesep=1.0      horizontal spacing between nodes
# -Granksep=3        vertical spacing between ranks
# -Nwidth=2.5        node width
# -Nheight=1         node height
# -Nfontsize=13      node label font size
# -Efontsize=13      edge label font size
# -Elabelfloat=false keep edge labels anchored
# -Gdecorate=false   cleaner edge decoration
# --------------------------------------------------------
sed 's/PACBIO_VARIANT_WORKFLOW://g' "$DOT_FILE" | \
    dot -Tsvg \
        -Grankdir=LR \
        -Gnodesep=1.0 \
        -Granksep=3 \
        -Nwidth=2.5 \
        -Nheight=1 \
        -Nfontsize=13 \
        -Efontsize=13 \
        -Elabelfloat=false \
        -Gdecorate=false \
    -o "$OUT_FILE"

echo "DAG rendered successfully: $OUT_FILE"
echo
