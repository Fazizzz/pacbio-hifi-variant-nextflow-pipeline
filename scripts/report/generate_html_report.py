#!/usr/bin/env python3

"""
Generate an HTML QC report for the PacBio HiFi variant pipeline.

Reads pre-computed metrics from a JSON file produced by qc_summary.sh
and formats them into a clean, human-readable HTML report with
embedded coverage and variant position plots.

Usage:
    python3 generate_html_report.py \
        --sample   sample_name \
        --metrics  reports/sample_metrics.json \
        --qc       reports/sample_qc_report.txt \
        --coverage reports/sample_coverage.tsv \
        --vcf      test_data/variants_norm.vcf.gz \
        --output   reports/sample_qc_report.html

Dependencies:
    matplotlib (pip install matplotlib)
"""

import argparse
import base64
import io
import json
import subprocess
from datetime import datetime
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend, no display required
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# --------------------------------------------------------
# Data loading
# --------------------------------------------------------

def load_metrics(metrics_path):
    """Load pre-computed metrics from JSON file."""
    with open(metrics_path) as f:
        return json.load(f)


def load_qc_text(qc_path):
    """Load raw QC text report."""
    with open(qc_path) as f:
        return f.read()


def load_coverage(coverage_tsv):
    """
    Load per-base (or downsampled) coverage from samtools depth output.

    Returns:
        positions: list of int
        depths:    list of int
    """
    positions = []
    depths = []

    with open(coverage_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                positions.append(int(parts[1]))
                depths.append(int(parts[2]))
            except ValueError:
                continue

    return positions, depths


def load_variants(vcf_path):
    """
    Extract variant positions from a compressed VCF using bcftools.

    Returns:
        list of (chrom, pos, filter) tuples
    """
    try:
        result = subprocess.check_output(
            f"bcftools view -H {vcf_path}",
            shell=True,
            text=True
        )
        variants = []
        for line in result.strip().splitlines():
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            chrom = parts[0]
            pos   = int(parts[1])
            filt  = parts[6]
            variants.append((chrom, pos, filt))
        return variants
    except subprocess.CalledProcessError:
        return []


# --------------------------------------------------------
# Plot generation
# --------------------------------------------------------

def fig_to_base64(fig):
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return encoded


def plot_coverage(positions, depths, sample, coverage_note=""):
    """
    Generate a coverage depth plot across the sequenced region.

    Handles both full per-base and downsampled coverage TSVs transparently.
    The plot appearance is identical regardless of downsampling factor.

    Returns base64-encoded PNG string.
    """
    fig, ax = plt.subplots(figsize=(12, 3.5))

    ax.fill_between(positions, depths, alpha=0.4, color="#1a6faf", linewidth=0)
    ax.plot(positions, depths, color="#1a6faf", linewidth=0.8)

    if depths:
        mean_depth = sum(depths) / len(depths)
        ax.axhline(
            mean_depth,
            color="#e05c2a",
            linewidth=1.2,
            linestyle="--",
            label=f"Mean depth: {mean_depth:.1f}x"
        )
        ax.legend(fontsize=9)

    # Show coverage note under x-axis if downsampled
    ax.set_title(f"{sample} — Coverage across region", fontsize=11, fontweight="bold")

    if coverage_note:
        fig.text(
            0.5, -0.02,
            f"Coverage note: {coverage_note}",
            ha="center",
            fontsize=8,
            color="#888",
            style="italic"
        )

    ax.set_xlabel("Genomic position", fontsize=10)
    ax.set_ylabel("Read depth", fontsize=10)
    ax.set_xlim(min(positions), max(positions)) if positions else None
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linewidth=0.5)

    fig.tight_layout()
    return fig_to_base64(fig)


def plot_variant_positions(variants, positions, sample):
    """
    Generate a variant position scatter plot across the reference region.

    Returns base64-encoded PNG string.
    """
    fig, ax = plt.subplots(figsize=(12, 2.2))

    if positions:
        region_start = min(positions)
        region_end   = max(positions)
    else:
        region_start = min(v[1] for v in variants) if variants else 0
        region_end   = max(v[1] for v in variants) if variants else 1

    # Reference region bar
    ax.barh(
        0, region_end - region_start,
        left=region_start,
        height=0.15,
        color="#cccccc",
        zorder=1
    )

    # Color variants by filter status
    pass_positions     = [v[1] for v in variants if v[2] == "PASS"]
    filtered_positions = [v[1] for v in variants if v[2] != "PASS"]

    if pass_positions:
        ax.scatter(
            pass_positions, [0] * len(pass_positions),
            color="#2e7d32", s=80, zorder=3,
            marker="|", linewidths=2.5, label="PASS"
        )

    if filtered_positions:
        ax.scatter(
            filtered_positions, [0] * len(filtered_positions),
            color="#c62828", s=80, zorder=3,
            marker="|", linewidths=2.5, label="Filtered (LowQual)"
        )

    ax.set_xlabel("Genomic position", fontsize=10)
    ax.set_title(
        f"{sample} — Variant positions ({len(variants)} total)",
        fontsize=11, fontweight="bold"
    )
    ax.set_xlim(region_start - 200, region_end + 200)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    if pass_positions or filtered_positions:
        ax.legend(fontsize=9, loc="upper right")

    fig.tight_layout()
    return fig_to_base64(fig)


# --------------------------------------------------------
# HTML generation
# --------------------------------------------------------

def generate_html(sample, metrics, qc_text, coverage_plot_b64,
                  variant_plot_b64, output):
    """Render self-contained HTML report."""

    mapping_rate  = metrics.get("mapping_rate",  "NA")
    variant_count = metrics.get("variant_count", "NA")
    titv          = metrics.get("titv",          "NA")
    genome_size   = metrics.get("genome_size",   "NA")
    coverage_note = metrics.get("coverage_note", "")
    generated     = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Format genome size with commas for readability
    if isinstance(genome_size, int):
        genome_size_fmt = f"{genome_size:,}bp"
    else:
        genome_size_fmt = str(genome_size)

    # Build plot sections — gracefully omit if plots unavailable
    coverage_section = ""
    if coverage_plot_b64:
        coverage_section = f"""
<h2>Coverage Across Region</h2>
<img src="data:image/png;base64,{coverage_plot_b64}"
     alt="Coverage plot" style="width:100%; max-width:960px;">
"""

    variant_section = ""
    if variant_plot_b64:
        variant_section = f"""
<h2>Variant Positions</h2>
<img src="data:image/png;base64,{variant_plot_b64}"
     alt="Variant positions plot" style="width:100%; max-width:960px;">
"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PacBio HiFi Pipeline QC Report — {sample}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 40px;
            max-width: 960px;
            color: #333;
        }}
        h1 {{
            color: #1a1a2e;
            border-bottom: 2px solid #1a1a2e;
            padding-bottom: 8px;
        }}
        h2 {{
            color: #16213e;
            margin-top: 36px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin-top: 8px;
        }}
        th {{
            background-color: #1a1a2e;
            color: white;
            padding: 8px 12px;
            text-align: left;
        }}
        td {{
            border: 1px solid #ddd;
            padding: 8px 12px;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        .metric-value {{
            font-weight: bold;
            color: #1a1a2e;
        }}
        .note {{
            font-size: 12px;
            color: #888;
            font-style: italic;
        }}
        img {{
            border: 1px solid #ddd;
            border-radius: 4px;
            margin-top: 8px;
        }}
        pre {{
            background: #f4f4f4;
            padding: 16px;
            overflow-x: auto;
            max-height: 500px;
            overflow-y: scroll;
            font-size: 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            white-space: pre-wrap;
            word-wrap: break-word;
        }}
        .footer {{
            margin-top: 48px;
            font-size: 12px;
            color: #888;
            border-top: 1px solid #ddd;
            padding-top: 12px;
        }}
    </style>
</head>
<body>

<h1>PacBio HiFi Variant Pipeline QC Report</h1>

<table>
    <tr><td><b>Sample</b></td><td class="metric-value">{sample}</td></tr>
    <tr><td><b>Generated</b></td><td>{generated}</td></tr>
    <tr><td><b>Genome size</b></td><td>{genome_size_fmt}</td></tr>
</table>

<h2>Alignment Summary</h2>
<table>
    <tr><th>Metric</th><th>Value</th></tr>
    <tr><td>Primary mapping rate</td><td class="metric-value">{mapping_rate}</td></tr>
</table>

<h2>Variant Summary</h2>
<table>
    <tr><th>Metric</th><th>Value</th></tr>
    <tr><td>Total variants (post-filter)</td><td class="metric-value">{variant_count}</td></tr>
    <tr><td>Ti/Tv ratio</td><td class="metric-value">{titv}</td></tr>
</table>

{coverage_section}
{"<p class='note'>Coverage note: " + coverage_note + "</p>" if coverage_note else ""}

{variant_section}

<h2>Full QC Report</h2>
<pre>{qc_text}</pre>

<div class="footer">
    Generated by PacBio HiFi Variant Pipeline &mdash; {generated}
</div>

</body>
</html>
"""

    with open(output, "w") as f:
        f.write(html)


# --------------------------------------------------------
# Main
# --------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description="Generate HTML QC report for PacBio HiFi variant pipeline."
    )

    parser.add_argument("--sample",   required=True, help="Sample name")
    parser.add_argument("--metrics",  required=True, help="Path to metrics JSON")
    parser.add_argument("--qc",       required=True, help="Path to QC text report")
    parser.add_argument("--coverage", required=True, help="Path to samtools depth TSV")
    parser.add_argument("--vcf",      required=True, help="Path to final VCF (compressed)")
    parser.add_argument("--output",   required=True, help="Path to write HTML report")

    args = parser.parse_args()

    # Validate inputs
    for label, path in [
        ("Metrics JSON",    args.metrics),
        ("QC text report",  args.qc),
        ("Coverage TSV",    args.coverage),
        ("VCF",             args.vcf),
    ]:
        if not Path(path).exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    # Load data
    metrics           = load_metrics(args.metrics)
    qc_text           = load_qc_text(args.qc)
    positions, depths = load_coverage(args.coverage)
    variants          = load_variants(args.vcf)
    coverage_note     = metrics.get("coverage_note", "")

    # Generate plots
    coverage_plot_b64 = None
    variant_plot_b64  = None

    if not HAS_MATPLOTLIB:
        print("Warning: matplotlib not found. Plots will be skipped.")
        print("Install with: pip install matplotlib")
    else:
        if positions:
            print("Generating coverage plot...")
            coverage_plot_b64 = plot_coverage(
                positions, depths, args.sample, coverage_note
            )
        else:
            print("Warning: No coverage data found. Coverage plot skipped.")

        if variants:
            print("Generating variant position plot...")
            variant_plot_b64 = plot_variant_positions(
                variants, positions, args.sample
            )
        else:
            print("Warning: No variants found. Variant plot skipped.")

    # Generate report
    generate_html(
        sample=args.sample,
        metrics=metrics,
        qc_text=qc_text,
        coverage_plot_b64=coverage_plot_b64,
        variant_plot_b64=variant_plot_b64,
        output=args.output
    )

    print(f"HTML report written to: {args.output}")


if __name__ == "__main__":
    main()
