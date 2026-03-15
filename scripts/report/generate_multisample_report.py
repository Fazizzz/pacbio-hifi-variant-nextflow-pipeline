#!/usr/bin/env python3

"""
Generate a harmonized multi-sample summary HTML report for the
PacBio HiFi variant calling pipeline.

Layout:
    1. Sticky top bar
    2. Legend bar
    3. Summary table (thead sticks below top bar + legend)
    4. Variant position track

Variant track behaviour:
    - Single chromosome <= 10Mb  : per-variant marker track
    - Single chromosome > 10Mb   : binned density track (500 windows)
    - Multiple chromosomes        : variant count bar chart with note

Usage:
    python3 generate_multisample_report.py \\
        --reports "results_multisample/reports/*_qc_report.txt" \\
        --output  results_multisample/reports/multisample_summary.html
"""

import argparse
import glob
import html as html_module
import re
from datetime import datetime
from pathlib import Path


# --------------------------------------------------------
# Parsing
# --------------------------------------------------------

def parse_qc_report(report_path):
    try:
        text = Path(report_path).read_text()
    except Exception as e:
        print(f"Warning: Could not read {report_path}: {e}")
        return None

    metrics = {
        "sample":         "unknown",
        "generated":      "unknown",
        "total_reads":    "NA",
        "mapped_reads":   "NA",
        "mapping_rate":   "NA",
        "mean_length":    "NA",
        "avg_quality":    "NA",
        "mean_depth":     "NA",
        "total_variants": "NA",
        "snps":           "NA",
        "indels":         "NA",
        "pass_variants":  "NA",
        "variants":       [],
        "report_path":    str(report_path),
        "raw_text":       html_module.escape(text),
    }

    for pattern, key in [
        (r"Sample\s*:\s*(.+)",          "sample"),
        (r"Generated\s*:\s*(.+)",       "generated"),
        (r"number of records:\s+(\d+)", "total_variants"),
        (r"number of SNPs:\s+(\d+)",    "snps"),
        (r"number of indels:\s+(\d+)",  "indels"),
    ]:
        m = re.search(pattern, text)
        if m:
            metrics[key] = m.group(1).strip()

    m = re.search(r"(\d+) \+ \d+ in total", text)
    if m:
        metrics["total_reads"] = m.group(1)

    m = re.search(r"(\d+) \+ \d+ primary mapped \((\d+\.\d+)%", text)
    if m:
        metrics["mapped_reads"] = m.group(1)
        metrics["mapping_rate"] = m.group(2)

    m = re.search(r"average length:\s+(\d+)", text)
    if m:
        metrics["mean_length"] = m.group(1)

    m = re.search(r"average quality:\s+([\d.]+)", text)
    if m:
        metrics["avg_quality"] = m.group(1)

    m = re.search(r"chr\S+\s+\d+\s+\d+\s+\d+\s+\d+\s+[\d.]+\s+([\d.]+)", text)
    if m:
        metrics["mean_depth"] = m.group(1)

    preview_section = re.search(
        r"VARIANT PREVIEW.*?\n(.*?)(?:\n\n|\Z)", text, re.DOTALL
    )
    if preview_section:
        pass_count = 0
        for line in preview_section.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 6 and parts[0].startswith("chr"):
                filt = parts[5]
                metrics["variants"].append({
                    "chrom":  parts[0],
                    "pos":    int(parts[1]),
                    "ref":    parts[2],
                    "alt":    parts[3],
                    "qual":   parts[4],
                    "filter": filt,
                })
                if filt == "PASS":
                    pass_count += 1
        metrics["pass_variants"] = str(pass_count)

    return metrics


def load_reports(report_paths):
    samples = []
    for path in report_paths:
        parsed = parse_qc_report(path)
        if parsed:
            samples.append(parsed)
    return sorted(samples, key=lambda x: x["sample"])


# --------------------------------------------------------
# Variant track SVG
# --------------------------------------------------------

LABEL_W = 150
BAR_W   = 640
SVG_W   = LABEL_W + BAR_W + 30
ROW_H   = 34
ROW_GAP = 10
TOP_PAD = 20
AXIS_H  = 44


def _svg_h(n):
    return TOP_PAD + n * (ROW_H + ROW_GAP) + AXIS_H


def _xp(pos, r_start, r_size):
    return LABEL_W + int((pos - r_start) / r_size * BAR_W)


def build_variant_track(samples):
    all_v = [v for s in samples for v in s["variants"]]
    if not all_v:
        return ""

    chroms = list({v["chrom"] for v in all_v})
    if len(chroms) > 1:
        return _count_bars(samples, chroms)

    positions = [v["pos"] for v in all_v]
    r_start_raw = min(positions)
    r_end_raw   = max(positions)
    r_size_raw  = max(r_end_raw - r_start_raw, 1)
    pad         = int(r_size_raw * 0.05)
    r_start     = max(0, r_start_raw - pad)
    r_end       = r_end_raw + pad
    r_size      = r_end - r_start
    chrom       = chroms[0]

    if r_size > 10_000_000:
        return _density_track(samples, chrom, r_start, r_end, r_size)
    return _marker_track(samples, chrom, r_start, r_end, r_size)


def _marker_track(samples, chrom, r_start, r_end, r_size):
    n    = len(samples)
    svgh = _svg_h(n)
    body = ""
    for i, s in enumerate(samples):
        mid  = TOP_PAD + i * (ROW_H + ROW_GAP) + ROW_H // 2
        top  = mid - ROW_H // 2 + 5
        bot  = mid + ROW_H // 2 - 5
        body += (
            f'<text x="{LABEL_W-10}" y="{mid+4}" text-anchor="end" '
            f'font-size="11" fill="#555" font-family="Arial">{html_module.escape(s["sample"])}</text>\n'
            f'<rect x="{LABEL_W}" y="{mid-2}" width="{BAR_W}" height="4" rx="2" fill="#ddd"/>\n'
        )
        for v in s["variants"]:
            vx    = _xp(v["pos"], r_start, r_size)
            color = "#2e7d32" if v["filter"] == "PASS" else "#e65100"
            tip   = html_module.escape(
                f'{v["chrom"]}:{v["pos"]} {v["ref"]}→{v["alt"]} ({v["filter"]})'
            )
            body += (
                f'<line x1="{vx}" y1="{top}" x2="{vx}" y2="{bot}" '
                f'stroke="{color}" stroke-width="2.5" stroke-linecap="round" opacity="0.85">'
                f'<title>{tip}</title></line>\n'
            )
    body += _ticks(r_start, r_end, r_size, svgh)
    body += _marker_legend(svgh, chrom, r_start, r_end)
    return _wrap(body, SVG_W, svgh)


def _density_track(samples, chrom, r_start, r_end, r_size):
    n      = len(samples)
    svgh   = _svg_h(n)
    n_bins = 500
    bin_px = BAR_W / n_bins
    body   = ""
    for i, s in enumerate(samples):
        mid = TOP_PAD + i * (ROW_H + ROW_GAP) + ROW_H // 2
        body += (
            f'<text x="{LABEL_W-10}" y="{mid+4}" text-anchor="end" '
            f'font-size="11" fill="#555" font-family="Arial">{html_module.escape(s["sample"])}</text>\n'
            f'<rect x="{LABEL_W}" y="{mid-2}" width="{BAR_W}" height="4" rx="2" fill="#eee"/>\n'
        )
        bins = [0] * n_bins
        for v in s["variants"]:
            idx = min(int((v["pos"] - r_start) / r_size * n_bins), n_bins - 1)
            bins[idx] += 1
        mx = max(bins) if any(bins) else 1
        for b, count in enumerate(bins):
            if not count:
                continue
            bx = LABEL_W + int(b * bin_px)
            bh = max(2, int((count / mx) * (ROW_H - 10)))
            by = mid + ROW_H // 2 - 5 - bh
            body += (
                f'<rect x="{bx}" y="{by}" width="{max(1,int(bin_px))}" height="{bh}" '
                f'fill="#1a6faf" opacity="0.7">'
                f'<title>{count} variant(s) window {b+1}</title></rect>\n'
            )
    body += _ticks(r_start, r_end, r_size, svgh)
    rmb   = r_size / 1_000_000
    ny    = svgh - 6
    body += (
        f'<text x="{LABEL_W}" y="{ny}" font-size="10" fill="#aaa" font-family="Arial">'
        f'Density · {chrom}:{r_start:,}–{r_end:,} ({rmb:.1f} Mb) · '
        f'bar height = relative density per window</text>'
    )
    return _wrap(body, SVG_W, svgh)


def _count_bars(samples, chroms):
    n    = len(samples)
    svgh = _svg_h(n) + 20
    mx   = max((len(s["variants"]) for s in samples), default=1)
    body = ""
    for i, s in enumerate(samples):
        mid   = TOP_PAD + i * (ROW_H + ROW_GAP) + ROW_H // 2
        count = len(s["variants"])
        bw    = int((count / max(mx, 1)) * BAR_W)
        body += (
            f'<text x="{LABEL_W-10}" y="{mid+4}" text-anchor="end" '
            f'font-size="11" fill="#555" font-family="Arial">{html_module.escape(s["sample"])}</text>\n'
            f'<rect x="{LABEL_W}" y="{mid-ROW_H//2+4}" width="{max(bw,2)}" '
            f'height="{ROW_H-8}" rx="3" fill="#1a6faf" opacity="0.75"/>\n'
            f'<text x="{LABEL_W+max(bw,2)+6}" y="{mid+4}" '
            f'font-size="11" fill="#555" font-family="Arial">{count}</text>\n'
        )
    cl   = html_module.escape(", ".join(sorted(chroms)))
    body += (
        f'<text x="{LABEL_W}" y="{svgh-10}" font-size="10" fill="#e65100" font-family="Arial">'
        f'⚠ Positional track unavailable — variants span multiple chromosomes: {cl}</text>'
    )
    return _wrap(body, SVG_W, svgh)


def _ticks(r_start, r_end, r_size, svgh):
    ay  = svgh - AXIS_H + 8
    out = ""
    for i in range(7):
        pos = r_start + int(i * r_size / 6)
        vx  = _xp(pos, r_start, r_size)
        lbl = (
            f"{pos/1_000_000:.1f}Mb" if pos >= 1_000_000
            else f"{pos//1000}kb"    if pos >= 1_000
            else str(pos)
        )
        out += (
            f'<line x1="{vx}" y1="{ay}" x2="{vx}" y2="{ay+6}" stroke="#bbb" stroke-width="0.5"/>\n'
            f'<text x="{vx}" y="{ay+17}" text-anchor="middle" font-size="9" '
            f'fill="#aaa" font-family="Arial">{lbl}</text>\n'
        )
    return out


def _marker_legend(svgh, chrom, r_start, r_end):
    ly = svgh - 10
    return (
        f'<circle cx="{LABEL_W}" cy="{ly}" r="4" fill="#2e7d32"/>'
        f'<text x="{LABEL_W+8}" y="{ly+4}" font-size="10" fill="#555" font-family="Arial">PASS</text>'
        f'<circle cx="{LABEL_W+55}" cy="{ly}" r="4" fill="#e65100"/>'
        f'<text x="{LABEL_W+63}" y="{ly+4}" font-size="10" fill="#555" font-family="Arial">LowQual</text>'
        f'<text x="{LABEL_W+130}" y="{ly+4}" font-size="10" fill="#aaa" font-family="Arial">'
        f'Hover for details · {chrom}:{r_start:,}–{r_end:,}</text>'
    )


def _wrap(content, w, h):
    return (
        f'<div style="overflow-x:auto;margin:0 0 8px 0;">'
        f'<svg width="100%" viewBox="0 0 {w} {h}" style="min-width:500px;display:block;">'
        f'{content}</svg></div>'
    )


# --------------------------------------------------------
# Quality helpers
# --------------------------------------------------------

def qc(value, metric):
    try:
        v = float(value)
    except (ValueError, TypeError):
        return "na"
    if metric == "mapping_rate":
        return "pass" if v >= 85 else ("warn" if v >= 70 else "fail")
    if metric == "avg_quality":
        return "pass" if v >= 20 else ("warn" if v >= 15 else "fail")
    if metric == "mean_depth":
        return "pass" if v >= 15 else ("warn" if v >= 5 else "fail")
    return "neutral"


def cell(value, metric=""):
    return f'<td class="q-{qc(value, metric)}">{html_module.escape(str(value))}</td>'


# --------------------------------------------------------
# HTML
# --------------------------------------------------------

def generate_html(samples, output_path):
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n         = len(samples)
    n_pass    = sum(1 for s in samples if qc(s["mapping_rate"], "mapping_rate") == "pass")
    n_warn    = sum(1 for s in samples if qc(s["mapping_rate"], "mapping_rate") == "warn")
    n_fail    = sum(1 for s in samples if qc(s["mapping_rate"], "mapping_rate") == "fail")
    track     = build_variant_track(samples)

    rows = ""
    for i, s in enumerate(samples):
        sid   = f"s{i}"
        sname = html_module.escape(s["sample"])
        rows += f"""<tr class="sr" data-detail="{sid}">
  <td class="sn">{sname} <span class="ei">▼</span></td>
  <td>{html_module.escape(str(s['total_reads']))}</td>
  {cell(s['mapped_reads'])}
  {cell(s['mapping_rate'], 'mapping_rate')}
  <td>{html_module.escape(str(s['mean_length']))}</td>
  {cell(s['avg_quality'], 'avg_quality')}
  {cell(s['mean_depth'], 'mean_depth')}
  <td>{html_module.escape(str(s['total_variants']))}</td>
  <td>{html_module.escape(str(s['snps']))}</td>
  <td>{html_module.escape(str(s['indels']))}</td>
  <td>{html_module.escape(str(s['pass_variants']))}</td>
</tr>
<tr id="{sid}" style="display:none;"><td colspan="11"><div class="db">
  <div class="dm"><b>Sample:</b> {sname} &nbsp;|&nbsp;
  <b>Generated:</b> {html_module.escape(s['generated'])}</div>
  <pre class="dp">{s['raw_text']}</pre>
</div></td></tr>
"""

    track_section = (
        f'<div class="sec"><div class="st">Variant positions — all samples</div>{track}</div>'
        if track else ""
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>PacBio HiFi — Multi-Sample QC</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:Arial,sans-serif;font-size:13px;color:#333;background:#f8f8f8}}

/* ---- Sticky top bar ---- */
#topbar{{
  background:#1a1a2e;color:#fff;
  padding:14px 24px;
  position:sticky;top:0;z-index:200;
  display:flex;align-items:center;justify-content:space-between;
}}
#topbar h1{{font-size:15px;font-weight:700}}
.hm{{font-size:11px;color:#aaa;margin-top:2px}}
.bdg{{display:flex;gap:8px}}
.b{{padding:3px 10px;border-radius:12px;font-size:11px;font-weight:700}}
.bp{{background:#2e7d32;color:#fff}}.bw{{background:#f57c00;color:#fff}}
.bf{{background:#c62828;color:#fff}}.bt{{background:#444;color:#fff}}

/* ---- Legend bar (sticky, sits below top bar) ---- */
#legbar{{
  background:#fff;border-bottom:1px solid #eee;
  padding:7px 24px;
  position:sticky;top:53px;z-index:150;
  display:flex;gap:14px;align-items:center;flex-wrap:wrap;
  font-size:11px;color:#666;
}}
.li{{display:flex;align-items:center;gap:4px}}
.ld{{width:9px;height:9px;border-radius:50%;display:inline-block}}

/* ---- Table wrapper — NO overflow, let page scroll ---- */
.tw{{padding:0 24px 16px;padding-top:8px}}

table{{
  width:100%;border-collapse:collapse;background:#fff;
  box-shadow:0 1px 4px rgba(0,0,0,.08);
}}

/* ---- Sticky thead — applied to <thead> not <th> ---- */
thead{{
  position:sticky;
  top:0;   /* JS sets this to topbar+legbar height after load */
  z-index:100;
  box-shadow:0 2px 4px rgba(0,0,0,.15);
}}
thead th{{
  background:#1a1a2e;color:#fff;
  padding:10px 12px;text-align:left;
  font-size:11px;font-weight:700;
  text-transform:uppercase;letter-spacing:.5px;white-space:nowrap;
}}

.sr{{cursor:pointer;border-bottom:1px solid #eee;transition:background .1s}}
.sr:hover{{background:#f0f4ff}}
td{{padding:8px 12px;white-space:nowrap}}
.sn{{font-weight:700;color:#1a1a2e;min-width:160px}}
.ei{{font-size:10px;color:#888;margin-left:4px}}
.q-pass{{color:#2e7d32;font-weight:700}}.q-warn{{color:#e65100;font-weight:700}}
.q-fail{{color:#c62828;font-weight:700}}.q-neutral{{color:#333}}.q-na{{color:#aaa}}

.db{{background:#f4f4f4;border-left:3px solid #1a1a2e;
  padding:12px 16px;margin:4px 8px 8px;border-radius:4px}}
.dm{{font-size:11px;color:#555;margin-bottom:8px}}
.dp{{font-size:11px;white-space:pre-wrap;word-wrap:break-word;
  max-height:350px;overflow-y:auto;background:#fff;
  padding:10px;border:1px solid #ddd;border-radius:4px}}

.sec{{padding:16px 24px 8px}}
.st{{font-size:11px;font-weight:700;text-transform:uppercase;
  letter-spacing:.5px;color:#888;margin-bottom:8px}}
.ft{{padding:16px 24px;font-size:11px;color:#aaa;
  border-top:1px solid #ddd;margin-top:8px}}
</style>
</head>
<body>

<!-- 1. Sticky top bar -->
<div id="topbar">
  <div>
    <h1>PacBio HiFi Variant Pipeline — Multi-Sample QC Summary</h1>
    <div class="hm">Generated: {generated}</div>
  </div>
  <div class="bdg">
    <span class="b bt">{n} samples</span>
    <span class="b bp">{n_pass} pass</span>
    <span class="b bw">{n_warn} warn</span>
    <span class="b bf">{n_fail} fail</span>
  </div>
</div>

<!-- 2. Sticky legend bar (below top bar) -->
<div id="legbar">
  <span>Mapping rate:</span>
  <span class="li"><span class="ld" style="background:#2e7d32"></span>≥85% pass</span>
  <span class="li"><span class="ld" style="background:#e65100"></span>70–85% warn</span>
  <span class="li"><span class="ld" style="background:#c62828"></span>&lt;70% fail</span>
  <span style="margin-left:8px;color:#aaa">Click any row to expand full QC report</span>
</div>

<!-- 3. Summary table (no scroll container — page scrolls) -->
<div class="tw">
<table id="t">
  <thead>
    <tr>
      <th>Sample</th><th>Total reads</th><th>Mapped</th><th>Mapping rate (%)</th>
      <th>Mean length (bp)</th><th>Avg quality</th><th>Mean depth</th>
      <th>Total variants</th><th>SNPs</th><th>Indels</th><th>PASS variants</th>
    </tr>
  </thead>
  <tbody>{rows}</tbody>
</table>
</div>

<!-- 4. Variant track (below table) -->
{track_section}

<div class="ft">
  Generated by PacBio HiFi Variant Pipeline &mdash; {generated} &mdash; {n} samples
</div>

<script>
// Set thead sticky offset to sit exactly below topbar + legbar.
// window.onload ensures all elements are fully rendered before measuring.
window.onload = function() {{
  var tb  = document.getElementById('topbar');
  var lb  = document.getElementById('legbar');
  var offset = (tb ? tb.offsetHeight : 53) + (lb ? lb.offsetHeight : 32);
  document.querySelector('#t thead').style.top = offset + 'px';
}};

// Row toggle — only fires on the .sr row itself, not the expanded detail
document.querySelectorAll('.sr').forEach(function(row) {{
  row.addEventListener('click', function(e) {{
    if (e.target.tagName === 'PRE' || e.target.tagName === 'CODE') return;
    var detail = document.getElementById(row.getAttribute('data-detail'));
    if (detail) {{
      detail.style.display = detail.style.display === 'none' ? 'table-row' : 'none';
    }}
  }});
}});
</script>
</body>
</html>
"""

    with open(output_path, "w") as f:
        f.write(html)


# --------------------------------------------------------
# Main
# --------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a multi-sample QC summary HTML report."
    )
    parser.add_argument("--reports", nargs="+", required=True,
                        help="QC text report files or glob patterns.")
    parser.add_argument("--output",  required=True,
                        help="Output HTML file path.")
    args = parser.parse_args()

    report_paths = []
    for pattern in args.reports:
        expanded = glob.glob(pattern, recursive=True)
        if expanded:
            report_paths.extend(expanded)
        elif Path(pattern).exists():
            report_paths.append(pattern)
        else:
            print(f"Warning: No files matched: {pattern}")

    if not report_paths:
        raise FileNotFoundError("No QC report files found.")

    report_paths = sorted(set(report_paths))
    print(f"Found {len(report_paths)} report(s).")

    samples = load_reports(report_paths)
    if not samples:
        raise ValueError("No reports could be parsed.")

    print(f"Parsed {len(samples)} sample(s): {[s['sample'] for s in samples]}")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    generate_html(samples, args.output)
    print(f"Report written to: {args.output}")


if __name__ == "__main__":
    main()
