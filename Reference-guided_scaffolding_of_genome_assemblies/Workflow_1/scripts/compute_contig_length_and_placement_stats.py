#!/usr/bin/env python3
# Author: Kopp K, Pristimantis euphronides genome project
#
# Compute contig placement and length statistics for the manuscript.
#
# Two outputs:
#   (1) Placement summary (tier counts, total length, mean length, percent of
#       ungapped scaffolded assembly) for the three groups:
#         contigs in scaffolds 1-13
#         contigs in scaffolds 14-31
#         unplaced contigs
#       This feeds the first prose paragraph of the Reference-guided
#       scaffolding section.
#   (2) Contig length distribution statistics for the two placed groups
#       (scaffolds 1-13 vs scaffolds 14-31), feeding Table SR6 / Figure SR4
#       in Additional file 3.
#
# Input files (read from Workflow_1/output/, two-column <contig_name> <length_bp>):
#   scaffolds_1-13.len
#   scaffolds_14-31.len
#   unplaced_contigs.len
#
# Output files:
#   Workflow_1/output/contig_placement_summary.tsv
#   Workflow_1/output/contig_length_and_placement_stats.tsv
#   Workflow_1/figures/contig_length_and_placement_distributions.png
#   Workflow_1/figures/contig_length_and_placement_distributions.tiff
#   Workflow_1/figures/contig_length_and_placement_distributions.pdf
#
# Figures sized for BMC Genomics full-page width (170 mm) at 600 dpi.
# Boxplot rendered black and white, showing only the elements named in the
# Figure SR4 caption: medians, interquartile ranges, whiskers, outliers.
#
# Usage:  python3 compute_contig_length_and_placement_stats.py

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR    = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "output"))
FIG_DIR    = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "figures"))

INPUTS = {
    "scaffolds_1-13":   os.path.join(OUT_DIR, "scaffolds_1-13.len"),
    "scaffolds_14-31":  os.path.join(OUT_DIR, "scaffolds_14-31.len"),
    "unplaced":         os.path.join(OUT_DIR, "unplaced_contigs.len"),
}
SUMMARY_OUT = os.path.join(OUT_DIR, "contig_placement_summary.tsv")
TSV_OUT     = os.path.join(OUT_DIR, "contig_length_and_placement_stats.tsv")
FIG_BASE    = os.path.join(FIG_DIR, "contig_length_and_placement_distributions")

os.makedirs(FIG_DIR, exist_ok=True)

def read_lengths(path):
    if not os.path.isfile(path):
        print(f"ERROR: input file not found: {path}", file=sys.stderr)
        sys.exit(1)
    lengths = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                lengths.append(int(parts[1]))
    if not lengths:
        print(f"ERROR: no contigs read from {path}", file=sys.stderr)
        sys.exit(2)
    return np.array(lengths)

g1 = read_lengths(INPUTS["scaffolds_1-13"])
g2 = read_lengths(INPUTS["scaffolds_14-31"])
g3 = read_lengths(INPUTS["unplaced"])

# Totals
n1, sum1 = len(g1), int(g1.sum())
n2, sum2 = len(g2), int(g2.sum())
n3, sum3 = len(g3), int(g3.sum())

n_placed   = n1 + n2
sum_placed = sum1 + sum2

n_total   = n_placed + n3
sum_total = sum_placed + sum3   # = total ungapped scaffolded assembly length

def pct(x, ref):
    return 100.0 * x / ref if ref > 0 else float("nan")

# Output 1: placement summary TSV
with open(SUMMARY_OUT, "w") as f:
    f.write("Group\tContigs\tTotal length (bp)\tMean length (bp)\t"
            "Percent of ungapped scaffolded assembly\n")
    f.write(f"Total input contigs\t{n_total}\t{sum_total}\t{int(sum_total/n_total)}\t100.000\n")
    f.write(f"All placed contigs\t{n_placed}\t{sum_placed}\t{int(sum_placed/n_placed)}\t{pct(sum_placed, sum_total):.3f}\n")
    f.write(f"Contigs in scaffolds 1-13\t{n1}\t{sum1}\t{int(sum1/n1)}\t{pct(sum1, sum_total):.3f}\n")
    f.write(f"Contigs in scaffolds 14-31\t{n2}\t{sum2}\t{int(sum2/n2)}\t{pct(sum2, sum_total):.3f}\n")
    f.write(f"Unplaced contigs\t{n3}\t{sum3}\t{int(sum3/n3)}\t{pct(sum3, sum_total):.3f}\n")

print(f"Wrote summary: {SUMMARY_OUT}", file=sys.stderr)

# Output 2: length distribution stats for Table SR6 (scaffolds 1-13 vs 14-31 only)
def stats(lengths):
    return {
        "n":               len(lengths),
        "Mean length":     int(np.mean(lengths)),
        "Median (Q2)":     int(np.median(lengths)),
        "Max length":      int(np.max(lengths)),
        "25th percentile": int(np.percentile(lengths, 25)),
        "75th percentile": int(np.percentile(lengths, 75)),
    }

s1 = stats(g1)
s2 = stats(g2)

col1 = f"Contigs in scaffolds 1-13 (n = {s1['n']:,})"
col2 = f"Contigs in scaffolds 14-31 (n = {s2['n']:,})"

metrics = ["Mean length", "Median (Q2)", "Max length", "25th percentile", "75th percentile"]

with open(TSV_OUT, "w") as f:
    f.write(f"Metric (bp)\t{col1}\t{col2}\n")
    for m in metrics:
        f.write(f"{m}\t{s1[m]}\t{s2[m]}\n")

print(f"Wrote stats: {TSV_OUT}", file=sys.stderr)

# Output 3: boxplot figure (scaffolds 1-13 vs 14-31 only)
# BMC Genomics figure requirements:
#   Single-column width: 85 mm  = 3.346 in
#   Full-page width:    170 mm  = 6.693 in
#   Raster resolution: 300 dpi color/greyscale, 600-1200 dpi line art
#   Vector preferred (PDF), TIFF with LZW otherwise
#   Sans-serif >= 8 pt, line weight > 0.25 pt
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,   # embed TrueType (editable text in PDF)
    "ps.fonttype": 42,
})

fig, ax = plt.subplots(figsize=(6.693, 4.5))
ax.boxplot(
    [g1, g2],
    tick_labels=["Scaffolds 1-13", "Scaffolds 14-31"],
    widths=0.6,
    showfliers=True,
    medianprops=dict(color="black", linewidth=1.0),
    boxprops=dict(color="black", linewidth=0.75),
    whiskerprops=dict(color="black", linewidth=0.75),
    capprops=dict(color="black", linewidth=0.75),
    flierprops=dict(marker="o", markersize=3, markerfacecolor="none",
                    markeredgecolor="black", linewidth=0.5),
)

ax.set_yscale("log")
ax.set_ylabel("Contig length (bp)")

plt.tight_layout()

png_path  = f"{FIG_BASE}.png"
tiff_path = f"{FIG_BASE}.tiff"
pdf_path  = f"{FIG_BASE}.pdf"

plt.savefig(png_path,  dpi=600, bbox_inches="tight")
print(f"Wrote figure: {png_path}", file=sys.stderr)

plt.savefig(tiff_path, dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print(f"Wrote figure: {tiff_path}", file=sys.stderr)

plt.savefig(pdf_path,  bbox_inches="tight")
print(f"Wrote figure: {pdf_path}", file=sys.stderr)

plt.close(fig)
