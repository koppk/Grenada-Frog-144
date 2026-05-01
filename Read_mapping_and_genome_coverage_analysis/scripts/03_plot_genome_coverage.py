#!/usr/bin/env python3
"""
03_plot_genome_coverage.py

Generates a genome-wide coverage plot for ONT reads across scaffolds 1-13
using absolute genome positions. Displays:
  - Thin black line: lightly smoothed coverage (2-bin uniform filter)
  - Thick blue line: broad trend line (1000-bin uniform filter)
  - Blue dashed horizontal line: size-weighted mean coverage
  - Blue dashed vertical lines: scaffold boundaries

Corresponds to Additional file 2, section "Read mapping and genome
coverage analysis", steps 4 (smoothing) and 5 (genome-wide plot).

Input:
  eup_cov.regions.bed         : mosdepth 10 kb windowed coverage (BED4)
  ragtag.scaffold.renamed.agp : RagTag v2.1.0 [Alonge, 2022] AGP file

Output:
  output_coverage_plots/genomewide_scaffolds1_13_with_trend_and_mean.png

Author: Kopp K, Pristimantis euphronides genome project
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy.ndimage import uniform_filter1d

# === Input files ===
COV_FILE = "/data/GrenadaFrog144/coverage/eup_cov.regions.bed"
AGP_FILE = "/data/GrenadaFrog144/coverage/ragtag.scaffold.renamed.agp"
OUTDIR = "/data/GrenadaFrog144/coverage/output_coverage_plots"
OUTFILE = os.path.join(OUTDIR, "genomewide_scaffolds1_13_with_trend_and_mean.png")

os.makedirs(OUTDIR, exist_ok=True)

# === Target scaffolds ===
scaffolds_to_plot = [f"scaffold_{i}" for i in range(1, 14)]

# === Load coverage and AGP data ===
cov = pd.read_csv(COV_FILE, sep="\t", header=None,
                  names=["scaffold", "start", "end", "depth"])
agp = pd.read_csv(AGP_FILE, sep="\t", comment="#", header=None,
                  names=["object", "obj_start", "obj_end", "part_number",
                         "component_type", "component_id_or_gap_length",
                         "component_beg_or_gap_type",
                         "component_end_or_linkage",
                         "orientation_or_evidence"])

# === Filter to scaffolds 1-13 ===
cov = cov[cov["scaffold"].isin(scaffolds_to_plot)].copy()
agp = agp[agp["object"].isin(scaffolds_to_plot)].copy()

# === Compute scaffold lengths and cumulative genome offsets ===
scaffold_lengths = agp.groupby("object")["obj_end"].max()
scaffold_offsets = {}
offset = 0
for s in scaffolds_to_plot:
    scaffold_offsets[s] = offset
    offset += scaffold_lengths[s]

# === Compute absolute midpoints and smoothed values ===
cov["mid"] = (cov["start"] + cov["end"]) // 2
cov["abs_mid"] = cov.apply(
    lambda row: scaffold_offsets[row["scaffold"]] + row["mid"], axis=1)
cov["smoothed_black"] = uniform_filter1d(cov["depth"].values, size=2)
cov["smoothed_blue"] = uniform_filter1d(cov["depth"].values, size=1000)

mean_depth = cov.groupby("scaffold")["smoothed_black"].mean().mean()

# === Plot ===
fig, ax = plt.subplots(figsize=(14, 5))

ax.plot(cov["abs_mid"], cov["smoothed_black"],
        color='black', linewidth=0.1)
ax.plot(cov["abs_mid"], cov["smoothed_blue"],
        color='blue', linewidth=0.8)
ax.axhline(mean_depth, color='blue', linestyle='--', linewidth=0.5)

ax.annotate(f"Mean: {mean_depth:.1f}\u00d7",
            xy=(1.01, mean_depth),
            xycoords=("axes fraction", "data"),
            fontsize=8, color='blue',
            ha='left', va='bottom', clip_on=False)

for s in scaffolds_to_plot:
    ax.axvline(scaffold_offsets[s], color='blue',
               linestyle='--', linewidth=0.5)

# === Secondary x-axis for scaffold labels ===
tick_positions = [scaffold_offsets[s] for s in scaffolds_to_plot]
tick_labels = [f"Scaffold {i+1}" for i in range(len(scaffolds_to_plot))]
ax2 = ax.secondary_xaxis('bottom')
ax2.set_xticks(tick_positions)
ax2.set_xticklabels(tick_labels, rotation=45, ha='right',
                    fontsize=8, color='blue')
ax2.tick_params(axis='x', pad=30)

# === Formatting ===
ax.set_xlabel("Genome position (bp)")
ax.set_ylabel("Coverage (\u00d7)")
ax.set_ylim(0, 60)
ax.set_title("ONT Read Coverage Across Scaffolds 1\u201313")

plt.tight_layout()
plt.savefig(OUTFILE, dpi=300)
print(f"Output: {OUTFILE}")
