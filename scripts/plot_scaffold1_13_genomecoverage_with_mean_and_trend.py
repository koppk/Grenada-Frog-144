#!/usr/bin/env python3
"""
plot_scaffold1_13_genomecoverage_with_mean_and_trend.py

Description:
------------
Generates a genome-wide coverage plot for ONT reads across scaffolds 1–13
using absolute genome positions. Adds:
- Thin black line: lightly smoothed base-level coverage
- Thick blue line: heavily smoothed trend line (broad signal)
- Blue dashed horizontal line: mean coverage across scaffolds 1–13
- Blue dashed vertical lines: scaffold boundaries
- Scaffold tick labels below x-axis (as secondary x-axis), rotated for readability

Input:
------
- eup_cov.regions.bed: mosdepth output (10 kb bins)
- ragtag.scaffold.renamed.agp: AGP file defining scaffold structure

Output:
-------
- genomewide_scaffolds1_13_with_trend_and_mean.png

Author: Katharina Kopp
Date: 2025-05-12
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d

# === Input files ===
cov_file = "eup_cov.regions.bed"
agp_file = "ragtag.scaffold.renamed.agp"

# === Define target scaffolds (first 13)
scaffolds_to_plot = [f"scaffold_{i}" for i in range(1, 14)]

# === Load coverage and AGP data ===
cov = pd.read_csv(cov_file, sep="\t", header=None,
                  names=["scaffold", "start", "end", "depth"])
agp = pd.read_csv(agp_file, sep="\t", comment="#", header=None,
                  names=["object", "obj_start", "obj_end", "part_number", "component_type",
                         "component_id_or_gap_length", "component_beg_or_gap_type",
                         "component_end_or_linkage", "orientation_or_evidence"])

# === Filter to scaffolds 1–13 only
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
cov["abs_mid"] = cov.apply(lambda row: scaffold_offsets[row["scaffold"]] + row["mid"], axis=1)
cov["smoothed_black"] = uniform_filter1d(cov["depth"].values, size=2)        # Fine smoothing
cov["smoothed_blue"] = uniform_filter1d(cov["depth"].values, size=1000)      # Broad trend smoothing
mean_depth = cov["smoothed_black"].mean()

# === Plot ===
fig, ax = plt.subplots(figsize=(14, 5))

# Black line: base coverage
ax.plot(cov["abs_mid"], cov["smoothed_black"], color='black', linewidth=0.1)

# Blue line: smooth coverage trend
ax.plot(cov["abs_mid"], cov["smoothed_blue"], color='blue', linewidth=0.8)

# Dashed horizontal mean line
ax.axhline(mean_depth, color='blue', linestyle='--', linewidth=0.5)

# Mean annotation slightly outside right edge
ax.annotate(f"Mean: {mean_depth:.1f}×",
            xy=(1.01, mean_depth),
            xycoords=("axes fraction", "data"),
            fontsize=8, color='blue',
            ha='left', va='bottom', clip_on=False)

# === Draw scaffold boundary lines ===
for s in scaffolds_to_plot:
    ax.axvline(scaffold_offsets[s], color='blue', linestyle='--', linewidth=0.5)

# === Secondary x-axis for scaffold labels ===
tick_positions = [scaffold_offsets[s] for s in scaffolds_to_plot]
tick_labels = [f"Scaffold {i+1}" for i in range(len(scaffolds_to_plot))]

ax2 = ax.secondary_xaxis('bottom')
ax2.set_xticks(tick_positions)
ax2.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=8, color='blue')
ax2.tick_params(axis='x', pad=30)

# === Final formatting ===
ax.set_xlabel("Genome position (bp)")
ax.set_ylabel("Coverage (×)")
ax.set_ylim(0, 60)
ax.set_title("ONT Read Coverage Across Scaffolds 1–13")
plt.tight_layout()
plt.savefig("genomewide_scaffolds1_13_with_trend_and_mean.png", dpi=300)

