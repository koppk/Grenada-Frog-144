#!/usr/bin/env python3
"""
plot_scaffolds_1to13_individually.py

Description:
------------
Generates individual per-scaffold coverage plots for ONT reads aligned to scaffolds 1–13.
Each figure shows:
- Lightly smoothed black coverage trace (window size = 2)
- Broad blue trend line (window size = 1000)
- Blue dashed horizontal line for scaffold mean coverage
- Mean coverage label outside the plot frame
- X-axis = absolute genome position (cumulative offset)
- One PNG image per scaffold

Input:
------
- eup_cov.regions.bed: coverage file (e.g., 10 kb bins from mosdepth)
- ragtag.scaffold.renamed.agp: AGP file defining scaffold positions and lengths

Output:
-------
- scaffold_plots/scaffold_1.png to scaffold_13.png

Author: Katharina Kopp
Date: 2025-05-12
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import numpy as np
import os
from scipy.ndimage import uniform_filter1d

# === Input files ===
cov_file = "eup_cov.regions.bed"
agp_file = "ragtag.scaffold.renamed.agp"

# === Target scaffolds
scaffolds = [f"scaffold_{i}" for i in range(1, 14)]

# === Load coverage and AGP
cov = pd.read_csv(cov_file, sep="\t", header=None,
                  names=["scaffold", "start", "end", "depth"])
agp = pd.read_csv(agp_file, sep="\t", comment="#", header=None,
                  names=["object", "obj_start", "obj_end", "part_number", "component_type",
                         "component_id_or_gap_length", "component_beg_or_gap_type",
                         "component_end_or_linkage", "orientation_or_evidence"])

# === Filter coverage and AGP to scaffolds 1–13
cov = cov[cov["scaffold"].isin(scaffolds)].copy()
agp = agp[agp["object"].isin(scaffolds)].copy()

# === Compute scaffold lengths and absolute genome offsets
scaffold_lengths = agp.groupby("object")["obj_end"].max()
scaffold_offsets = {}
offset = 0
for s in scaffolds:
    scaffold_offsets[s] = offset
    offset += scaffold_lengths[s]

# === Create output directory
os.makedirs("scaffold_plots", exist_ok=True)

# === Process and plot each scaffold
for i, s in enumerate(scaffolds):
    cov_s = cov[cov["scaffold"] == s].copy()
    cov_s["mid"] = (cov_s["start"] + cov_s["end"]) // 2
    cov_s["abs_mid"] = cov_s["mid"] + scaffold_offsets[s]

    # Smoothed coverage lines
    cov_s["smoothed"] = uniform_filter1d(cov_s["depth"].values, size=2)
    smoothed_blue = uniform_filter1d(cov_s["depth"].values, size=1000)

    # Compute mean coverage for scaffold
    avg_depth = cov_s["smoothed"].mean()

    # === Plot
    fig, ax = plt.subplots(figsize=(14, 4))

    # Black line: lightly smoothed coverage
    ax.plot(cov_s["abs_mid"], cov_s["smoothed"], color='black', linewidth=0.1)

    # Blue line: broad smoothed trend
    ax.plot(cov_s["abs_mid"], smoothed_blue, color='blue', linewidth=0.8)

    # Blue dashed horizontal mean line
    ax.axhline(avg_depth, color='blue', linestyle='--', linewidth=0.5)

    # Annotate mean outside the plot area
    ax.annotate(f"Mean: {avg_depth:.1f}×",
                xy=(1.01, avg_depth),
                xycoords=("axes fraction", "data"),
                fontsize=8, color='blue',
                ha='left', va='bottom', clip_on=False)

    # Axes formatting
    ax.set_xlabel("Genome position (Gb)")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: f'{x * 1e-9:.2f}'))
    ax.set_ylabel("Coverage (×)")
    ax.set_ylim(0, 60)
    ax.set_title(f"ONT Read Coverage – Scaffold {i+1}")
    plt.tight_layout()
    fig.savefig(f"scaffold_plots/{s}.png", dpi=300)
    plt.close()
