#!/usr/bin/env python3
"""
04_plot_individual_scaffolds.py

Generates individual per-scaffold coverage plots for ONT reads aligned
to scaffolds 1-13. Each figure shows:
  - Lightly smoothed black coverage trace (2-bin uniform filter)
  - Broad blue trend line (1000-bin uniform filter)
  - Blue dashed horizontal line: per-scaffold mean coverage
  - Mean coverage label outside the plot frame

Corresponds to Additional file 2, section "Read mapping and genome
coverage analysis", step 4 (smoothing) and Additional file 9.

Input:
  eup_cov.regions.bed         : mosdepth 10 kb windowed coverage (BED4)
  ragtag.scaffold.renamed.agp : RagTag v2.1.0 [Alonge, 2022] AGP file

Output:
  output_scaffold_plots/scaffold_1.png to scaffold_13.png

Author: Kopp K, Pristimantis euphronides genome project
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import numpy as np
import os
from scipy.ndimage import uniform_filter1d

# === Input files ===
COV_FILE = "/data/GrenadaFrog144/coverage/eup_cov.regions.bed"
AGP_FILE = "/data/GrenadaFrog144/coverage/ragtag.scaffold.renamed.agp"
OUTDIR = "/data/GrenadaFrog144/coverage/output_scaffold_plots"

os.makedirs(OUTDIR, exist_ok=True)

# === Target scaffolds ===
scaffolds = [f"scaffold_{i}" for i in range(1, 14)]

# === Load coverage and AGP ===
cov = pd.read_csv(COV_FILE, sep="\t", header=None,
                  names=["scaffold", "start", "end", "depth"])
agp = pd.read_csv(AGP_FILE, sep="\t", comment="#", header=None,
                  names=["object", "obj_start", "obj_end", "part_number",
                         "component_type", "component_id_or_gap_length",
                         "component_beg_or_gap_type",
                         "component_end_or_linkage",
                         "orientation_or_evidence"])

# === Filter to scaffolds 1-13 ===
cov = cov[cov["scaffold"].isin(scaffolds)].copy()
agp = agp[agp["object"].isin(scaffolds)].copy()

# === Compute scaffold lengths and absolute genome offsets ===
scaffold_lengths = agp.groupby("object")["obj_end"].max()
scaffold_offsets = {}
offset = 0
for s in scaffolds:
    scaffold_offsets[s] = offset
    offset += scaffold_lengths[s]

# === Process and plot each scaffold ===
for i, s in enumerate(scaffolds):
    cov_s = cov[cov["scaffold"] == s].copy()
    cov_s["mid"] = (cov_s["start"] + cov_s["end"]) // 2
    cov_s["abs_mid"] = cov_s["mid"] + scaffold_offsets[s]

    cov_s["smoothed"] = uniform_filter1d(cov_s["depth"].values, size=2)
    smoothed_blue = uniform_filter1d(cov_s["depth"].values, size=1000)
    avg_depth = cov_s["smoothed"].mean()

    fig, ax = plt.subplots(figsize=(14, 4))

    ax.plot(cov_s["abs_mid"], cov_s["smoothed"],
            color='black', linewidth=0.1)
    ax.plot(cov_s["abs_mid"], smoothed_blue,
            color='blue', linewidth=0.8)
    ax.axhline(avg_depth, color='blue', linestyle='--', linewidth=0.5)

    ax.annotate(f"Mean: {avg_depth:.1f}\u00d7",
                xy=(1.01, avg_depth),
                xycoords=("axes fraction", "data"),
                fontsize=8, color='blue',
                ha='left', va='bottom', clip_on=False)

    ax.set_xlabel("Genome position (Gb)")
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, pos: f'{x * 1e-9:.2f}'))
    ax.set_ylabel("Coverage (\u00d7)")
    ax.set_ylim(0, 60)
    ax.set_title(f"ONT Read Coverage \u2013 Scaffold {i+1}")

    plt.tight_layout()
    outfile = os.path.join(OUTDIR, f"{s}.png")
    fig.savefig(outfile, dpi=300)
    plt.close()
    print(f"  {outfile}")

print(f"\nDone. {len(scaffolds)} scaffold plots written to {OUTDIR}")
