#!/usr/bin/env python3
"""
detect_smoothed_coverage_anomalies.py

Description:
------------
Identifies low-coverage and high-coverage anomalies across all scaffolds using
smoothed coverage profiles (1000 × 10 kb windows = 10 Mb resolution), consistent
with the trend line in coverage plots.

Criteria:
---------
- Spike = smoothed depth > mean + 5×
- Dip   = smoothed depth < mean - 5×
- Minimum region length = 300,000 bp

Input:
------
- eup_cov.regions.bed : 10 kb coverage windows from mosdepth

Output:
-------
- coverage_anomalies_smoothed.tsv

Author: Katharina Kopp
Date: 2025-05-12
"""

import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d

# === Parameters
bed_file = "eup_cov.regions.bed"
smooth_window = 1000  # 1000 × 10 kb = 10 Mb
threshold_margin = 5  # mean ± 5×
min_len_bp = 300_000  # minimum region length

# === Load coverage data
df = pd.read_csv(bed_file, sep="\t", header=None,
                 names=["scaffold", "start", "end", "depth"])
df["mid"] = (df["start"] + df["end"]) // 2

# === Smooth and scan per scaffold
anomalies = []
for scaffold, group in df.groupby("scaffold"):
    group = group.reset_index(drop=True)
    smoothed = uniform_filter1d(group["depth"].values, size=smooth_window)
    mean_depth = group["depth"].mean()
    lower = mean_depth - threshold_margin
    upper = mean_depth + threshold_margin
    group["smoothed"] = smoothed
    group["status"] = np.where(smoothed < lower, "low",
                         np.where(smoothed > upper, "high", "normal"))

    # Extract contiguous regions by status
    current = None
    for _, row in group.iterrows():
        if row["status"] != "normal":
            if current is None:
                current = {
                    "scaffold": scaffold,
                    "start": row["start"],
                    "end": row["end"],
                    "depths": [row["smoothed"]],
                    "type": row["status"]
                }
            else:
                current["end"] = row["end"]
                current["depths"].append(row["smoothed"])
        else:
            if current:
                length = current["end"] - current["start"]
                if length >= min_len_bp:
                    anomalies.append({
                        "scaffold": current["scaffold"],
                        "start": current["start"],
                        "end": current["end"],
                        "length_bp": length,
                        "type": current["type"],
                        "mean_smoothed_depth": round(np.mean(current["depths"]), 1),
                        "scaffold_mean": round(mean_depth, 1),
                        "threshold_used": round(lower if current["type"] == "low" else upper, 1)
                    })
                current = None
    if current:
        length = current["end"] - current["start"]
        if length >= min_len_bp:
            anomalies.append({
                "scaffold": current["scaffold"],
                "start": current["start"],
                "end": current["end"],
                "length_bp": length,
                "type": current["type"],
                "mean_smoothed_depth": round(np.mean(current["depths"]), 1),
                "scaffold_mean": round(mean_depth, 1),
                "threshold_used": round(lower if current["type"] == "low" else upper, 1)
            })

# === Save output
out_df = pd.DataFrame(anomalies)
out_df = out_df.sort_values(["scaffold", "start"])
out_df.to_csv("coverage_anomalies_smoothed.tsv", sep="\t", index=False)

print("Done. Output saved to 'coverage_anomalies_smoothed.tsv'")

