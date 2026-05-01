#!/usr/bin/env python3
"""
05_detect_coverage_anomalies.py

Identifies low-coverage and high-coverage anomalies across all scaffolds
using smoothed coverage profiles (1000 x 10 kb windows = 10 Mb resolution),
consistent with the trend line in coverage plots (scripts 03, 04).

Criteria:
  - High = smoothed depth > scaffold mean + 5x
  - Low  = smoothed depth < scaffold mean - 5x
  - Minimum contiguous region length = 300,000 bp

Corresponds to Additional file 2, section "Read mapping and genome
coverage analysis", step 6 (anomaly detection).

Input:
  eup_cov.regions.bed : mosdepth v0.3.11 [Pedersen & Quinlan, 2018]
                        10 kb windowed coverage (BED4)

Output:
  output_coverage_anomalies/coverage_anomalies_smoothed.tsv

Author: Kopp K, Pristimantis euphronides genome project
"""

import pandas as pd
import numpy as np
import os
from scipy.ndimage import uniform_filter1d

# === Parameters ===
BED_FILE = "/data/GrenadaFrog144/coverage/eup_cov.regions.bed"
OUTDIR = "/data/GrenadaFrog144/coverage/output_coverage_anomalies"
OUTFILE = os.path.join(OUTDIR, "coverage_anomalies_smoothed.tsv")

SMOOTH_WINDOW = 1000    # 1000 x 10 kb = 10 Mb
THRESHOLD_MARGIN = 5    # scaffold mean +/- 5x
MIN_LEN_BP = 300_000    # minimum region length

os.makedirs(OUTDIR, exist_ok=True)

# === Load coverage data ===
df = pd.read_csv(BED_FILE, sep="\t", header=None,
                 names=["scaffold", "start", "end", "depth"])
df["mid"] = (df["start"] + df["end"]) // 2

# === Smooth and scan per scaffold ===
anomalies = []

for scaffold, group in df.groupby("scaffold"):
    group = group.reset_index(drop=True)
    smoothed = uniform_filter1d(group["depth"].values, size=SMOOTH_WINDOW)
    mean_depth = group["depth"].mean()

    lower = mean_depth - THRESHOLD_MARGIN
    upper = mean_depth + THRESHOLD_MARGIN

    group["smoothed"] = smoothed
    group["status"] = np.where(
        smoothed < lower, "low",
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
                if length >= MIN_LEN_BP:
                    anomalies.append({
                        "scaffold": current["scaffold"],
                        "start": current["start"],
                        "end": current["end"],
                        "length_bp": length,
                        "type": current["type"],
                        "mean_smoothed_depth": round(
                            np.mean(current["depths"]), 1),
                        "scaffold_mean": round(mean_depth, 1),
                        "threshold_used": round(
                            lower if current["type"] == "low"
                            else upper, 1)
                    })
                current = None
    # Handle final region at scaffold end
    if current:
        length = current["end"] - current["start"]
        if length >= MIN_LEN_BP:
            anomalies.append({
                "scaffold": current["scaffold"],
                "start": current["start"],
                "end": current["end"],
                "length_bp": length,
                "type": current["type"],
                "mean_smoothed_depth": round(
                    np.mean(current["depths"]), 1),
                "scaffold_mean": round(mean_depth, 1),
                "threshold_used": round(
                    lower if current["type"] == "low"
                    else upper, 1)
            })

# === Save output ===
out_df = pd.DataFrame(anomalies)
out_df = out_df.sort_values(["scaffold", "start"])
out_df.to_csv(OUTFILE, sep="\t", index=False)

# === Print summary ===
print(f"Anomalies detected: {len(out_df)}")
print(f"  Low-coverage:  {len(out_df[out_df['type']=='low'])}")
print(f"  High-coverage: {len(out_df[out_df['type']=='high'])}")
print(f"\nOutput: {OUTFILE}")
