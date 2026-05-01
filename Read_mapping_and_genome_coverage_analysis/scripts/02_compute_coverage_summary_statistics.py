#!/usr/bin/env python3
"""
02_compute_coverage_summary_statistics.py

Computes coverage summary statistics from mosdepth 10 kb windowed output
at two scopes:
  (A) Full assembly (all scaffolds and contigs)
  (B) 13 longest scaffolds only (scaffold_1 through scaffold_13)

For each scope, reports:
  - Number of 10 kb windows
  - Total bases covered
  - Mean coverage (size-weighted)
  - Percentage of 10 kb windows with depth >= 20x
  - Percentage of 10 kb windows with depth < 10x

Additionally reports per-scaffold statistics for scaffolds 1-13:
  - Per-scaffold mean coverage
  - Per-scaffold window counts at >= 20x and < 10x
  - Size-weighted mean across scaffolds 1-13
  - Mean-of-per-scaffold-means (equal weight per scaffold)

Corresponds to Additional file 2, section "Read mapping and genome
coverage analysis", steps 2b-2d and 3.

Input:
  eup_cov.regions.bed : mosdepth v0.3.11 [Pedersen & Quinlan, 2018]
                        output (--by 10000), 4-column BED

Output:
  coverage_summary_statistics.tsv : tab-separated summary table

Author: Kopp K, Pristimantis euphronides genome project
"""

import sys
import os
from collections import defaultdict

# === Configuration ===
BED_FILE = "/data/GrenadaFrog144/coverage/eup_cov.regions.bed"
OUTPUT_FILE = "/data/GrenadaFrog144/coverage/output_coverage_summary/coverage_summary_statistics.tsv"
SCAFFOLDS_TOP13 = [f"scaffold_{i}" for i in range(1, 14)]

# === Verify input ===
if not os.path.exists(BED_FILE):
    sys.exit(f"ERROR: Input file not found: {BED_FILE}")

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

# === Parse BED file ===
scaffold_data = defaultdict(list)
all_depths = []

with open(BED_FILE) as f:
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) < 4:
            continue
        scaffold = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        depth = float(fields[3])
        scaffold_data[scaffold].append((start, end, depth))
        all_depths.append((scaffold, depth, end - start))

print(f"Loaded {len(all_depths)} windows from {len(scaffold_data)} sequences")
print()

# === Helper function ===
def compute_stats(depths_with_sizes, label):
    """Compute summary statistics for a set of (depth, window_size) tuples."""
    n = len(depths_with_sizes)
    if n == 0:
        return None

    total_bases = sum(ws for _, ws in depths_with_sizes)
    mean_depth = sum(d * ws for d, ws in depths_with_sizes) / total_bases
    n_ge20 = sum(1 for d, _ in depths_with_sizes if d >= 20)
    n_lt10 = sum(1 for d, _ in depths_with_sizes if d < 10)
    pct_ge20 = 100.0 * n_ge20 / n
    pct_lt10 = 100.0 * n_lt10 / n

    return {
        "label": label,
        "n_windows": n,
        "total_bases": total_bases,
        "mean_depth": mean_depth,
        "n_ge20": n_ge20,
        "pct_ge20": pct_ge20,
        "n_lt10": n_lt10,
        "pct_lt10": pct_lt10,
    }


# === Scope A: Full assembly (all entries) ===
all_ds = [(d, ws) for _, d, ws in all_depths]
stats_all = compute_stats(all_ds, "Full assembly (all scaffolds + contigs)")

# === Scope B: Scaffolds 1-13 only ===
top13_ds = [(d, ws) for s, d, ws in all_depths if s in SCAFFOLDS_TOP13]
stats_top13 = compute_stats(top13_ds, "Scaffolds 1-13 only")

# === Per-scaffold statistics for scaffolds 1-13 ===
per_scaffold = []
for s in SCAFFOLDS_TOP13:
    windows = scaffold_data.get(s, [])
    if not windows:
        continue
    ds = [(d, e - st) for st, e, d in windows]
    n = len(ds)
    total_bp = sum(ws for _, ws in ds)
    mean_d = sum(d * ws for d, ws in ds) / total_bp
    n_ge20 = sum(1 for d, _ in ds if d >= 20)
    n_lt10 = sum(1 for d, _ in ds if d < 10)
    per_scaffold.append({
        "scaffold": s,
        "n_windows": n,
        "total_bp": total_bp,
        "mean_depth": mean_d,
        "n_ge20": n_ge20,
        "pct_ge20": 100.0 * n_ge20 / n,
        "n_lt10": n_lt10,
        "pct_lt10": 100.0 * n_lt10 / n,
    })

mean_of_means = sum(ps["mean_depth"] for ps in per_scaffold) / len(per_scaffold)
size_weighted_top13 = stats_top13["mean_depth"]

# === Print report ===
print("=" * 80)
print("COVERAGE SUMMARY STATISTICS")
print("=" * 80)

for stats in [stats_all, stats_top13]:
    print(f"\n--- {stats['label']} ---")
    print(f"  Windows:          {stats['n_windows']:,}")
    print(f"  Total bases:      {stats['total_bases']:,} ({stats['total_bases']/1e9:.3f} Gb)")
    print(f"  Mean depth:       {stats['mean_depth']:.2f}x")
    print(f"  Windows >= 20x:   {stats['n_ge20']:,} / {stats['n_windows']:,} = {stats['pct_ge20']:.2f}%")
    print(f"  Windows < 10x:    {stats['n_lt10']:,} / {stats['n_windows']:,} = {stats['pct_lt10']:.2f}%")

print(f"\n--- Per-scaffold statistics (scaffolds 1-13) ---")
print(f"{'Scaffold':<14} {'Windows':>8} {'Mean depth':>12} {'>=20x (%)':>10} {'<10x (%)':>10}")
print("-" * 60)
for ps in per_scaffold:
    print(f"{ps['scaffold']:<14} {ps['n_windows']:>8,} {ps['mean_depth']:>11.2f}x {ps['pct_ge20']:>9.2f}% {ps['pct_lt10']:>9.2f}%")
print("-" * 60)
print(f"{'Size-weighted mean':<24} {size_weighted_top13:>11.2f}x")
print(f"{'Mean-of-scaffold-means':<24} {mean_of_means:>11.2f}x")
print(f"{'Min scaffold mean':<24} {min(ps['mean_depth'] for ps in per_scaffold):>11.2f}x  ({min(per_scaffold, key=lambda x: x['mean_depth'])['scaffold']})")
print(f"{'Max scaffold mean':<24} {max(ps['mean_depth'] for ps in per_scaffold):>11.2f}x  ({max(per_scaffold, key=lambda x: x['mean_depth'])['scaffold']})")

# === Write TSV output ===
with open(OUTPUT_FILE, "w") as out:
    out.write("# Coverage summary statistics\n")
    out.write("# Input: eup_cov.regions.bed (mosdepth v0.3.11, --by 10000)\n")
    out.write("# Window size: 10 kb non-overlapping\n")
    out.write("#\n")
    out.write("scope\tn_windows\ttotal_bases\tmean_depth\tn_ge20x\tpct_ge20x\tn_lt10x\tpct_lt10x\n")
    for stats in [stats_all, stats_top13]:
        out.write(f"{stats['label']}\t{stats['n_windows']}\t{stats['total_bases']}\t"
                  f"{stats['mean_depth']:.2f}\t{stats['n_ge20']}\t{stats['pct_ge20']:.2f}\t"
                  f"{stats['n_lt10']}\t{stats['pct_lt10']:.2f}\n")

    out.write("#\n")
    out.write("# Per-scaffold statistics (scaffolds 1-13)\n")
    out.write("scaffold\tn_windows\ttotal_bp\tmean_depth\tn_ge20x\tpct_ge20x\tn_lt10x\tpct_lt10x\n")
    for ps in per_scaffold:
        out.write(f"{ps['scaffold']}\t{ps['n_windows']}\t{ps['total_bp']}\t"
                  f"{ps['mean_depth']:.2f}\t{ps['n_ge20']}\t{ps['pct_ge20']:.2f}\t"
                  f"{ps['n_lt10']}\t{ps['pct_lt10']:.2f}\n")

    out.write("#\n")
    out.write("# Derived statistics for scaffolds 1-13\n")
    out.write(f"# Size-weighted mean (all bins): {size_weighted_top13:.2f}x\n")
    out.write(f"# Mean-of-per-scaffold-means: {mean_of_means:.2f}x\n")
    out.write(f"# Min scaffold mean: {min(ps['mean_depth'] for ps in per_scaffold):.2f}x "
              f"({min(per_scaffold, key=lambda x: x['mean_depth'])['scaffold']})\n")
    out.write(f"# Max scaffold mean: {max(ps['mean_depth'] for ps in per_scaffold):.2f}x "
              f"({max(per_scaffold, key=lambda x: x['mean_depth'])['scaffold']})\n")

print(f"\nOutput written to: {OUTPUT_FILE}")
