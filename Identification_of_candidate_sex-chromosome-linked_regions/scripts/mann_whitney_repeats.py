#!/usr/bin/env python3
"""
mann_whitney_repeats.py

Mann-Whitney U test comparing per-Mb repeat density between
Z-candidate and autosomal bins.

Two comparisons per metric:
  1. Z-candidate vs ALL autosomal bins (genome-wide)
  2. Z-candidate vs same-scaffold flanking autosomal bins only

Uses normal approximation with continuity correction (large n makes
exact enumeration infeasible).

Usage:
    python3 mann_whitney_repeats.py <scaffold_repeat_per_mb.tsv> [output_dir]

Output:
    <output_dir>/mann_whitney_repeats.tsv

Author: Kopp K., Pristimantis euphronides genome project
"""

import sys
import os
import math

infile = sys.argv[1]
outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(infile)
outfile = os.path.join(outdir, "mann_whitney_repeats.tsv")

# ── Read per-Mb data ───────────────────────────────────────────
z_data = {}        # field -> list of values
a_data = {}        # all autosomal
f_data = {}        # flanking autosomal (same scaffold as Z regions)

fields_frac = []
header = None

# First pass: identify which scaffolds contain Z-candidate bins
z_scaffolds = set()
with open(infile) as f:
    header = f.readline().strip().split("\t")
    for line in f:
        parts = line.strip().split("\t")
        row = dict(zip(header, parts))
        if row["type"] == "Z_candidate":
            z_scaffolds.add(row["scaffold"])

# Second pass: read data into groups
with open(infile) as f:
    header = f.readline().strip().split("\t")

    for col in header:
        if col.endswith("_frac"):
            fields_frac.append(col)
            z_data[col] = []
            a_data[col] = []
            f_data[col] = []

    for line in f:
        parts = line.strip().split("\t")
        row = dict(zip(header, parts))
        bt = row["type"]
        scaff = row["scaffold"]

        if bt == "Z_candidate":
            for col in fields_frac:
                z_data[col].append(float(row[col]))
        else:
            for col in fields_frac:
                a_data[col].append(float(row[col]))
            if scaff in z_scaffolds:
                for col in fields_frac:
                    f_data[col].append(float(row[col]))

n_z = len(z_data[fields_frac[0]])
n_a = len(a_data[fields_frac[0]])
n_f = len(f_data[fields_frac[0]])


def mann_whitney_normal_approx(x, y):
    """One-tailed Mann-Whitney U (H1: x > y).

    Returns U statistic and p-value using normal approximation
    with continuity correction.
    """
    nx = len(x)
    ny = len(y)

    U = 0
    for xv in x:
        for yv in y:
            if xv > yv:
                U += 1
            elif xv == yv:
                U += 0.5

    mu = nx * ny / 2.0
    sigma = math.sqrt(nx * ny * (nx + ny + 1) / 12.0)

    if sigma == 0:
        return U, nx * ny, 1.0

    z = (U - mu - 0.5) / sigma
    p = 0.5 * math.erfc(z / math.sqrt(2))

    return U, nx * ny, p


# ── Display names ──────────────────────────────────────────────
display_names = {
    "wm_masked_frac": "WindowMasker total",
    "rm_masked_frac": "RepeatMasker total",
    "rm_SINE_frac": "SINE",
    "rm_LINE_frac": "LINE",
    "rm_LTR_frac": "LTR",
    "rm_DNA_frac": "DNA transposons",
    "rm_Satellite_frac": "Satellite",
    "rm_Simple_repeat_frac": "Simple repeat",
    "rm_Other_repeat_frac": "Other/Unclassified",
}

# ── Run tests ──────────────────────────────────────────────────
results = []  # (comparison, col, name, z_mean, comp_mean, n_comp, U, U_max, p, sig)

for col in fields_frac:
    z_vals = z_data[col]
    a_vals = a_data[col]
    f_vals = f_data[col]
    name = display_names.get(col, col)

    z_mean = sum(z_vals) / len(z_vals)

    # Comparison 1: Z vs all autosomal
    a_mean = sum(a_vals) / len(a_vals)
    U1, U1_max, p1 = mann_whitney_normal_approx(z_vals, a_vals)
    sig1 = "***" if p1 < 0.001 else "**" if p1 < 0.01 else "*" if p1 < 0.05 else "ns"
    results.append(("Z_vs_all_autosomal", col, name,
                     z_mean, a_mean, n_a, U1, U1_max, p1, sig1))

    # Comparison 2: Z vs same-scaffold flanks
    if n_f > 0:
        f_mean = sum(f_vals) / len(f_vals)
        U2, U2_max, p2 = mann_whitney_normal_approx(z_vals, f_vals)
        sig2 = "***" if p2 < 0.001 else "**" if p2 < 0.01 else "*" if p2 < 0.05 else "ns"
        results.append(("Z_vs_same_scaffold_flanks", col, name,
                         z_mean, f_mean, n_f, U2, U2_max, p2, sig2))

# ── Write output file ─────────────────────────────────────────
with open(outfile, "w") as fh:
    fh.write("comparison\tmetric\tdisplay_name\t"
             "Z_candidate_mean\tAutosomal_mean\t"
             "Z_n\tAutosomal_n\tU\tU_max\tp\tsignificance\n")
    for comp, col, name, zm, am, na, U, Um, p, sig in results:
        fh.write("%s\t%s\t%s\t%.6f\t%.6f\t%d\t%d\t%.0f\t%.0f\t%.2e\t%s\n" % (
            comp, col, name, zm, am, n_z, na, U, Um, p, sig))

# ── Print to stdout ───────────────────────────────────────────
print("=== Mann-Whitney U tests: repeat density (Layer 5) ===")
print("")
print("  H0: Z-candidate and autosomal Mb bins have the same")
print("      repeat density distribution.")
print("  H1: Z-candidate bins have higher repeat density")
print("      than autosomal bins (one-tailed).")
print("")

print("  --- Comparison 1: Z vs all autosomal ---")
print("  Z-candidate bins: n = %d" % n_z)
print("  Autosomal bins:   n = %d" % n_a)
print("")
print("  %-25s  %8s  %8s  %10s  %8s" % (
    "Metric", "Z_mean%", "Auto%", "U", "p"))
print("  " + "-" * 70)
for comp, col, name, zm, am, na, U, Um, p, sig in results:
    if comp == "Z_vs_all_autosomal":
        print("  %-25s  %7.2f%%  %7.2f%%  %10.0f  %.2e %s" % (
            name, 100 * zm, 100 * am, U, p, sig))

print("")
print("  --- Comparison 2: Z vs same-scaffold flanks ---")
print("  Z-candidate bins: n = %d" % n_z)
print("  Flanking bins:    n = %d  (scaffolds: %s)" % (
    n_f, ", ".join(sorted(z_scaffolds))))
print("")
print("  %-25s  %8s  %8s  %10s  %8s" % (
    "Metric", "Z_mean%", "Flank%", "U", "p"))
print("  " + "-" * 70)
for comp, col, name, zm, am, na, U, Um, p, sig in results:
    if comp == "Z_vs_same_scaffold_flanks":
        print("  %-25s  %7.2f%%  %7.2f%%  %10.0f  %.2e %s" % (
            name, 100 * zm, 100 * am, U, p, sig))

print("")
print("  Method: normal approximation with continuity correction")
print("  Output: %s" % outfile)
print("")
