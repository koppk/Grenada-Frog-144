#!/usr/bin/env python3
"""
mann_whitney_identity.py

Exact one-tailed Mann-Whitney U test comparing haplotype identity
between Z-candidate and autosomal regions.

Two comparisons:
  1. Z-candidate vs ALL autosomal regions (genome-wide)
  2. Z-candidate vs same-scaffold flanking autosomal regions only

Reads haplotype_identity.tsv (Layer 4 output), splits values by Type
column (Z_candidate vs Autosomal). Flanking regions are identified
by matching chromosome number between Z and autosomal entries.

Usage:
    python3 mann_whitney_identity.py <haplotype_identity.tsv> <output_dir>

Output:
    <output_dir>/mann_whitney_identity.tsv

Author: Kopp K., Pristimantis euphronides genome project
"""

import sys
import os
from itertools import combinations

infile = sys.argv[1]
outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(infile)
outfile = os.path.join(outdir, "mann_whitney_identity.tsv")

z_vals = []
z_chrs = []       # chromosome numbers of Z regions
a_vals = []
a_regions = []    # (chr_number, identity) for autosomal

with open(infile) as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3 or parts[2] == "NA":
            continue
        region = parts[0]
        identity = float(parts[2])

        # Extract chromosome number from region name (e.g. Chr2_... -> 2)
        chr_num = ""
        if region.startswith("Chr"):
            num_part = region[3:].split("_")[0]
            chr_num = num_part

        if parts[1] == "Z_candidate":
            z_vals.append(identity)
            z_chrs.append(chr_num)
        else:
            a_vals.append(identity)
            a_regions.append((chr_num, identity))

n_z = len(z_vals)
n_a = len(a_vals)
n_total = n_z + n_a

if n_z == 0 or n_a == 0:
    print("  Cannot test: need both Z-candidate and autosomal values.")
    sys.exit(0)

# Identify flanking autosomal regions (same chromosome as any Z region)
z_chr_set = set(z_chrs)
flank_vals = [v for c, v in a_regions if c in z_chr_set]
other_vals = [v for c, v in a_regions if c not in z_chr_set]

z_mean = sum(z_vals) / n_z
a_mean = sum(a_vals) / n_a
f_mean = sum(flank_vals) / len(flank_vals) if flank_vals else 0


def exact_mw_test(group1, group2):
    """Exact one-tailed Mann-Whitney U (H1: group1 < group2).

    Returns U, U_max, p_exact.
    """
    n1 = len(group1)
    n2 = len(group2)

    U_obs = 0
    for g1v in group1:
        for g2v in group2:
            if g1v < g2v:
                U_obs += 1
            elif g1v == g2v:
                U_obs += 0.5

    U_max = n1 * n2

    all_vals = group1 + group2
    n_total = n1 + n2
    count_ge = 0
    total = 0
    for combo in combinations(range(n_total), n1):
        total += 1
        grp1 = [all_vals[i] for i in combo]
        grp2 = [all_vals[i] for i in range(n_total) if i not in combo]
        u = 0
        for gv in grp1:
            for av in grp2:
                if gv < av:
                    u += 1
                elif gv == av:
                    u += 0.5
        if u >= U_obs:
            count_ge += 1

    p = count_ge / total
    return U_obs, U_max, count_ge, total, p


def sig_label(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


# ── Test 1: Z vs all autosomal ────────────────────────────────
U1, U1_max, cnt1, tot1, p1 = exact_mw_test(z_vals, a_vals)
sig1 = sig_label(p1)

# ── Test 2: Z vs same-scaffold flanks ─────────────────────────
if len(flank_vals) >= 1:
    U2, U2_max, cnt2, tot2, p2 = exact_mw_test(z_vals, flank_vals)
    sig2 = sig_label(p2)
else:
    U2, U2_max, cnt2, tot2, p2, sig2 = 0, 0, 0, 0, 1.0, "ns"

# ── Write output file ─────────────────────────────────────────
with open(outfile, "w") as fh:
    fh.write("comparison\tmetric\tgroup1\tgroup1_mean\tgroup1_n\t"
             "group2\tgroup2_mean\tgroup2_n\t"
             "U\tU_max\tp_exact\tsignificance\n")
    fh.write("Z_vs_all_autosomal\thaplotype_identity\t"
             "Z_candidate\t%.6f\t%d\t"
             "Autosomal\t%.6f\t%d\t"
             "%g\t%d\t%.4e\t%s\n" % (
                 z_mean, n_z, a_mean, n_a,
                 U1, U1_max, p1, sig1))
    if flank_vals:
        fh.write("Z_vs_same_scaffold_flanks\thaplotype_identity\t"
                 "Z_candidate\t%.6f\t%d\t"
                 "Flanking_autosomal\t%.6f\t%d\t"
                 "%g\t%d\t%.4e\t%s\n" % (
                     z_mean, n_z, f_mean, len(flank_vals),
                     U2, U2_max, p2, sig2))

# ── Print to stdout ───────────────────────────────────────────
print("=== Mann-Whitney U test (Layer 4) ===")
print("  H0: Z-candidate and autosomal regions have the same")
print("      haplotype identity distribution.")
print("  H1: Z-candidate regions have lower haplotype identity")
print("      than autosomal regions (one-tailed).")
print("")

print("  --- Comparison 1: Z vs all autosomal ---")
print("  Z-candidate:  mean = %.4f  (n = %d)" % (z_mean, n_z))
print("  Autosomal:    mean = %.4f  (n = %d)" % (a_mean, n_a))
print("  Difference:   %.4f" % (a_mean - z_mean))
print("  U = %g / %d  (observed / maximum)" % (U1, U1_max))
print("  Exact one-tailed p = %d / %d = %.4e  %s" % (cnt1, tot1, p1, sig1))
print("")

if flank_vals:
    print("  --- Comparison 2: Z vs same-scaffold flanks ---")
    print("  Z-candidate:       mean = %.4f  (n = %d)" % (z_mean, n_z))
    print("  Flanking autosomal: mean = %.4f  (n = %d)" % (f_mean, len(flank_vals)))
    print("  Difference:        %.4f" % (f_mean - z_mean))
    print("  U = %g / %d  (observed / maximum)" % (U2, U2_max))
    print("  Exact one-tailed p = %d / %d = %.4e  %s" % (cnt2, tot2, p2, sig2))
    print("")

print("  Output: %s" % outfile)
print("")
