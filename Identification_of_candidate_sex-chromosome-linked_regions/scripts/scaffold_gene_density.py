#!/usr/bin/env python3
"""
scaffold_gene_density.py

Layer 6: Per-Mb gene density across scaffolds 1-13.

Reads gene coordinates from HANNO BESTMODELS bedDB, bins genes into
1 Mb windows, compares Z-candidate vs autosomal gene density.

Two comparisons:
  1. Z-candidate vs ALL autosomal bins (genome-wide)
  2. Z-candidate vs same-scaffold flanking autosomal bins only

Gene assignment: each gene is assigned to the Mb bin containing its
midpoint (start + end) / 2.

Input:
    beddb               HANNO BESTMODELS-FINAL.bedDB
    z_regions_tsv       z_candidate_regions.tsv (Layer 3)
    output_dir          Directory for output files

Output:
    scaffold_gene_density_per_mb.tsv    Per-Mb gene counts
    scaffold_gene_density_summary.tsv   Z-candidate vs autosomal comparison
    mann_whitney_gene_density.tsv       Statistical tests

Usage:
    python3 scaffold_gene_density.py <BESTMODELS.bedDB> \
        <z_candidate_regions.tsv> <output_dir>

Author: Kopp K., Pristimantis euphronides genome project
"""

import sys
import os
import math
from collections import defaultdict

beddb = sys.argv[1]
z_regions = sys.argv[2]
outdir = sys.argv[3]

MB = 1000000

# ── Read Z-candidate regions ──────────────────────────────────
z_regions_list = []
z_scaffolds = set()
with open(z_regions) as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            z_regions_list.append((parts[0], int(parts[1]), int(parts[2])))
            z_scaffolds.add(parts[0])


def bin_type(scaff, mb):
    for zs, zstart, zend in z_regions_list:
        if scaff == zs and zstart <= mb <= zend:
            return "Z_candidate"
    return "Autosomal"


# ── Target scaffolds ──────────────────────────────────────────
target_scaffolds = set("scaffold_%d" % i for i in range(1, 14))

# ── Read genes from bedDB ─────────────────────────────────────
print("=== Gene density analysis (Layer 6) ===")
print("  Input: %s" % beddb)
print("")

scaff_max_mb = defaultdict(int)
gene_counts = defaultdict(int)
total_genes = 0
skipped = 0

# Per-gene lengths grouped by region type
z_gene_lengths = []
a_gene_lengths = []
f_gene_lengths = []    # flanking autosomal on same scaffold as Z

with open(beddb) as f:
    header = f.readline().strip().split("\t")
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue

        scaff = parts[0]
        if scaff not in target_scaffolds:
            skipped += 1
            continue

        start = int(parts[1])
        end = int(parts[2])
        gene_len = end - start
        midpoint = (start + end) // 2
        mb = midpoint // MB

        gene_counts[(scaff, mb)] += 1
        if mb > scaff_max_mb[scaff]:
            scaff_max_mb[scaff] = mb
        total_genes += 1

        bt = bin_type(scaff, mb)
        if bt == "Z_candidate":
            z_gene_lengths.append(gene_len)
        else:
            a_gene_lengths.append(gene_len)
            if scaff in z_scaffolds:
                f_gene_lengths.append(gene_len)

print("  Genes on scaffolds 1-13: %d" % total_genes)
print("  Genes on other scaffolds (skipped): %d" % skipped)
print("")

# ── Build sorted scaffold list ─────────────────────────────────
scaffolds = sorted(scaff_max_mb.keys(),
                   key=lambda s: int(s.split("_")[1]))

# ── Write per-Mb output ────────────────────────────────────────
out_per_mb = os.path.join(outdir, "scaffold_gene_density_per_mb.tsv")
print("  Writing per-Mb profiles: %s" % out_per_mb)

z_densities = []
a_densities = []
f_densities = []    # flanking autosomal on same scaffold as Z

with open(out_per_mb, "w") as f:
    f.write("scaffold\tMb\ttype\tgene_count\n")
    for scaff in scaffolds:
        max_mb = scaff_max_mb[scaff]
        for mb in range(max_mb + 1):
            bt = bin_type(scaff, mb)
            count = gene_counts.get((scaff, mb), 0)
            f.write("%s\t%d\t%s\t%d\n" % (scaff, mb, bt, count))

            if bt == "Z_candidate":
                z_densities.append(count)
            else:
                a_densities.append(count)
                if scaff in z_scaffolds:
                    f_densities.append(count)

# ── Summary ────────────────────────────────────────────────────
out_summary = os.path.join(outdir, "scaffold_gene_density_summary.tsv")
print("  Writing summary: %s" % out_summary)

z_total = sum(z_densities)
a_total = sum(a_densities)
f_total = sum(f_densities)
n_z = len(z_densities)
n_a = len(a_densities)
n_f = len(f_densities)
z_mean = z_total / n_z if n_z > 0 else 0
a_mean = a_total / n_a if n_a > 0 else 0
f_mean = f_total / n_f if n_f > 0 else 0

with open(out_summary, "w") as f:
    f.write("type\tn_bins\ttotal_genes\tmean_genes_per_Mb\n")
    f.write("Z_candidate\t%d\t%d\t%.2f\n" % (n_z, z_total, z_mean))
    f.write("Autosomal_all\t%d\t%d\t%.2f\n" % (n_a, a_total, a_mean))
    f.write("Flanking_same_scaffold\t%d\t%d\t%.2f\n" % (n_f, f_total, f_mean))

print("")
print("  === Gene density summary ===")
print("")
print("  Z_candidate:              %d bins, %d genes, %.2f genes/Mb" % (
    n_z, z_total, z_mean))
print("  Autosomal (all):          %d bins, %d genes, %.2f genes/Mb" % (
    n_a, a_total, a_mean))
print("  Flanking (same scaffold): %d bins, %d genes, %.2f genes/Mb" % (
    n_f, f_total, f_mean))
print("")


# ── Mann-Whitney U test ───────────────────────────────────────
def mann_whitney_normal_approx(x, y):
    """One-tailed Mann-Whitney U (H1: x < y, i.e. y stochastically greater)."""
    nx = len(x)
    ny = len(y)

    U = 0
    for xv in x:
        for yv in y:
            if xv < yv:
                U += 1
            elif xv == yv:
                U += 0.5

    mu = nx * ny / 2.0
    sigma = math.sqrt(nx * ny * (nx + ny + 1) / 12.0)

    if sigma == 0:
        return U, nx * ny, 1.0

    z_stat = (U - mu - 0.5) / sigma
    p = 0.5 * math.erfc(z_stat / math.sqrt(2))

    return U, nx * ny, p


def sig_label(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


# Test 1: Z vs all autosomal
U1, U1_max, p1 = mann_whitney_normal_approx(z_densities, a_densities)
sig1 = sig_label(p1)

# Test 2: Z vs same-scaffold flanks
if n_f > 0:
    U2, U2_max, p2 = mann_whitney_normal_approx(z_densities, f_densities)
    sig2 = sig_label(p2)
else:
    U2, U2_max, p2, sig2 = 0, 0, 1.0, "ns"

# Write test output
out_test = os.path.join(outdir, "mann_whitney_gene_density.tsv")
with open(out_test, "w") as f:
    f.write("comparison\tmetric\tgroup1\tgroup1_mean\tgroup1_n\t"
            "group2\tgroup2_mean\tgroup2_n\t"
            "U\tU_max\tp\tsignificance\n")
    f.write("Z_vs_all_autosomal\tgene_density\t"
            "Z_candidate\t%.2f\t%d\t"
            "Autosomal\t%.2f\t%d\t"
            "%.0f\t%d\t%.2e\t%s\n" % (
                z_mean, n_z, a_mean, n_a, U1, U1_max, p1, sig1))
    if n_f > 0:
        f.write("Z_vs_same_scaffold_flanks\tgene_density\t"
                "Z_candidate\t%.2f\t%d\t"
                "Flanking_autosomal\t%.2f\t%d\t"
                "%.0f\t%d\t%.2e\t%s\n" % (
                    z_mean, n_z, f_mean, n_f, U2, U2_max, p2, sig2))

print("  === Mann-Whitney U tests ===")
print("  H0: Z-candidate and autosomal bins have the same")
print("      gene density distribution.")
print("  H1: Z-candidate bins have lower gene density")
print("      than autosomal bins (one-tailed).")
print("")

print("  --- Comparison 1: Z vs all autosomal ---")
print("  Z-candidate:  mean = %.2f genes/Mb  (n = %d)" % (z_mean, n_z))
print("  Autosomal:    mean = %.2f genes/Mb  (n = %d)" % (a_mean, n_a))
print("  U = %.0f / %d" % (U1, U1_max))
print("  p = %.2e  %s" % (p1, sig1))
print("")

if n_f > 0:
    print("  --- Comparison 2: Z vs same-scaffold flanks ---")
    print("  Z-candidate:       mean = %.2f genes/Mb  (n = %d)" % (z_mean, n_z))
    print("  Flanking autosomal: mean = %.2f genes/Mb  (n = %d)" % (f_mean, n_f))
    print("  U = %.0f / %d" % (U2, U2_max))
    print("  p = %.2e  %s" % (p2, sig2))
    print("")

print("  Method: normal approximation with continuity correction")
print("  Output: %s" % out_test)
print("")

# ── Gene body length comparison ────────────────────────────────
# Tests whether Z genes have larger gene bodies (end - start)
# than autosomal genes, which would indicate intronic repeat
# expansion rather than intergenic accumulation.

def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return 0
    if n % 2 == 1:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2


def mann_whitney_greater(x, y):
    """One-tailed Mann-Whitney U (H1: x > y)."""
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

    z_stat = (U - mu - 0.5) / sigma
    p = 0.5 * math.erfc(z_stat / math.sqrt(2))

    return U, nx * ny, p


nz_g = len(z_gene_lengths)
na_g = len(a_gene_lengths)
nf_g = len(f_gene_lengths)

z_gl_mean = sum(z_gene_lengths) / nz_g if nz_g > 0 else 0
a_gl_mean = sum(a_gene_lengths) / na_g if na_g > 0 else 0
f_gl_mean = sum(f_gene_lengths) / nf_g if nf_g > 0 else 0
z_gl_med = median(z_gene_lengths)
a_gl_med = median(a_gene_lengths)
f_gl_med = median(f_gene_lengths)

# Test 1: Z vs all autosomal gene lengths
UL1, UL1_max, pL1 = mann_whitney_greater(z_gene_lengths, a_gene_lengths)
sigL1 = sig_label(pL1)

# Test 2: Z vs same-scaffold flank gene lengths
if nf_g > 0:
    UL2, UL2_max, pL2 = mann_whitney_greater(z_gene_lengths, f_gene_lengths)
    sigL2 = sig_label(pL2)
else:
    UL2, UL2_max, pL2, sigL2 = 0, 0, 1.0, "ns"

# Write gene length results to output files
out_gene_length = os.path.join(outdir, "gene_body_length_summary.tsv")
with open(out_gene_length, "w") as f:
    f.write("type\tn_genes\tmean_length_bp\tmedian_length_bp\n")
    f.write("Z_candidate\t%d\t%.0f\t%.0f\n" % (nz_g, z_gl_mean, z_gl_med))
    f.write("Autosomal_all\t%d\t%.0f\t%.0f\n" % (na_g, a_gl_mean, a_gl_med))
    f.write("Flanking_same_scaffold\t%d\t%.0f\t%.0f\n" % (
        nf_g, f_gl_mean, f_gl_med))

# Append gene length tests to mann_whitney output
with open(out_test, "a") as f:
    f.write("Z_vs_all_autosomal\tgene_body_length\t"
            "Z_candidate\t%.0f\t%d\t"
            "Autosomal\t%.0f\t%d\t"
            "%.0f\t%d\t%.2e\t%s\n" % (
                z_gl_mean, nz_g, a_gl_mean, na_g,
                UL1, UL1_max, pL1, sigL1))
    if nf_g > 0:
        f.write("Z_vs_same_scaffold_flanks\tgene_body_length\t"
                "Z_candidate\t%.0f\t%d\t"
                "Flanking_autosomal\t%.0f\t%d\t"
                "%.0f\t%d\t%.2e\t%s\n" % (
                    z_gl_mean, nz_g, f_gl_mean, nf_g,
                    UL2, UL2_max, pL2, sigL2))

# Print
print("  === Gene body length comparison ===")
print("  H1: Z-candidate genes have larger gene bodies")
print("      than autosomal genes (one-tailed).")
print("")
print("  %-30s  %8s  %8s  %8s" % ("Group", "n", "mean bp", "median bp"))
print("  " + "-" * 60)
print("  %-30s  %8d  %8.0f  %8.0f" % ("Z_candidate", nz_g, z_gl_mean, z_gl_med))
print("  %-30s  %8d  %8.0f  %8.0f" % ("Autosomal (all)", na_g, a_gl_mean, a_gl_med))
print("  %-30s  %8d  %8.0f  %8.0f" % ("Flanking (same scaffold)", nf_g, f_gl_mean, f_gl_med))
print("")

print("  --- Z vs all autosomal ---")
print("  U = %.0f / %d,  p = %.2e  %s" % (UL1, UL1_max, pL1, sigL1))
print("")

if nf_g > 0:
    print("  --- Z vs same-scaffold flanks ---")
    print("  U = %.0f / %d,  p = %.2e  %s" % (UL2, UL2_max, pL2, sigL2))
    print("")

print("  Output: %s" % out_gene_length)
print("")

