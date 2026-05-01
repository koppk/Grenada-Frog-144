#!/usr/bin/env python3
"""
identify_z_candidate_regions.py

Layer 1: Statistical identification of scaffolds with significant
enrichment of half-coverage 1 Mb bins.

Aggregates mosdepth 10 kb windowed coverage into 1 Mb bins per scaffold.
Computes genome-wide background rate of half-coverage bins (< 0.625 x
median scaffold mean). Tests each scaffold for enrichment using a
one-sided binomial test with Bonferroni correction.

Output:
    scaffold_screening_summary.tsv   Per-scaffold test results
    scaffold_screening.log           Full log

Usage:
    python3 identify_z_candidate_regions.py <eup_cov.regions.bed.gz> <output_dir>
Author: Kopp K, Pristimantis euphronides genome project
"""

import gzip
import math
import sys
import os
from collections import defaultdict
from datetime import datetime

HEMI_STRICT = 0.625
ALPHA = 0.05
N_SCAFFOLDS = 13


def binomial_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p). One-sided enrichment p-value."""
    if p <= 0:
        return 0.0 if k > 0 else 1.0
    if p >= 1:
        return 1.0 if k <= n else 0.0
    if k <= 0:
        return 1.0
    log_p = math.log(p)
    log_1mp = math.log(1.0 - p)
    cumulative = 0.0
    for i in range(k):
        log_pmf = (math.lgamma(n + 1) - math.lgamma(i + 1)
                   - math.lgamma(n - i + 1)
                   + i * log_p + (n - i) * log_1mp)
        cumulative += math.exp(log_pmf)
        if cumulative >= 1.0:
            return 0.0
    return max(0.0, 1.0 - cumulative)


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 identify_z_candidate_regions.py "
              "<eup_cov.regions.bed.gz> <output_dir>")
        sys.exit(1)

    bedgz = sys.argv[1]
    outdir = sys.argv[2]

    if not os.path.isfile(bedgz):
        print("ERROR: File not found: " + bedgz, file=sys.stderr)
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    logpath = os.path.join(outdir, "scaffold_screening.log")
    logfh = open(logpath, "w")

    def log(msg=""):
        print(msg)
        logfh.write(msg + "\n")

    log("=== identify_z_candidate_regions.py (Layer 1) ===")
    log("Start: " + str(datetime.now()))
    log("Input: " + bedgz)
    log("Output: " + outdir)
    log()
    log("Parameters:")
    log("  HEMI_STRICT = " + str(HEMI_STRICT) + "  (half-cov threshold, fraction of reference)")
    log("  ALPHA       = " + str(ALPHA) + "  (significance level, Bonferroni-corrected)")
    log()

    # Step 1: Read 10 kb windows, aggregate to 1 Mb bins
    log("[Step 1] Reading 10 kb windows and aggregating to 1 Mb bins...")

    scaff_data = defaultdict(lambda: defaultdict(list))
    valid_scaffolds = set("scaffold_" + str(i) for i in range(1, N_SCAFFOLDS + 1))

    opener = gzip.open if bedgz.endswith(".gz") else open
    with opener(bedgz, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            scaff = parts[0]
            if scaff not in valid_scaffolds:
                continue
            start = int(parts[1])
            cov = float(parts[3])
            mb_bin = start // 1000000
            scaff_data[scaff][mb_bin].append(cov)

    per_mb = {}
    for scaff in sorted(scaff_data, key=lambda s: int(s.split("_")[1])):
        bins = []
        for mb in sorted(scaff_data[scaff]):
            vals = scaff_data[scaff][mb]
            bins.append((mb, sum(vals) / len(vals), len(vals)))
        per_mb[scaff] = bins

    total_bins = sum(len(v) for v in per_mb.values())
    log("  Total 1 Mb bins across %d scaffolds: %d" % (len(per_mb), total_bins))
    log()

    # Step 2: Genome-wide reference
    log("[Step 2] Computing genome-wide reference coverage...")

    scaff_means = {}
    for scaff, bins in per_mb.items():
        mean_cov = sum(c for _, c, _ in bins) / len(bins)
        scaff_means[scaff] = (mean_cov, len(bins))

    sorted_means = sorted(scaff_means[s][0] for s in scaff_means)
    n = len(sorted_means)
    if n % 2 == 1:
        reference = sorted_means[n // 2]
    else:
        reference = (sorted_means[n // 2 - 1] + sorted_means[n // 2]) / 2.0

    strict_cutoff = reference * HEMI_STRICT

    log("  Per-scaffold mean coverages:")
    for scaff in sorted(scaff_means, key=lambda s: int(s.split("_")[1])):
        mc, nb = scaff_means[scaff]
        log("    %-14s %.2fx  (%d Mb)" % (scaff, mc, nb))
    log()
    log("  Genome-wide reference (median of scaffold means): %.2fx" % reference)
    log("  Half-coverage cutoff (%.3f x %.2f): < %.2fx" % (HEMI_STRICT, reference, strict_cutoff))
    log()

    # Step 3: Binomial test per scaffold
    log("[Step 3] Testing each scaffold for enrichment of half-coverage bins...")
    log()
    log("  Method: One-sided binomial test per scaffold.")
    log("    H0: Number of half-coverage bins is consistent with")
    log("        the genome-wide background rate.")
    log("    H1: This scaffold has significantly more half-coverage")
    log("        bins than expected by chance.")
    log("    Correction: Bonferroni (%d tests)." % len(per_mb))
    log()

    scaff_half = {}
    scaff_ntotal = {}
    for scaff, bins in per_mb.items():
        total = len(bins)
        half_count = sum(1 for _, c, _ in bins if c < strict_cutoff)
        scaff_half[scaff] = half_count
        scaff_ntotal[scaff] = total

    total_half = sum(scaff_half.values())
    total_all = sum(scaff_ntotal.values())
    bg_rate = total_half / total_all if total_all > 0 else 0

    log("  Genome-wide half-coverage bins: %d / %d (background rate p = %.6f)" % (total_half, total_all, bg_rate))

    n_tests = len(per_mb)
    log("  Bonferroni correction: %d tests, alpha_corrected = %.4f" % (n_tests, ALPHA / n_tests))
    log()

    screening_path = os.path.join(outdir, "scaffold_screening_summary.tsv")
    screening_rows = []

    for scaff in sorted(per_mb, key=lambda s: int(s.split("_")[1])):
        n_bins = scaff_ntotal[scaff]
        k = scaff_half[scaff]
        expected = n_bins * bg_rate
        frac = k / n_bins if n_bins > 0 else 0
        pval_raw = binomial_sf(k, n_bins, bg_rate)
        pval_corrected = min(1.0, pval_raw * n_tests)
        flagged = pval_corrected < ALPHA
        screening_rows.append((scaff, n_bins, k, frac, expected, pval_raw, pval_corrected, flagged))

    with open(screening_path, "w") as fh:
        fh.write("scaffold\tn_bins\thalf_cov_bins\tfraction\texpected\tp_value_raw\tp_value_bonferroni\tflagged\n")
        for row in screening_rows:
            scaff, nb, k, frac, exp, praw, pcorr, flag = row
            fh.write("%s\t%d\t%d\t%.4f\t%.2f\t%.2e\t%.2e\t%s\n" % (
                scaff, nb, k, frac, exp, praw, pcorr, "YES" if flag else "no"))

    log("  Results:")
    log("    %-14s %4s  %4s  %5s  %6s  %10s  %10s" % ("Scaffold", "Bins", "Obs", "Exp", "Frac", "p(raw)", "p(Bonf)"))
    log("    %s %s  %s  %s  %s  %s  %s" % ("-"*14, "-"*4, "-"*4, "-"*5, "-"*6, "-"*10, "-"*10))
    for row in screening_rows:
        scaff, nb, k, frac, exp, praw, pcorr, flag = row
        marker = " ***" if flag else ""
        log("    %-14s %4d  %4d  %5.1f  %5.1f%%  %10.2e  %10.2e%s" % (
            scaff, nb, k, exp, frac*100, praw, pcorr, marker))
    log()

    flagged_scaffolds = [s for s, _, _, _, _, _, _, f in screening_rows if f]
    log("  Flagged scaffolds (%d): %s" % (len(flagged_scaffolds), " ".join(flagged_scaffolds)))
    log("  Output: " + screening_path)
    log()
    log("Done: " + str(datetime.now()))
    logfh.close()


if __name__ == "__main__":
    main()
