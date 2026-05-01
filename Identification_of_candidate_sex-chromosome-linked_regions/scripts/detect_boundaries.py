#!/usr/bin/env python3
"""
detect_boundaries.py

Layer 3: Detect half-coverage block boundaries on a single scaffold.

Method:
  1. Cap outlier bins at 2x median coverage
  2. Optimal binary segmentation (2-segment vs 3-segment, BIC)
     to identify the core half-coverage block
  3. Extend boundaries outward into adjacent bins falling more
     than 2 SD below the flanking-segment mean

Reports: start_Mb, end_Mb (inclusive), length = end - start + 1.

Usage:
    python3 detect_boundaries.py <scaffold> <per_mb_tsv> <output_tsv>

Appends one line per detected region to the output TSV.
Author: Kopp K, Pristimantis euphronides genome project
"""

import sys
import math

scaff = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

# Read per-Mb coverage
covs = []
with open(infile) as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            covs.append((int(parts[0]), float(parts[1])))

n = len(covs)
if n < 15:
    print("  %s: too few bins (%d), skipping" % (scaff, n))
    sys.exit(0)

mbs = [m for m, _ in covs]
vals_raw = [c for _, c in covs]

# Cap outliers at 2x median before segmentation
srt = sorted(vals_raw)
median_cov = srt[n // 2]
cap = 2.0 * median_cov
vals = [min(v, cap) for v in vals_raw]
ncapped = sum(1 for v in vals_raw if v > cap)
print("  %s: %d Mb bins, median=%.1fx, capped %d outlier bins at %.1fx" % (
    scaff, n, median_cov, ncapped, cap))


def seg_cost(v, s, e):
    """Sum of squared deviations from segment mean."""
    seg = v[s:e]
    if not seg:
        return 0.0
    m = sum(seg) / len(seg)
    return sum((x - m) ** 2 for x in seg)


def seg_mean(v, s, e):
    seg = v[s:e]
    return sum(seg) / len(seg) if seg else 0.0


# 2-segment model: one break point
best2_cost = float("inf")
best2_bp = 0
for bp in range(5, n - 5):
    c = seg_cost(vals, 0, bp) + seg_cost(vals, bp, n)
    if c < best2_cost:
        best2_cost = c
        best2_bp = bp

bic2 = n * math.log(best2_cost / n + 1e-10) + 3 * math.log(n)

# 3-segment model: two break points
best3_cost = float("inf")
best3_bp1 = 0
best3_bp2 = 0
for b1 in range(5, n - 10):
    cl = seg_cost(vals, 0, b1)
    for b2 in range(b1 + 5, n - 5):
        c = cl + seg_cost(vals, b1, b2) + seg_cost(vals, b2, n)
        if c < best3_cost:
            best3_cost = c
            best3_bp1 = b1
            best3_bp2 = b2

bic3 = n * math.log(best3_cost / n + 1e-10) + 5 * math.log(n)

print("    2-segment: break at Mb %d  (BIC=%.1f)" % (mbs[best2_bp], bic2))
print("    3-segment: breaks at Mb %d, %d  (BIC=%.1f)" % (
    mbs[best3_bp1], mbs[best3_bp2], bic3))

# Determine core block indices (half-open: [idx_start, idx_end))
if bic3 < bic2:
    model = "3-segment"
    means = [
        seg_mean(vals, 0, best3_bp1),
        seg_mean(vals, best3_bp1, best3_bp2),
        seg_mean(vals, best3_bp2, n),
    ]
    mi = means.index(min(means))
    if mi == 0:
        idx_start, idx_end = 0, best3_bp1
    elif mi == 1:
        idx_start, idx_end = best3_bp1, best3_bp2
    else:
        idx_start, idx_end = best3_bp2, n
    flank_mean = sum(m for i, m in enumerate(means) if i != mi) / 2
    region_mean = means[mi]
    # Collect flanking-segment values for SD calculation
    if mi == 0:
        flank_vals = vals[best3_bp1:best3_bp2] + vals[best3_bp2:]
    elif mi == 1:
        flank_vals = vals[0:best3_bp1] + vals[best3_bp2:]
    else:
        flank_vals = vals[0:best3_bp1] + vals[best3_bp1:best3_bp2]
else:
    model = "2-segment"
    mean_left = seg_mean(vals, 0, best2_bp)
    mean_right = seg_mean(vals, best2_bp, n)
    if mean_left < mean_right:
        idx_start, idx_end = 0, best2_bp
        region_mean, flank_mean = mean_left, mean_right
        flank_vals = vals[best2_bp:]
    else:
        idx_start, idx_end = best2_bp, n
        region_mean, flank_mean = mean_right, mean_left
        flank_vals = vals[0:best2_bp]

print("    Core block: Mb %d - %d  (indices %d:%d)" % (
    mbs[idx_start], mbs[idx_end - 1], idx_start, idx_end))

# Extend boundaries using flanking-segment statistics.
# A bin is included if it falls more than 2 SD below the
# flanking-segment mean (capped values), i.e. it is unlikely
# to be autosomal (p < 0.023, one-tailed).
flank_sd = math.sqrt(sum((x - flank_mean) ** 2 for x in flank_vals) / len(flank_vals))
threshold = flank_mean - 2.0 * flank_sd

print("    Flanking: mean=%.2fx, SD=%.2f, threshold (mean-2SD)=%.2fx" % (
    flank_mean, flank_sd, threshold))

while idx_start > 0 and vals_raw[idx_start - 1] < threshold:
    idx_start -= 1

while idx_end < n and vals_raw[idx_end] < threshold:
    idx_end += 1

# Report: start and end are inclusive Mb labels
s = mbs[idx_start]
e = mbs[idx_end - 1]
length = e - s + 1
ratio = region_mean / flank_mean if flank_mean > 0 else 0

with open(outfile, "a") as fh:
    fh.write("%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.3f\t%s\n" % (
        scaff, s, e, length, region_mean, flank_mean, ratio, model))

print("    Final block: Mb %d - %d  (%d Mb)" % (s, e, length))
print("    Region mean: %.1fx, flanking mean: %.1fx, ratio: %.2f" % (
    region_mean, flank_mean, ratio))
print("")
