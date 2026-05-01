#!/usr/bin/env python3
"""
parse_scaffold_repeats.py

Parse WindowMasker intervals and RepeatMasker .out files for scaffolds
1-13 into per-Mb repeat density profiles. Compare Z-candidate regions
vs autosomal regions.

Input:
    scaffolds_fasta     FASTA with scaffold_1..13 (for lengths via .fai)
    wm_intervals        WindowMasker interval output
    rm_out              RepeatMasker .out annotation table
    z_regions_tsv       z_candidate_regions.tsv (Layer 3)
    output_dir          Directory for output files

Output:
    scaffold_repeat_per_mb.tsv      Per-Mb bins with WM and RM masked bp
    scaffold_repeat_summary.tsv     Z-candidate vs autosomal comparison

Usage:
    python3 parse_scaffold_repeats.py <scaffolds.fasta> <wm_intervals.txt> \
        <rm.out> <z_candidate_regions.tsv> <output_dir>

Author: Kopp K., Pristimantis euphronides genome project
"""

import sys
import os
from collections import defaultdict

scaffolds_fa = sys.argv[1]
wm_intervals = sys.argv[2]
rm_out_file = sys.argv[3]
z_regions = sys.argv[4]
outdir = sys.argv[5]

MB = 1000000

# ── Read scaffold lengths from FAI ─────────────────────────────
fai = scaffolds_fa + ".fai"
scaff_lengths = {}
with open(fai) as f:
    for line in f:
        parts = line.strip().split("\t")
        scaff_lengths[parts[0]] = int(parts[1])

scaffolds = sorted(scaff_lengths.keys(),
                   key=lambda s: int(s.split("_")[1]))

print("  Scaffolds: %d" % len(scaffolds))
print("  Total: %.1f Mb" % (sum(scaff_lengths.values()) / 1e6))

# ── Read Z-candidate regions ──────────────────────────────────
z_regions_list = []
with open(z_regions) as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            z_regions_list.append((parts[0], int(parts[1]), int(parts[2])))

def bin_type(scaff, mb):
    for zs, zstart, zend in z_regions_list:
        if scaff == zs and zstart <= mb <= zend:
            return "Z_candidate"
    return "Autosomal"


# ── Initialize per-Mb bins ─────────────────────────────────────
rm_classes = ["SINE", "LINE", "LTR", "DNA", "Satellite",
              "Simple_repeat", "Other_repeat"]

bins = {}
for scaff in scaffolds:
    n_mb = scaff_lengths[scaff] // MB
    for mb in range(n_mb + 1):
        d = {"wm_masked_bp": 0, "rm_total_bp": 0}
        for c in rm_classes:
            d[c] = 0
        bins[(scaff, mb)] = d


def bin_size(scaff, mb):
    remaining = scaff_lengths[scaff] - mb * MB
    return min(MB, remaining)


def distribute_interval(scaff, start_bp, end_bp, field, amount=None):
    """Add overlap of [start_bp, end_bp) to each Mb bin it spans.
    If amount is None, adds actual overlap in bp.
    """
    mb_start = start_bp // MB
    mb_end = (end_bp - 1) // MB if end_bp > start_bp else mb_start

    for mb in range(mb_start, mb_end + 1):
        bin_lo = mb * MB
        bin_hi = (mb + 1) * MB
        overlap_start = max(start_bp, bin_lo)
        overlap_end = min(end_bp, bin_hi)
        overlap = overlap_end - overlap_start
        if overlap > 0:
            key = (scaff, mb)
            if key in bins:
                bins[key][field] += overlap


# ── Parse WindowMasker intervals ──────────────────────────────
print("  Parsing WindowMasker intervals...")
current_scaff = None
wm_count = 0
with open(wm_intervals) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current_scaff = line[1:].strip()
            continue
        if " - " in line:
            parts = line.split(" - ")
            start = int(parts[0].strip())
            end = int(parts[1].strip()) + 1  # WM intervals are inclusive
            wm_count += 1

            if current_scaff in scaff_lengths:
                distribute_interval(current_scaff, start, end,
                                    "wm_masked_bp")

print("    %d intervals processed" % wm_count)


# ── Parse RepeatMasker .out ────────────────────────────────────
print("  Parsing RepeatMasker .out file...")

def classify_rm(class_family):
    cl = class_family.split("/")[0]
    if cl == "SINE":
        return "SINE"
    elif cl == "LINE":
        return "LINE"
    elif cl in ("LTR", "Retroposon"):
        return "LTR"
    elif cl == "DNA":
        return "DNA"
    elif cl.startswith("Satellite"):
        return "Satellite"
    elif cl in ("Simple_repeat", "Low_complexity"):
        return "Simple_repeat"
    else:
        return "Other_repeat"

rm_count = 0

if os.path.exists(rm_out_file):
    with open(rm_out_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                _ = int(parts[0])
            except ValueError:
                continue

            scaff = parts[4]
            start = int(parts[5]) - 1  # RM is 1-based, convert to 0-based
            end = int(parts[6])        # now [start, end) in 0-based
            class_family = parts[10]
            rm_count += 1

            if scaff not in scaff_lengths:
                continue

            cat = classify_rm(class_family)
            distribute_interval(scaff, start, end, "rm_total_bp")
            distribute_interval(scaff, start, end, cat)

    print("    %d annotations processed" % rm_count)
else:
    print("    WARNING: RepeatMasker .out file not found: %s" % rm_out_file)
    print("    (RepeatMasker columns will be zero — rerun after RM finishes)")


# ── Write per-Mb output ────────────────────────────────────────
out_per_mb = os.path.join(outdir, "scaffold_repeat_per_mb.tsv")
print("  Writing per-Mb profiles: %s" % out_per_mb)

with open(out_per_mb, "w") as f:
    header = ["scaffold", "Mb", "bin_size_bp", "type",
              "wm_masked_bp", "wm_masked_frac",
              "rm_total_bp", "rm_masked_frac"]
    header += ["rm_%s_bp" % c for c in rm_classes]
    header += ["rm_%s_frac" % c for c in rm_classes]
    f.write("\t".join(header) + "\n")

    for scaff in scaffolds:
        n_mb = scaff_lengths[scaff] // MB
        for mb in range(n_mb + 1):
            key = (scaff, mb)
            b = bins[key]
            bs = bin_size(scaff, mb)
            bt = bin_type(scaff, mb)

            wm_frac = b["wm_masked_bp"] / bs if bs > 0 else 0
            rm_frac = b["rm_total_bp"] / bs if bs > 0 else 0

            row = [scaff, str(mb), str(bs), bt,
                   str(b["wm_masked_bp"]), "%.4f" % wm_frac,
                   str(b["rm_total_bp"]), "%.4f" % rm_frac]

            for c in rm_classes:
                row.append(str(b[c]))
            for c in rm_classes:
                row.append("%.4f" % (b[c] / bs if bs > 0 else 0))

            f.write("\t".join(row) + "\n")


# ── Summary: Z-candidate vs Autosomal ──────────────────────────
out_summary = os.path.join(outdir, "scaffold_repeat_summary.tsv")
print("  Writing summary: %s" % out_summary)

agg = {}
for bt in ["Z_candidate", "Autosomal"]:
    agg[bt] = {"total_bp": 0, "wm_bp": 0, "rm_bp": 0, "n_bins": 0}
    for c in rm_classes:
        agg[bt][c] = 0

for scaff in scaffolds:
    n_mb = scaff_lengths[scaff] // MB
    for mb in range(n_mb + 1):
        key = (scaff, mb)
        b = bins[key]
        bs = bin_size(scaff, mb)
        bt = bin_type(scaff, mb)

        agg[bt]["total_bp"] += bs
        agg[bt]["wm_bp"] += b["wm_masked_bp"]
        agg[bt]["rm_bp"] += b["rm_total_bp"]
        agg[bt]["n_bins"] += 1
        for c in rm_classes:
            agg[bt][c] += b[c]

with open(out_summary, "w") as f:
    header = ["type", "n_bins", "total_Mb",
              "wm_masked_frac", "rm_masked_frac"]
    header += ["rm_%s_frac" % c for c in rm_classes]
    f.write("\t".join(header) + "\n")

    for bt in ["Z_candidate", "Autosomal"]:
        a = agg[bt]
        total = a["total_bp"]
        row = [bt, str(a["n_bins"]),
               "%.1f" % (total / 1e6),
               "%.4f" % (a["wm_bp"] / total if total > 0 else 0),
               "%.4f" % (a["rm_bp"] / total if total > 0 else 0)]
        for c in rm_classes:
            row.append("%.4f" % (a[c] / total if total > 0 else 0))
        f.write("\t".join(row) + "\n")

# Print summary to stdout
print("")
print("  === Repeat density summary ===")
print("")
for bt in ["Z_candidate", "Autosomal"]:
    a = agg[bt]
    total = a["total_bp"]
    print("  %s  (%d bins, %.1f Mb):" % (bt, a["n_bins"], total / 1e6))
    if total > 0:
        print("    WindowMasker:  %.2f%% masked" % (100 * a["wm_bp"] / total))
        print("    RepeatMasker:  %.2f%% masked" % (100 * a["rm_bp"] / total))
        for c in rm_classes:
            frac = 100 * a[c] / total
            if frac > 0.01:
                print("      %-15s %.2f%%" % (c, frac))
    print("")
