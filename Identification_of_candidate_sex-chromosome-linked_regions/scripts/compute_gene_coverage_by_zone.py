#!/usr/bin/env python3
# Author: Kopp K, Pristimantis euphronides genome project
#
# compute_gene_coverage_by_zone.py
# Compute gene coverage classification by compositional zone within
# Z-candidate regions of scaffold_2 and scaffold_8 (Table SR28).
#
# Reads Z-region segment boundaries from the pixel-based segment
# TSVs produced by plot_z_regions_banding.py (Tables SR26/SR27),
# then cross-references with placed_genes.tsv to tabulate coverage
# classes per segment.
#
# Output contains one row per Z-region segment (10 rows total).
# Table SR28 in the manuscript groups adjacent scaffold_2 segments
# by rendered class for presentation (3 grouped zones on scaffold_2,
# 2 segments on scaffold_8 = 5 rows).
#
# Usage:
#   python3 compute_gene_coverage_by_zone.py [BASE_DIR]
#
# BASE_DIR defaults to the server working directory.

import csv
import os
import sys
from collections import defaultdict

BASE = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"

PLACED_GENES = os.path.join(BASE, "gametolog_discovery_hanno7", "placed_genes.tsv")
FIG_DIR = os.path.join(BASE, "gametolog_discovery_hanno7", "figures")
OUT = os.path.join(FIG_DIR, "gene_coverage_by_zone.tsv")

SEGMENT_FILES = [
    os.path.join(FIG_DIR, "z_banding_scaffold_2_segments.tsv"),
    os.path.join(FIG_DIR, "z_banding_scaffold_8_segments.tsv"),
]

COV_CLASSES = ["Hemi_0.5x", "Auto_1.0x", "High_Coverage/Repeat",
               "Low_Coverage/Other"]

# ---- Read Z-region segments from TSVs ----
zones = []
for seg_file in SEGMENT_FILES:
    print(f"Reading {os.path.basename(seg_file)} ...")
    with open(seg_file) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["in_Z_region"] != "Z":
                continue
            zones.append({
                "scaffold": r["scaffold"],
                "start": float(r["seg_start_Mb"]) * 1e6,
                "end": float(r["seg_end_Mb"]) * 1e6,
                "start_mb": float(r["seg_start_Mb"]),
                "end_mb": float(r["seg_end_Mb"]),
                "cls": r["rendered_class"],
            })

print(f"  {len(zones)} Z-region segments")
for z in zones:
    span = z["end_mb"] - z["start_mb"]
    print(f"    {z['scaffold']} {z['start_mb']:.2f}-{z['end_mb']:.2f} Mb "
          f"({span:.2f} Mb): {z['cls']}")

# ---- Assign genes to segments ----
zone_counts = [defaultdict(int) for _ in zones]

print(f"\nReading {os.path.basename(PLACED_GENES)} ...")
n_total = 0
n_assigned = 0
with open(PLACED_GENES) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        n_total += 1
        if not r.get("z_flag", "").startswith("Z-candidate"):
            continue
        scf = r["scaffold"]
        mid = (int(r["scf_start"]) + int(r["scf_end"])) / 2
        cov = r["cov_class"]
        for i, z in enumerate(zones):
            if scf == z["scaffold"] and z["start"] <= mid < z["end"]:
                zone_counts[i]["total"] += 1
                zone_counts[i][cov] += 1
                n_assigned += 1
                break

print(f"  {n_total} genes read, {n_assigned} assigned to Z-region segments")

# ---- Write output ----
print(f"\nWriting {os.path.basename(OUT)} ...")
with open(OUT, "w") as f:
    hdr = ["scaffold", "seg_start_Mb", "seg_end_Mb", "rendered_class",
           "span_Mb", "total_genes"]
    for cc in COV_CLASSES:
        safe = cc.replace("/", "_")
        hdr.extend([safe, safe + "_pct"])
    f.write("\t".join(hdr) + "\n")

    for i, z in enumerate(zones):
        c = zone_counts[i]
        total = c["total"]
        span = z["end_mb"] - z["start_mb"]
        row = [z["scaffold"], f"{z['start_mb']:.2f}", f"{z['end_mb']:.2f}",
               z["cls"], f"{span:.2f}", str(total)]
        for cc in COV_CLASSES:
            n = c[cc]
            pct = f"{n / total * 100:.1f}" if total > 0 else "0.0"
            row.extend([str(n), pct])
        f.write("\t".join(row) + "\n")
        print(f"  {z['scaffold']} {z['start_mb']:.1f}-{z['end_mb']:.1f}: "
              f"{total} genes (Hemi={c['Hemi_0.5x']}, Auto={c['Auto_1.0x']})")

print("\nDone.")
