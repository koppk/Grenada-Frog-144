#!/usr/bin/env python3
"""
plot_tier1a_gametolog_synteny.py
================================
Reads tier 1a gametolog data from the gametolog_discovery_hanno7 pipeline
output files and generates synteny figures for scaffold_2 and scaffold_8.

Input files (read automatically):
  - gene_summary.tsv         → tier classifications, contig assignments
  - placed_genes.tsv          → Z-side gene coordinates on scaffolds
  - unplaced_genes.tsv        → W-side gene coordinates on contigs
  - blastn/gametolog_blastn_results.tsv → pairwise identity and coverage

Output:
  - figures/tier1a_synteny_scaffold2.{pdf,png,svg,tiff}
  - figures/tier1a_synteny_scaffold8.{pdf,png,svg,tiff}

Usage:
  python3 plot_tier1a_gametolog_synteny.py [BASE_DIR]

  BASE_DIR defaults to:
    /data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/gametolog_discovery_hanno7

Author: K. Kopp, P. euphronides genome project
Date: 2026-03-03
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import csv
from collections import defaultdict

# ============================================================
# PATHS
# ============================================================

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/gametolog_discovery_hanno7"

GENE_SUMMARY   = os.path.join(BASE_DIR, "gene_summary.tsv")
PLACED_GENES   = os.path.join(BASE_DIR, "placed_genes.tsv")
UNPLACED_GENES = os.path.join(BASE_DIR, "unplaced_genes.tsv")
BLASTN_RESULTS = os.path.join(BASE_DIR, "blastn", "gametolog_blastn_results.tsv")
OUT_DIR        = os.path.join(BASE_DIR, "figures")

os.makedirs(OUT_DIR, exist_ok=True)

for f in [GENE_SUMMARY, PLACED_GENES, UNPLACED_GENES, BLASTN_RESULTS]:
    if not os.path.isfile(f):
        sys.exit(f"ERROR: required file not found: {f}")


# ============================================================
# 1. READ GENE SUMMARY → tier1a gene names + contig pairs
# ============================================================
# Columns (tab-separated, no header):
#   gene_name  [7 count cols: indices 1-7]  scaffold[8]  z_contig_info[9]  w_contig_info[10]  tier[11]
# z_contig_info: contig_N:scaffold_X:start-end:flag
# w_contig_info: contig_N:coverage

tier1a_genes = {}

print("Reading gene_summary.tsv ...")
with open(GENE_SUMMARY) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        tier = cols[-1]
        if not tier.startswith("tier1a"):
            continue

        gene_name = cols[0]
        scaffold  = cols[8]
        z_info    = cols[9]   # contig_N:scaffold_X:start-end:flag
        w_info    = cols[10]  # contig_N:coverage

        z_parts    = z_info.split(":")
        z_contig   = z_parts[0]
        z_coords_s = z_parts[2].split("-")
        z_scf_start = int(z_coords_s[0])
        z_scf_end   = int(z_coords_s[1])

        w_parts  = w_info.split(":")
        w_contig = w_parts[0]
        w_cov    = float(w_parts[1])

        tier1a_genes[gene_name] = {
            "scaffold":    scaffold,
            "z_contig":    z_contig,
            "z_scf_start": z_scf_start,
            "z_scf_end":   z_scf_end,
            "w_contig":    w_contig,
            "w_cov":       w_cov,
        }

print(f"  Found {len(tier1a_genes)} tier1a genes: {sorted(tier1a_genes.keys())}")


# ============================================================
# 2. READ PLACED_GENES → Z-side gene coordinates
# ============================================================
# Header: contig  gene_start  gene_end  gene_id  strand  gene_name  source
#         orflen  mRNAlen  scaffold  scf_start  scf_end  mean_cov  cov_class  z_flag

z_gene_coords = {}

print("Reading placed_genes.tsv ...")
with open(PLACED_GENES) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gname = row["gene_name"]
        if gname in tier1a_genes:
            z_gene_coords[gname] = {
                "contig":     row["contig"],
                "gene_start": int(row["gene_start"]),
                "gene_end":   int(row["gene_end"]),
                "strand":     row["strand"],
            }

print(f"  Found Z coordinates for {len(z_gene_coords)} genes")


# ============================================================
# 3. READ UNPLACED_GENES → W-side gene coordinates
# ============================================================
# Header: contig  gene_start  gene_end  gene_id  strand  gene_name  source
#         orflen  mRNAlen  mean_cov  cov_class

w_gene_coords = {}

print("Reading unplaced_genes.tsv ...")
with open(UNPLACED_GENES) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gname = row["gene_name"]
        if gname in tier1a_genes:
            w_gene_coords[gname] = {
                "contig":     row["contig"],
                "gene_start": int(row["gene_start"]),
                "gene_end":   int(row["gene_end"]),
                "strand":     row["strand"],
            }

print(f"  Found W coordinates for {len(w_gene_coords)} genes")


# ============================================================
# 4. READ BLASTN RESULTS → identity and query coverage
# ============================================================
# Header: gene_name  tier  Z_contig ... pct_identity ... query_coverage

blastn = {}

print("Reading gametolog_blastn_results.tsv ...")
with open(BLASTN_RESULTS) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if not row["tier"].startswith("tier1a"):
            continue
        blastn[row["gene_name"]] = {
            "pct_identity":   float(row["pct_identity"]),
            "query_coverage": float(row["query_coverage"]),
        }

print(f"  Found blastn results for {len(blastn)} genes")


# ============================================================
# 5. ASSEMBLE CONTIG PAIRS
# ============================================================

def build_pairs(scaffold_filter):
    """Group tier1a genes into contig pairs for a given scaffold."""
    pair_key_to_genes = defaultdict(list)

    for gname, info in tier1a_genes.items():
        if info["scaffold"] != scaffold_filter:
            continue
        if gname not in z_gene_coords or gname not in w_gene_coords or gname not in blastn:
            print(f"  WARNING: skipping {gname} — missing data in one or more files")
            continue

        pk = (info["z_contig"], info["w_contig"])
        z = z_gene_coords[gname]
        b = blastn[gname]
        w = w_gene_coords[gname]

        pair_key_to_genes[pk].append({
            "name":  gname,
            "z_s":   z["gene_start"],
            "z_e":   z["gene_end"],
            "zstr":  z["strand"],
            "wstr":  w["strand"],
            "ident": b["pct_identity"],
            "qcov":  b["query_coverage"],
        })

    pairs = []
    for (zc, wc), genes in pair_key_to_genes.items():
        info = tier1a_genes[genes[0]["name"]]
        genes.sort(key=lambda g: g["z_s"])
        pairs.append({
            "z_contig":    zc,
            "w_contig":    wc,
            "scaffold":    scaffold_filter,
            "z_scf_start": info["z_scf_start"],
            "z_scf_end":   info["z_scf_end"],
            "w_cov":       info["w_cov"],
            "genes":       genes,
        })

    pairs.sort(key=lambda p: p["z_scf_start"])
    return pairs


scf2_pairs = build_pairs("scaffold_2")
scf8_pairs = build_pairs("scaffold_8")

print(f"\nscaffold_2: {len(scf2_pairs)} contig pairs, "
      f"{sum(len(p['genes']) for p in scf2_pairs)} genes")
for p in scf2_pairs:
    gnames = [g["name"] for g in p["genes"]]
    print(f"  {p['z_contig']} <-> {p['w_contig']} "
          f"({p['z_scf_start']/1e6:.2f}-{p['z_scf_end']/1e6:.2f} Mb): {gnames}")

print(f"scaffold_8: {len(scf8_pairs)} contig pairs, "
      f"{sum(len(p['genes']) for p in scf8_pairs)} genes")
for p in scf8_pairs:
    gnames = [g["name"] for g in p["genes"]]
    print(f"  {p['z_contig']} <-> {p['w_contig']} "
          f"({p['z_scf_start']/1e6:.2f}-{p['z_scf_end']/1e6:.2f} Mb): {gnames}")


# ============================================================
# 6. PLOTTING FUNCTIONS
# ============================================================

def draw_vertical_gene(ax, x, y_top, y_bot, w, strand, color, min_h=None):
    """Vertical arrow: + = tip down (increasing Mb), - = tip up."""
    if min_h and abs(y_bot - y_top) < min_h:
        mid = (y_top + y_bot) / 2
        y_top, y_bot = mid - min_h / 2, mid + min_h / 2
    hw = w / 2
    tip_h = abs(y_bot - y_top) * 0.18
    if strand == "+":
        xs = [x - hw, x - hw, x - hw * 1.4, x, x + hw * 1.4, x + hw, x + hw]
        ys = [y_top, y_bot - tip_h, y_bot - tip_h, y_bot,
              y_bot - tip_h, y_bot - tip_h, y_top]
    else:
        xs = [x - hw, x - hw, x - hw * 1.4, x, x + hw * 1.4, x + hw, x + hw]
        ys = [y_bot, y_top + tip_h, y_top + tip_h, y_top,
              y_top + tip_h, y_top + tip_h, y_bot]
    ax.fill(xs, ys, color=color, edgecolor="black", linewidth=0.7,
            alpha=0.85, zorder=3)


def populate_panel(ax, pairs, y_lo, y_hi, min_gene_h_mb=None,
                   show_zw_headers=True):
    """Draw genes for all contig pairs within [y_lo, y_hi] Mb."""
    ax.set_ylim(y_hi, y_lo)
    z_x, w_x, gw = 2.0, 8.0, 0.30
    if min_gene_h_mb is None:
        min_gene_h_mb = (y_hi - y_lo) * 0.02

    palette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
               "#a65628", "#f781bf", "#66c2a5", "#fc8d62", "#8da0cb",
               "#e78ac3", "#a6d854", "#ffd92f"]

    ax.axvline(x=z_x, color="#2166ac", linewidth=1.5, alpha=0.2, zorder=0)
    ax.axvline(x=w_x, color="#b2182b", linewidth=1.5, alpha=0.2, zorder=0)

    if show_zw_headers:
        ax.text(z_x, y_lo - (y_hi - y_lo) * 0.03, "Z",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color="#2166ac")
        ax.text(w_x, y_lo - (y_hi - y_lo) * 0.03, "W",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color="#b2182b")

    # Global gene index for consistent colours
    all_genes_sorted = []
    for p in pairs:
        for g in p["genes"]:
            scf_mid = (p["z_scf_start"] + g["z_s"] + p["z_scf_start"] + g["z_e"]) / 2e6
            all_genes_sorted.append((scf_mid, g["name"]))
    all_genes_sorted.sort()
    gidx = {name: i for i, (_, name) in enumerate(all_genes_sorted)}

    for p in pairs:
        scf_s = p["z_scf_start"]
        scf_e = p["z_scf_end"]
        ct_mb, cb_mb = scf_s / 1e6, scf_e / 1e6

        if cb_mb < y_lo or ct_mb > y_hi:
            continue

        # Z contig bracket
        ct_c, cb_c = max(ct_mb, y_lo), min(cb_mb, y_hi)
        ax.add_patch(plt.Rectangle(
            (z_x - 0.45, ct_c), 0.9, cb_c - ct_c,
            facecolor="#2166ac", edgecolor="none", alpha=0.06, zorder=1))
        ax.text(z_x - 0.6, (ct_c + cb_c) / 2, p["z_contig"],
                ha="right", va="center", fontsize=6.5, color="#2166ac",
                fontstyle="italic")

        w_ys = []

        for gene in p["genes"]:
            zt_mb = (scf_s + gene["z_s"]) / 1e6
            zb_mb = (scf_s + gene["z_e"]) / 1e6
            zm = (zt_mb + zb_mb) / 2

            if zm < y_lo or zm > y_hi:
                continue

            gi = gidx.get(gene["name"], 0)
            gc = palette[gi % len(palette)]

            draw_vertical_gene(ax, z_x, zt_mb, zb_mb, gw, gene["zstr"], gc,
                              min_h=min_gene_h_mb)
            ax.text(z_x + gw / 2 + 0.1, zm, gene["name"],
                    ha="left", va="center", fontsize=7.5, fontweight="bold",
                    fontstyle="italic")

            draw_vertical_gene(ax, w_x, zt_mb, zb_mb, gw, gene["wstr"], gc,
                              min_h=min_gene_h_mb)
            ax.text(w_x - gw / 2 - 0.1, zm, gene["name"],
                    ha="right", va="center", fontsize=7.5, fontweight="bold",
                    fontstyle="italic")

            ax.plot([z_x + gw / 2 + 0.02, w_x - gw / 2 - 0.02],
                    [zm, zm], color="grey", linewidth=0.5, alpha=0.4, zorder=2)

            mx = (z_x + w_x) / 2
            ax.text(mx, zm,
                    f"{gene['ident']:.1f}% id\n{gene['qcov']:.0f}% qcov",
                    ha="center", va="center", fontsize=6.5,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="grey", alpha=0.95, linewidth=0.5),
                    zorder=5)

            w_ys.append((zt_mb, zb_mb))

        if w_ys:
            wt = min(y[0] for y in w_ys)
            wb = max(y[1] for y in w_ys)
            m = min_gene_h_mb * 0.5
            ax.add_patch(plt.Rectangle(
                (w_x - 0.45, wt - m), 0.9, wb - wt + 2 * m,
                facecolor="#b2182b", edgecolor="none", alpha=0.06, zorder=1))
            ax.text(w_x + 0.6, (wt + wb) / 2, p["w_contig"],
                    ha="left", va="center", fontsize=6.5, color="#b2182b",
                    fontstyle="italic")

    ax.set_xlim(-0.5, 10.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False, labelsize=9)


def save_figure(fig, basename):
    for fmt in ["pdf", "png", "svg", "tiff"]:
        kwargs = {"dpi": 300, "bbox_inches": "tight"}
        if fmt == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        outpath = os.path.join(OUT_DIR, f"{basename}.{fmt}")
        fig.savefig(outpath, **kwargs)
        print(f"Saved: {outpath}")


def compute_windows(pairs, padding_frac=0.3):
    """Compute zoom windows per contig pair, padded around gene extent."""
    windows = []
    for p in pairs:
        gene_starts = [(p["z_scf_start"] + g["z_s"]) / 1e6 for g in p["genes"]]
        gene_ends   = [(p["z_scf_start"] + g["z_e"]) / 1e6 for g in p["genes"]]
        lo, hi = min(gene_starts), max(gene_ends)
        span = hi - lo
        pad = max(span * padding_frac, 0.003)
        windows.append((lo - pad, hi + pad))
    return windows


# ============================================================
# 7. FIGURE A: scaffold_2
# ============================================================

if scf2_pairs:
    print("\nGenerating scaffold_2 figure ...")
    windows = compute_windows(scf2_pairs)
    ranges_mb = [hi - lo for lo, hi in windows]
    gaps_mb = [windows[i+1][0] - windows[i][1] for i in range(len(windows) - 1)]

    n_panels = len(windows)
    n_gaps = n_panels - 1
    n_rows = n_panels + n_gaps

    gap_h = 0.004
    h_ratios = []
    for i, r in enumerate(ranges_mb):
        h_ratios.append(r)
        if i < n_gaps:
            h_ratios.append(gap_h)

    fig_a = plt.figure(figsize=(10, 22))
    gs_a = gridspec.GridSpec(n_rows, 1, height_ratios=h_ratios, hspace=0.05)

    min_h = 0.002

    for i, (y_lo, y_hi) in enumerate(windows):
        row_idx = i * 2

        ax = fig_a.add_subplot(gs_a[row_idx])
        populate_panel(ax, scf2_pairs, y_lo, y_hi, min_gene_h_mb=min_h,
                       show_zw_headers=(i == 0))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.3f}"))

        if i == 0:
            ax.text(-0.02, 1.05, "Mb", ha="center", va="bottom", fontsize=9,
                    fontweight="bold", transform=ax.transAxes)

        if i < n_gaps:
            ax_g = fig_a.add_subplot(gs_a[row_idx + 1])
            ax_g.set_xlim(0, 10)
            ax_g.set_ylim(0, 1)
            ax_g.plot([0, 10], [0.5, 0.5], ':', color="grey", linewidth=1.0)
            gv = gaps_mb[i]
            gl = f"~{gv:.1f} Mb" if gv >= 1 else f"~{gv:.2f} Mb"
            ax_g.text(5.0, 0.55, gl, ha="center", va="bottom",
                      fontsize=7.5, color="grey", fontstyle="italic")
            ax_g.axis("off")

    fig_a.text(0.02, 0.5, "scaffold_2", ha="center", va="center",
               fontsize=10, rotation=90)

    save_figure(fig_a, "tier1a_synteny_scaffold2")
    plt.close(fig_a)


# ============================================================
# 8. FIGURE B: scaffold_8
# ============================================================

if scf8_pairs:
    print("\nGenerating scaffold_8 figure ...")
    windows_8 = compute_windows(scf8_pairs)
    y_lo, y_hi = windows_8[0]

    fig_b = plt.figure(figsize=(10, 8))
    ax_b = fig_b.add_subplot(111)
    populate_panel(ax_b, scf8_pairs, y_lo, y_hi, min_gene_h_mb=0.003)
    ax_b.set_ylabel("scaffold_8", fontsize=10)
    ax_b.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax_b.text(-0.02, 1.05, "Mb", ha="center", va="bottom", fontsize=9,
              fontweight="bold", transform=ax_b.transAxes)

    save_figure(fig_b, "tier1a_synteny_scaffold8")
    plt.close(fig_b)

print("\nDone.")
