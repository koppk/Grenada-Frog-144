#!/usr/bin/env python3
"""
plot_tier1a_gametolog_banding.py
================================
Banding-profile figures for tier 1a gametolog contig pairs, matching
the layout of plot_tier1a_gametolog_synteny.py: vertical contig strips,
Z left, W right, y-axis = scaffold Mb position.

Contig rectangles are painted with compositional banding (GC content
+ WindowMasker masking density). Gene positions are marked as black
ticks at the same y-position on Z and W (scaffold coordinates).

BMC Genomics: full page width 170 mm, max height 225 mm (incl. legend),
  TIFF LZW, 300 DPI, ≤10 MB, Arial, closely cropped.

Usage:
  python3 plot_tier1a_gametolog_banding.py [BASE_DIR] [ASSEMBLY_FASTA]

Author: Kopp K. Pristimantis euphronides genome project.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.gridspec as gridspec
import numpy as np
import os
import sys
import csv
from collections import defaultdict

# ============================================================
# PATHS
# ============================================================

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/gametolog_discovery_hanno7"

ASM_FASTA = sys.argv[2] if len(sys.argv) > 2 else \
    "/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"

WM_PLACED = os.path.join(os.path.dirname(BASE_DIR),
                         "windowmasker", "wm_placed_intervals.txt")
WM_UNPLACED = os.path.join(os.path.dirname(BASE_DIR),
                           "windowmasker", "wm_unplaced_intervals.txt")

GENE_SUMMARY   = os.path.join(BASE_DIR, "gene_summary.tsv")
PLACED_GENES   = os.path.join(BASE_DIR, "placed_genes.tsv")
UNPLACED_GENES = os.path.join(BASE_DIR, "unplaced_genes.tsv")
AGP_FILE       = "/data/GrenadaFrog144/ragtag.scaffold.renamed.agp"
OUT_DIR        = os.path.join(BASE_DIR, "figures")

os.makedirs(OUT_DIR, exist_ok=True)

for f in [GENE_SUMMARY, PLACED_GENES, UNPLACED_GENES, ASM_FASTA]:
    if not os.path.isfile(f):
        sys.exit(f"ERROR: required file not found: {f}")

# ============================================================
# BMC Genomics figure dimensions
# ============================================================
# Full page width: 170 mm = 6.693 in
# Max height (figure + legend): 225 mm = 8.858 in

BMC_WIDTH_MM  = 170
BMC_HEIGHT_MM = 225
BMC_WIDTH_IN  = BMC_WIDTH_MM / 25.4
BMC_HEIGHT_IN = BMC_HEIGHT_MM / 25.4

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 8,
    'axes.titlesize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'lines.linewidth': 0.8,
    'axes.linewidth': 0.5,
})

Z_COL = "#2166ac"
W_COL = "#b2182b"


# ============================================================
# SEQUENCE ANALYSIS
# ============================================================

def sliding_window_gc(seq, window_size, step_size):
    n = len(seq)
    gc_arr = np.array([1.0 if b in 'GC' else 0.0 for b in seq], dtype=np.float32)
    valid_arr = np.array([0.0 if b == 'N' else 1.0 for b in seq], dtype=np.float32)
    gc_cum = np.cumsum(gc_arr)
    valid_cum = np.cumsum(valid_arr)
    gc_values, positions = [], []
    for i in range(0, n - window_size + 1, step_size):
        end = i + window_size
        gc_sum = gc_cum[end-1] - (gc_cum[i-1] if i > 0 else 0)
        valid_sum = valid_cum[end-1] - (valid_cum[i-1] if i > 0 else 0)
        gc_values.append(gc_sum / valid_sum if valid_sum > 0 else 0.5)
        positions.append(i + window_size // 2)
    return np.array(positions), np.array(gc_values)


def composite_ideogram_color(gc, wm_density):
    """Three-colour composite banding using WindowMasker density as
    repeat indicator and GC for AT/GC classification.
    Consistent with master table classification method."""
    if wm_density is None or wm_density <= 0.4:
        return (0.93, 0.93, 0.93)       # Non-repetitive
    if gc < 0.50:
        return (0.75, 0.92, 0.15)       # AT-rich repetitive
    else:
        return (0.85, 0.65, 0.1)        # GC-rich repetitive


def parse_windowmasker_intervals(filepath):
    if not os.path.isfile(filepath):
        return {}
    masked = defaultdict(list)
    cur = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                cur = line[1:].split()[0]
            elif ' - ' in line and cur:
                parts = line.split(' - ')
                try:
                    masked[cur].append((int(parts[0].strip()), int(parts[1].strip())))
                except ValueError:
                    continue
    return dict(masked)


def windowmasker_to_density(intervals, seq_len, positions, window_size):
    mask = np.zeros(seq_len, dtype=np.float32)
    for s, e in intervals:
        mask[max(0, s):min(e+1, seq_len)] = 1.0
    hw = window_size // 2
    return np.array([np.mean(mask[max(0, p-hw):min(seq_len, p+hw)])
                     for p in positions])


# ============================================================
# FASTA + AGP
# ============================================================

def extract_contigs(fasta_path, names):
    seqs = {}
    target = set(names)
    cur, buf = None, []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if cur in target:
                    seqs[cur] = ''.join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
    if cur in target:
        seqs[cur] = ''.join(buf)
    return seqs


def read_agp_orientations(agp_path):
    orient = {}
    if not os.path.isfile(agp_path):
        return orient
    with open(agp_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.strip().split('\t')
            if len(c) >= 9 and c[4] == 'W':
                orient[c[5]] = c[8]
    return orient


# ============================================================
# READ GAMETOLOG DATA
# ============================================================

tier1a_genes = {}
print("Reading gene_summary.tsv ...")
with open(GENE_SUMMARY) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if not cols[-1].startswith("tier1a"):
            continue
        gn = cols[0]
        z_parts = cols[9].split(":")
        w_parts = cols[10].split(":")
        tier1a_genes[gn] = {
            "scaffold": cols[8],
            "z_contig": z_parts[0],
            "z_scf_start": int(z_parts[2].split("-")[0]),
            "z_scf_end":   int(z_parts[2].split("-")[1]),
            "w_contig": w_parts[0],
            "w_cov": float(w_parts[1]),
        }
print(f"  {len(tier1a_genes)} tier1a genes")

z_gene_coords, w_gene_coords = {}, {}
print("Reading placed_genes.tsv ...")
with open(PLACED_GENES) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["gene_name"] in tier1a_genes:
            z_gene_coords[row["gene_name"]] = {
                "contig": row["contig"],
                "start": int(row["gene_start"]),
                "end":   int(row["gene_end"]),
                "strand": row["strand"],
            }

print("Reading unplaced_genes.tsv ...")
with open(UNPLACED_GENES) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["gene_name"] in tier1a_genes:
            w_gene_coords[row["gene_name"]] = {
                "contig": row["contig"],
                "start": int(row["gene_start"]),
                "end":   int(row["gene_end"]),
                "strand": row["strand"],
            }

print(f"  Z coords: {len(z_gene_coords)}, W coords: {len(w_gene_coords)}")


# ============================================================
# BUILD PAIRS
# ============================================================

def build_pairs(scaffold_filter):
    pk2g = defaultdict(list)
    for gn, info in tier1a_genes.items():
        if info["scaffold"] != scaffold_filter:
            continue
        if gn not in z_gene_coords or gn not in w_gene_coords:
            continue
        z, w = z_gene_coords[gn], w_gene_coords[gn]
        pk2g[(info["z_contig"], info["w_contig"])].append({
            "name": gn,
            "z_s": z["start"], "z_e": z["end"], "zstr": z["strand"],
            "w_s": w["start"], "w_e": w["end"], "wstr": w["strand"],
        })
    pairs = []
    for (zc, wc), genes in pk2g.items():
        info = tier1a_genes[genes[0]["name"]]
        genes.sort(key=lambda g: g["z_s"])
        pairs.append({
            "z_contig": zc, "w_contig": wc,
            "scaffold": scaffold_filter,
            "z_scf_start": info["z_scf_start"],
            "z_scf_end": info["z_scf_end"],
            "w_cov": info["w_cov"],
            "genes": genes,
        })
    pairs.sort(key=lambda p: p["z_scf_start"])
    return pairs


scf2_pairs = build_pairs("scaffold_2")
scf8_pairs = build_pairs("scaffold_8")
print(f"Pairs: scaffold_2={len(scf2_pairs)}, scaffold_8={len(scf8_pairs)}")


# ============================================================
# EXTRACT SEQUENCES + WINDOWMASKER
# ============================================================

needed = set()
for p in scf2_pairs + scf8_pairs:
    needed.add(p["z_contig"])
    needed.add(p["w_contig"])

print(f"\nExtracting {len(needed)} contigs ...")
sequences = extract_contigs(ASM_FASTA, needed)
agp_orient = read_agp_orientations(AGP_FILE)

wm_intervals = {}
for wf in [WM_PLACED, WM_UNPLACED]:
    wm_intervals.update(parse_windowmasker_intervals(wf))
print(f"  WindowMasker: {len(wm_intervals)} contigs")


# ============================================================
# BANDING COMPUTATION
# ============================================================

WIN, STP = 2000, 400


def compute_banding(contig_name):
    seq = sequences.get(contig_name)
    if seq is None:
        return None
    clen = len(seq)
    ew = max(200, min(WIN, clen // 5))
    es = max(100, min(STP, ew // 2))
    pos, gc = sliding_window_gc(seq, ew, es)
    wm = None
    if contig_name in wm_intervals:
        wm = windowmasker_to_density(wm_intervals[contig_name], clen, pos, ew)
    return {"length": clen, "pos": pos, "gc": gc, "wm": wm}


# ============================================================
# DRAWING
# ============================================================

def paint_banding_rect(ax, banding, x_center, rect_width,
                       y_top_mb, y_bot_mb, side_color, is_reverse=False):
    clen = banding["length"]
    pos = banding["pos"]
    gc = banding["gc"]
    wm = banding["wm"]
    n = len(pos)
    if n < 2:
        return

    hw = rect_width / 2
    patches, colors = [], []
    for i in range(n - 1):
        frac = pos[i] / clen
        frac_next = pos[i+1] / clen
        if is_reverse:
            frac = 1.0 - frac
            frac_next = 1.0 - frac_next
            if frac > frac_next:
                frac, frac_next = frac_next, frac
        y_t = y_top_mb + frac * (y_bot_mb - y_top_mb)
        y_b = y_top_mb + frac_next * (y_bot_mb - y_top_mb)
        rect = plt.Rectangle((x_center - hw, y_t), rect_width, y_b - y_t)
        patches.append(rect)
        wm_d = wm[i] if wm is not None else 0.0
        colors.append(composite_ideogram_color(gc[i], wm_d))

    ax.add_collection(PatchCollection(patches, facecolors=colors,
                                      edgecolors='none', linewidths=0))

    ax.add_patch(mpatches.FancyBboxPatch(
        (x_center - hw, y_top_mb), rect_width, y_bot_mb - y_top_mb,
        boxstyle="round,pad=0.001",
        linewidth=1.0, edgecolor=side_color, facecolor='none', zorder=10))


# ============================================================
# MAIN PANEL
# ============================================================

def populate_banding_panel(ax, pairs, y_lo, y_hi, show_zw_headers=True):
    ax.set_ylim(y_hi, y_lo)

    z_x, w_x = 2.0, 8.0
    strip_w = 0.9

    ax.axvline(x=z_x, color=Z_COL, linewidth=1.5, alpha=0.15, zorder=0)
    ax.axvline(x=w_x, color=W_COL, linewidth=1.5, alpha=0.15, zorder=0)

    if show_zw_headers:
        ax.text(z_x, y_lo - (y_hi - y_lo) * 0.03, "Z",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color=Z_COL)
        ax.text(w_x, y_lo - (y_hi - y_lo) * 0.03, "W",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color=W_COL)

    for p in pairs:
        scf_s = p["z_scf_start"]
        scf_e = p["z_scf_end"]
        ct_mb, cb_mb = scf_s / 1e6, scf_e / 1e6

        if cb_mb < y_lo or ct_mb > y_hi:
            continue

        ct_c = max(ct_mb, y_lo)
        cb_c = min(cb_mb, y_hi)

        z_name = p["z_contig"]
        w_name = p["w_contig"]
        z_reverse = (agp_orient.get(z_name, '+') == '-')

        z_band = compute_banding(z_name)
        w_band = compute_banding(w_name)

        if z_band is None or w_band is None:
            print(f"  Skipping {z_name}/{w_name} — not in FASTA")
            continue

        # Z banding
        paint_banding_rect(ax, z_band, z_x, strip_w,
                           ct_c, cb_c, Z_COL, is_reverse=z_reverse)
        ax.text(z_x - 0.6, (ct_c + cb_c) / 2, z_name,
                ha="right", va="center", fontsize=6.5,
                color=Z_COL, fontstyle="italic")

        # W banding
        paint_banding_rect(ax, w_band, w_x, strip_w,
                           ct_c, cb_c, W_COL, is_reverse=False)
        ax.text(w_x + 0.6, (ct_c + cb_c) / 2, w_name,
                ha="left", va="center", fontsize=6.5,
                color=W_COL, fontstyle="italic")

        # Gene positions — unfilled black rectangles spanning gene length
        min_gene_h_mb = (cb_c - ct_c) * 0.02  # minimum visible height
        for g in p["genes"]:
            zt_mb = (scf_s + g["z_s"]) / 1e6
            zb_mb = (scf_s + g["z_e"]) / 1e6
            zm = (zt_mb + zb_mb) / 2

            if zm < y_lo or zm > y_hi:
                continue

            # Ensure minimum visible height
            if abs(zb_mb - zt_mb) < min_gene_h_mb:
                mid = (zt_mb + zb_mb) / 2
                zt_mb = mid - min_gene_h_mb / 2
                zb_mb = mid + min_gene_h_mb / 2

            gene_h = zb_mb - zt_mb

            # Z gene — unfilled rectangle on the Z strip
            ax.add_patch(plt.Rectangle(
                (z_x - strip_w/2, zt_mb), strip_w, gene_h,
                facecolor='none', edgecolor='black',
                linewidth=1.5, linestyle='--', zorder=15))
            ax.text(z_x + strip_w/2 + 0.1, zm, g["name"],
                    ha="left", va="center", fontsize=7,
                    fontstyle="italic", color="black")

            # W gene — unfilled rectangle on the W strip (same y-position)
            ax.add_patch(plt.Rectangle(
                (w_x - strip_w/2, zt_mb), strip_w, gene_h,
                facecolor='none', edgecolor='black',
                linewidth=1.5, linestyle='--', zorder=15))
            ax.text(w_x - strip_w/2 - 0.1, zm, g["name"],
                    ha="right", va="center", fontsize=7,
                    fontstyle="italic", color="black")

    ax.set_xlim(-0.5, 10.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False, labelsize=9)


def compute_windows(pairs, padding_frac=0.3):
    windows = []
    for p in pairs:
        gene_starts = [(p["z_scf_start"] + g["z_s"]) / 1e6 for g in p["genes"]]
        gene_ends   = [(p["z_scf_start"] + g["z_e"]) / 1e6 for g in p["genes"]]
        lo, hi = min(gene_starts), max(gene_ends)
        span = hi - lo
        pad = max(span * padding_frac, 0.003)
        windows.append((lo - pad, hi + pad))
    return windows


def save_figure(fig, basename):
    for fmt in ["pdf", "png", "tiff"]:
        kwargs = {"dpi": 300, "bbox_inches": "tight"}
        if fmt == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        outpath = os.path.join(OUT_DIR, f"{basename}.{fmt}")
        fig.savefig(outpath, **kwargs)
        print(f"  Saved: {outpath}")


# ============================================================
# LEGEND (reusable)
# ============================================================

LEGEND_PATCHES = [
    mpatches.Patch(color=(0.75, 0.92, 0.15), label='AT-rich repetitive'),
    mpatches.Patch(color=(0.85, 0.65, 0.1), label='GC-rich repetitive'),
    mpatches.Patch(facecolor=(0.93, 0.93, 0.93), edgecolor='#ccc',
                   label='Non-repetitive'),
]


# ============================================================
# SCAFFOLD 2
# ============================================================

if scf2_pairs:
    print("\nGenerating scaffold_2 banding figure ...")
    windows = compute_windows(scf2_pairs)
    ranges_mb = [hi - lo for lo, hi in windows]
    gaps_mb = [windows[i+1][0] - windows[i][1]
               for i in range(len(windows) - 1)]

    n_panels = len(windows)
    n_gaps = n_panels - 1
    n_rows = n_panels + n_gaps
    gap_h = 0.004

    h_ratios = []
    for i, r in enumerate(ranges_mb):
        h_ratios.append(r)
        if i < n_gaps:
            h_ratios.append(gap_h)

    fig_a = plt.figure(figsize=(BMC_WIDTH_IN, BMC_HEIGHT_IN))
    gs_a = gridspec.GridSpec(n_rows, 1, height_ratios=h_ratios, hspace=0.05)

    for i, (y_lo, y_hi) in enumerate(windows):
        row_idx = i * 2
        ax = fig_a.add_subplot(gs_a[row_idx])
        populate_banding_panel(ax, scf2_pairs, y_lo, y_hi,
                               show_zw_headers=(i == 0))
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{x:.3f}"))
        if i == 0:
            ax.text(-0.02, 1.05, "Mb", ha="center", va="bottom",
                    fontsize=9, fontweight="bold", transform=ax.transAxes)

        if i < n_gaps:
            ax_g = fig_a.add_subplot(gs_a[row_idx + 1])
            ax_g.set_xlim(0, 10); ax_g.set_ylim(0, 1)
            ax_g.plot([0, 10], [0.5, 0.5], ':', color="grey", lw=1.0)
            gv = gaps_mb[i]
            gl = f"~{gv:.1f} Mb" if gv >= 1 else f"~{gv:.2f} Mb"
            ax_g.text(5.0, 0.55, gl, ha="center", va="bottom",
                      fontsize=7.5, color="grey", fontstyle="italic")
            ax_g.axis("off")

    fig_a.text(0.01, 0.5, "scaffold_2", ha="center", va="center",
               fontsize=10, rotation=90)

    fig_a.legend(handles=LEGEND_PATCHES, loc="lower center", ncol=3,
                 fontsize=6, frameon=True, borderpad=0.3,
                 handletextpad=0.4, columnspacing=0.8)

    plt.subplots_adjust(bottom=0.045, top=0.97, left=0.12, right=0.95)

    save_figure(fig_a, "tier1a_banding_scaffold2")
    plt.close(fig_a)


# ============================================================
# SCAFFOLD 8
# ============================================================

if scf8_pairs:
    print("\nGenerating scaffold_8 banding figure ...")
    windows_8 = compute_windows(scf8_pairs)
    y_lo, y_hi = windows_8[0]

    # Single pair — half-page height is enough
    fig_h = min(BMC_HEIGHT_IN, BMC_WIDTH_IN * 0.6)
    fig_b = plt.figure(figsize=(BMC_WIDTH_IN, fig_h))
    ax_b = fig_b.add_subplot(111)
    populate_banding_panel(ax_b, scf8_pairs, y_lo, y_hi)
    ax_b.set_ylabel("scaffold_8", fontsize=10)
    ax_b.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax_b.text(-0.02, 1.05, "Mb", ha="center", va="bottom",
              fontsize=9, fontweight="bold", transform=ax_b.transAxes)

    fig_b.legend(handles=LEGEND_PATCHES, loc="lower center", ncol=3,
                 fontsize=6, frameon=True, borderpad=0.3,
                 handletextpad=0.4, columnspacing=0.8)

    plt.subplots_adjust(bottom=0.10, top=0.95, left=0.12, right=0.95)

    save_figure(fig_b, "tier1a_banding_scaffold8")
    plt.close(fig_b)


print("\nDone.")
