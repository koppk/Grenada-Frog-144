#!/usr/bin/env python3
"""
plot_w_genic_banding.py
========================
Composite banding of W-genic contigs (tier 1a–c).

Same three-colour scheme, WindowMasker intervals, 2 kb / 400 bp
windows, and WM_THRESHOLD=0.4 as all other banding figures.

Panels A (tier1a, 5), B (tier1b, 6 longest of 29), C (tier1c, 6).

Author: Kopp K. Pristimantis euphronides genome project.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import sys
import csv
from collections import defaultdict

# ============================================================
# PATHS
# ============================================================

ASM_FASTA = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"

MASTER_TABLE = sys.argv[2] if len(sys.argv) > 2 else \
    "/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/master_table/contig_master_table.tsv"

BASE = "/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"
WM_PLACED = os.path.join(BASE, "windowmasker", "wm_placed_intervals.txt")
WM_UNPLACED = os.path.join(BASE, "windowmasker", "wm_unplaced_intervals.txt")

OUT_DIR = os.path.join(BASE, "gametolog_discovery_hanno7", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

WINDOW = 2000
STEP   = 400
WM_THRESHOLD = 0.4

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
})

COL_AT_REP = (0.75, 0.92, 0.15)
COL_GC_REP = (0.85, 0.65, 0.1)
COL_REP_POOR = (0.93, 0.93, 0.93)


def composite_colour(gc, wm_density):
    if wm_density is None or wm_density <= WM_THRESHOLD:
        return COL_REP_POOR
    return COL_AT_REP if gc < 0.50 else COL_GC_REP


def classify_windows(gc_arr, wm_arr):
    """Classify each window into three banding categories.
    Returns dict with counts, fractions, and means."""
    n = len(gc_arr)
    n_at = n_gc = n_nr = 0
    for i in range(n):
        wm_d = wm_arr[i] if wm_arr is not None else 0.0
        if wm_d <= WM_THRESHOLD:
            n_nr += 1
        elif gc_arr[i] < 0.50:
            n_at += 1
        else:
            n_gc += 1
    return {
        "n_windows": n,
        "gc_mean": float(np.mean(gc_arr)) if n > 0 else 0.0,
        "wm_mean": float(np.mean(wm_arr)) if (wm_arr is not None and n > 0) else 0.0,
        "n_at_rep": n_at, "n_gc_rep": n_gc, "n_nonrep": n_nr,
        "frac_at_rep": n_at / n if n > 0 else 0.0,
        "frac_gc_rep": n_gc / n if n > 0 else 0.0,
        "frac_nonrep": n_nr / n if n > 0 else 0.0,
    }


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
                    masked[cur].append((int(parts[0].strip()),
                                        int(parts[1].strip())))
                except ValueError:
                    continue
    return dict(masked)


def windowmasker_to_density(intervals, seq_len, positions, window_size):
    mask = np.zeros(seq_len, dtype=np.float32)
    for s, e in intervals:
        mask[max(0, s):min(e + 1, seq_len)] = 1.0
    hw = window_size // 2
    return np.array([np.mean(mask[max(0, p - hw):min(seq_len, p + hw)])
                     for p in positions])


def extract_contigs(fasta_path, contig_names):
    sequences = {}
    target = set(contig_names)
    cur, buf = None, []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if cur in target:
                    sequences[cur] = ''.join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
    if cur in target:
        sequences[cur] = ''.join(buf)
    return sequences


def sliding_window_gc(seq, window_size, step_size):
    n = len(seq)
    gc_arr = np.array([1.0 if b in 'GC' else 0.0 for b in seq], dtype=np.float32)
    valid_arr = np.array([0.0 if b == 'N' else 1.0 for b in seq], dtype=np.float32)
    gc_cum = np.cumsum(gc_arr)
    valid_cum = np.cumsum(valid_arr)
    gc_values, positions = [], []
    for i in range(0, n - window_size + 1, step_size):
        end = i + window_size
        gc_sum = gc_cum[end - 1] - (gc_cum[i - 1] if i > 0 else 0)
        valid_sum = valid_cum[end - 1] - (valid_cum[i - 1] if i > 0 else 0)
        gc_values.append(gc_sum / valid_sum if valid_sum > 0 else 0.5)
        positions.append(i + window_size // 2)
    return np.array(positions), np.array(gc_values)


def render_banding_image(pos, gc, wm, seq_len, n_pixels):
    gc_sum = np.zeros(n_pixels, dtype=np.float64)
    wm_sum = np.zeros(n_pixels, dtype=np.float64)
    counts = np.zeros(n_pixels, dtype=np.float64)
    for i in range(len(pos)):
        pix = min(int(pos[i] / seq_len * n_pixels), n_pixels - 1)
        gc_sum[pix] += gc[i]
        wm_sum[pix] += (wm[i] if wm is not None else 0.0)
        counts[pix] += 1
    img = np.full((n_pixels, 3), COL_REP_POOR, dtype=np.float32)
    for p in range(n_pixels):
        if counts[p] > 0:
            img[p] = composite_colour(gc_sum[p] / counts[p],
                                       wm_sum[p] / counts[p])
    return img


# ============================================================
# SELECT CONTIGS
# ============================================================

def select_contigs():
    tier_contigs = {'tier1a': [], 'tier1b': [], 'tier1c': []}
    with open(MASTER_TABLE) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            uc = row.get("unplaced_class", "").strip()
            conf = row.get("unplaced_confidence", "").strip()
            if uc == "W-genic" and conf in tier_contigs:
                tier_contigs[conf].append({
                    "contig": row["contig"],
                    "length": int(row["length"]),
                })
    for tier in tier_contigs:
        tier_contigs[tier].sort(key=lambda x: x["length"], reverse=True)
    tier_contigs["tier1b"] = tier_contigs["tier1b"][:6]
    for tier in ["tier1a", "tier1b", "tier1c"]:
        n = len(tier_contigs[tier])
        bp = sum(c["length"] for c in tier_contigs[tier])
        print(f"  {tier}: {n} contigs, {bp / 1000:.1f} kb")
    return tier_contigs


# ============================================================
# FIGURE — manual positioning for clean layout
# ============================================================

print("Selecting W-genic contigs ...")
tier_contigs = select_contigs()

all_contig_names = []
for tier in ["tier1a", "tier1b", "tier1c"]:
    for c in tier_contigs[tier]:
        all_contig_names.append(c["contig"])

print(f"\nExtracting {len(all_contig_names)} contigs ...")
sequences = extract_contigs(ASM_FASTA, all_contig_names)
print(f"  Got {len(sequences)} contigs")

print("Loading WindowMasker intervals ...")
wm_intervals = {}
for wf in [WM_PLACED, WM_UNPLACED]:
    wm_intervals.update(parse_windowmasker_intervals(wf))
print(f"  {len(wm_intervals)} contigs with WM data")

print("Computing profiles ...")
profiles = {}
for cname in all_contig_names:
    if cname not in sequences:
        print(f"  WARNING: {cname} not in FASTA")
        continue
    seq = sequences[cname]
    clen = len(seq)
    ew = max(200, min(WINDOW, clen // 5))
    es = max(100, min(STEP, ew // 2))
    pos, gc = sliding_window_gc(seq, ew, es)
    wm = None
    if cname in wm_intervals:
        wm = windowmasker_to_density(wm_intervals[cname], clen, pos, ew)
    else:
        print(f"  WARNING: {cname} no WM data")
    profiles[cname] = {"pos": pos, "gc": gc, "wm": wm, "length": clen}

panels = []
for tier in ["tier1a", "tier1b", "tier1c"]:
    panels.append([c for c in tier_contigs[tier] if c["contig"] in profiles])

panel_labels = ['A', 'B', 'C']
total = sum(len(p) for p in panels)

# ---- TSV summary ----
tier_names = ["tier1a", "tier1b", "tier1c"]
tsv_path = os.path.join(OUT_DIR, "w_genic_banding_summary.tsv")
print(f"Writing summary: {tsv_path}")
with open(tsv_path, "w") as fout:
    hdr = ["tier", "contig", "length_bp", "n_windows", "gc_mean",
           "wm_mean", "frac_AT_rep", "frac_GC_rep", "frac_nonrep",
           "n_AT_rep", "n_GC_rep", "n_nonrep"]
    fout.write("\t".join(hdr) + "\n")
    for tier, panel_list in zip(tier_names, panels):
        tot = {"nw": 0, "bp": 0, "n_at": 0, "n_gc": 0, "n_nr": 0,
               "gc_wt": 0.0, "wm_wt": 0.0}
        for cinfo in panel_list:
            cname = cinfo["contig"]
            prof = profiles[cname]
            st = classify_windows(prof["gc"], prof["wm"])
            fout.write("\t".join([
                tier, cname, str(prof["length"]), str(st["n_windows"]),
                f"{st['gc_mean']:.4f}", f"{st['wm_mean']:.4f}",
                f"{st['frac_at_rep']:.4f}", f"{st['frac_gc_rep']:.4f}",
                f"{st['frac_nonrep']:.4f}",
                str(st["n_at_rep"]), str(st["n_gc_rep"]),
                str(st["n_nonrep"]),
            ]) + "\n")
            tot["nw"] += st["n_windows"]
            tot["bp"] += prof["length"]
            tot["n_at"] += st["n_at_rep"]
            tot["n_gc"] += st["n_gc_rep"]
            tot["n_nr"] += st["n_nonrep"]
            tot["gc_wt"] += st["gc_mean"] * st["n_windows"]
            tot["wm_wt"] += st["wm_mean"] * st["n_windows"]
        tn = tot["nw"]
        nc = len(panel_list)
        fout.write("\t".join([
            f"{tier}_TOTAL", f"({nc}_contigs)", str(tot["bp"]), str(tn),
            f"{tot['gc_wt']/tn:.4f}" if tn > 0 else "0",
            f"{tot['wm_wt']/tn:.4f}" if tn > 0 else "0",
            f"{tot['n_at']/tn:.4f}" if tn > 0 else "0",
            f"{tot['n_gc']/tn:.4f}" if tn > 0 else "0",
            f"{tot['n_nr']/tn:.4f}" if tn > 0 else "0",
            str(tot["n_at"]), str(tot["n_gc"]), str(tot["n_nr"]),
        ]) + "\n")

print(f"\nGenerating figure ({total} contigs) ...")

fig_w = 170 / 25.4
fig_h = 225 / 25.4
fig = plt.figure(figsize=(fig_w, fig_h))
fig.patch.set_facecolor('white')

n_pixels_h = 400

# Manual layout: 3 panel rows
# Each panel occupies a vertical band of the figure
# Panel heights proportional to max contigs (equal here since 5,6,6)
panel_tops = [0.93, 0.62, 0.31]
panel_heights = [0.27, 0.27, 0.27]
strip_bottom_margin = 0.02
label_size_margin = 0.015

for pi, (panel_contigs, plabel) in enumerate(zip(panels, panel_labels)):
    n_c = len(panel_contigs)
    if n_c == 0:
        continue

    p_top = panel_tops[pi]
    p_h = panel_heights[pi]
    strip_top = p_top - 0.05   # space for title
    strip_h = p_h - 0.07       # space for title

    # Evenly space contigs across figure width
    left_margin = 0.10
    right_margin = 0.95
    total_w = right_margin - left_margin
    strip_w = total_w / n_c * 0.65
    gap = total_w / n_c * 0.35

    for ci, cinfo in enumerate(panel_contigs):
        cname = cinfo["contig"]
        prof = profiles[cname]

        x_left = left_margin + ci * (strip_w + gap)

        # Banding strip
        ax = fig.add_axes([x_left, strip_top - strip_h, strip_w, strip_h])
        img = render_banding_image(prof["pos"], prof["gc"], prof["wm"],
                                    prof["length"], n_pixels_h)
        ax.imshow(img.reshape(n_pixels_h, 1, 3), aspect='auto',
                  extent=[0, 1, 1, 0], interpolation='nearest')
        ax.add_patch(mpatches.FancyBboxPatch(
            (0.02, 0.0), 0.96, 1.0, boxstyle="round,pad=0.005",
            linewidth=0.4, edgecolor='black', facecolor='none', zorder=10))
        ax.set_xlim(0, 1); ax.set_ylim(1, 0)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

        # Contig name + size above
        size_kb = prof["length"] / 1000
        ax.set_title(f"{cname}  ({size_kb:.0f} kb)", fontsize=5,
                     color='black', fontstyle='italic', pad=4)

    # Panel label — aligned to left margin
    fig.text(left_margin - 0.06, p_top - 0.01, plabel,
             fontsize=13, fontweight='bold', color='black',
             va='top', ha='left')

# Legend — boxed, at bottom with space
ax_leg = fig.add_axes([0.10, 0.01, 0.80, 0.035])
ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
ax_leg.set_facecolor('white')
for s in ax_leg.spines.values():
    s.set_visible(True)
    s.set_color('#cccccc')
    s.set_linewidth(0.4)
ax_leg.set_xticks([]); ax_leg.set_yticks([])

for i, (col, label) in enumerate([
    (COL_AT_REP, 'AT-rich repetitive'),
    (COL_GC_REP, 'GC-rich repetitive'),
    (COL_REP_POOR, 'Non-repetitive'),
]):
    x = 0.03 + i * 0.33
    ec = '#cccccc' if col == COL_REP_POOR else 'none'
    ax_leg.add_patch(plt.Rectangle((x, 0.20), 0.04, 0.60,
                     facecolor=col, edgecolor=ec, linewidth=0.4))
    ax_leg.text(x + 0.06, 0.50, label, fontsize=5.5, va='center',
                color='black')

# Save
basename = "w_genic_banding"
for fmt, kwargs in [
    ("png",  {"dpi": 300, "bbox_inches": "tight", "facecolor": "white"}),
    ("tiff", {"dpi": 300, "bbox_inches": "tight", "facecolor": "white",
              "pil_kwargs": {"compression": "tiff_lzw"}}),
    ("pdf",  {"bbox_inches": "tight", "facecolor": "white"}),
]:
    out = os.path.join(OUT_DIR, f"{basename}.{fmt}")
    fig.savefig(out, **kwargs)
    print(f"  Saved: {out}")

plt.close(fig)
print("\nDone.")
