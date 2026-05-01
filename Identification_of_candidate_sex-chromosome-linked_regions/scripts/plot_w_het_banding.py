#!/usr/bin/env python3
"""
plot_w_het_banding.py
======================
Composite banding of the 60 largest W-heterochromatin (repeat-strong)
contigs.

Same three-colour scheme, WindowMasker intervals, 2 kb / 400 bp
windows, and WM_THRESHOLD=0.4 as all other banding figures.

Author: Kopp K. Pristimantis euphronides genome project.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
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

N_CONTIGS = 60

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


def select_contigs():
    candidates = []
    with open(MASTER_TABLE) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            uc = row.get("unplaced_class", "").strip()
            conf = row.get("unplaced_confidence", "").strip()
            length = int(row["length"])
            if uc == "W-heterochromatin" and conf == "repeat-strong":
                candidates.append({"contig": row["contig"], "length": length})
    candidates.sort(key=lambda x: x["length"], reverse=True)
    return candidates[:N_CONTIGS]


# ============================================================
# FIGURE
# ============================================================

print("Selecting contigs ...")
candidates = select_contigs()
print(f"  {len(candidates)} contigs")

contig_names = [c["contig"] for c in candidates]
print("Extracting ...")
sequences = extract_contigs(ASM_FASTA, contig_names)
print(f"  {len(sequences)} contigs")

print("Loading WindowMasker ...")
wm_intervals = {}
for wf in [WM_PLACED, WM_UNPLACED]:
    wm_intervals.update(parse_windowmasker_intervals(wf))
print(f"  {len(wm_intervals)} contigs")

print("Computing profiles ...")
profiles = {}
for cinfo in candidates:
    cname = cinfo["contig"]
    if cname not in sequences:
        continue
    seq = sequences[cname]
    clen = len(seq)
    if clen < WINDOW:
        continue
    ew = max(200, min(WINDOW, clen // 5))
    es = max(100, min(STEP, ew // 2))
    pos, gc = sliding_window_gc(seq, ew, es)
    wm = None
    if cname in wm_intervals:
        wm = windowmasker_to_density(wm_intervals[cname], clen, pos, ew)
    profiles[cname] = {"pos": pos, "gc": gc, "wm": wm, "length": clen}

plot_contigs = [c for c in candidates if c["contig"] in profiles]
n = len(plot_contigs)

n_cols = 4
n_rows = (n + n_cols - 1) // n_cols
n_pixels_h = 300

fig_w = 170 / 25.4
fig_h = min(225 / 25.4, n_rows * 0.88 + 0.8)

fig = plt.figure(figsize=(fig_w, fig_h))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(n_rows, n_cols, figure=fig,
                       wspace=0.15, hspace=0.55)

for idx, cinfo in enumerate(plot_contigs):
    row = idx // n_cols
    col = idx % n_cols
    cname = cinfo["contig"]
    prof = profiles[cname]

    ax = fig.add_subplot(gs[row, col])
    ax.set_facecolor('white')

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

    # Single-line title: name + size
    size_kb = prof["length"] / 1000
    ax.set_title(f"{cname}  ({size_kb:.0f} kb)", fontsize=5,
                 color='black', fontstyle='italic', pad=3)

# Legend — boxed
ax_leg = fig.add_axes([0.10, 0.003, 0.80, 0.028])
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
    ax_leg.add_patch(plt.Rectangle((x, 0.15), 0.04, 0.70,
                     facecolor=col, edgecolor=ec, linewidth=0.4))
    ax_leg.text(x + 0.06, 0.50, label, fontsize=5.5, va='center',
                color='black')

plt.subplots_adjust(bottom=0.04, top=0.98, left=0.02, right=0.98)

basename = "w_het_banding"
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
