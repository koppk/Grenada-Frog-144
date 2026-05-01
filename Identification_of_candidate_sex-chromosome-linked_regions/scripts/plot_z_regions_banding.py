#!/usr/bin/env python3
"""
plot_z_regions_banding.py
==========================
Composite banding of scaffold_2 and scaffold_8 with Z-candidate
regions marked.

Same three-colour scheme, WindowMasker intervals, 2 kb / 400 bp
windows, and WM_THRESHOLD=0.4 as all other banding figures.

One figure per scaffold showing the full scaffold as a single
vertical strip with the Z-region boundary marked in blue.

Author: Kopp K. Pristimantis euphronides genome project.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import sys
from collections import defaultdict

# ============================================================
# PATHS
# ============================================================

ASM_FASTA = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/GrenadaFrog144/Pristimantis_euphronides.genome.fasta"

BASE = "/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"
WM_PLACED = os.path.join(BASE, "windowmasker", "wm_placed_intervals.txt")
AGP_FILE = "/data/GrenadaFrog144/ragtag.scaffold.renamed.agp"

OUT_DIR = os.path.join(BASE, "gametolog_discovery_hanno7", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

SCAFFOLDS = {
    "scaffold_2": {"length": 176_323_946, "z_start": 0, "z_end": 42_000_000},
    "scaffold_8": {"length": 81_027_173, "z_start": 47_000_000, "z_end": 64_000_000},
}

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


CLASS_NAMES = ["AT-rich_repetitive", "GC-rich_repetitive", "non-repetitive"]


def classify_pixel(rgb):
    """Match a rendered pixel RGB to one of the three classes."""
    if tuple(rgb) == COL_AT_REP:
        return 0
    elif tuple(rgb) == COL_GC_REP:
        return 1
    else:
        return 2


def pixel_segments(img, seq_len, n_pixels, min_segment_mb=1.0):
    """Walk the rendered pixel array, find segments of same class,
    convert to Mb coordinates. Merge segments shorter than min_segment_mb
    into their longer neighbour."""
    classes = np.array([classify_pixel(img[p]) for p in range(n_pixels)])

    # Raw segments
    raw = []
    seg_start = 0
    seg_class = classes[0]
    for p in range(1, n_pixels):
        if classes[p] != seg_class:
            raw.append((seg_start, p, int(seg_class)))
            seg_start = p
            seg_class = classes[p]
    raw.append((seg_start, n_pixels, int(seg_class)))

    # Convert pixel boundaries to bp
    segments = []
    for (ps, pe, cls) in raw:
        bp_s = int(ps / n_pixels * seq_len)
        bp_e = int(pe / n_pixels * seq_len)
        segments.append([bp_s, bp_e, cls])

    # Merge short segments into longer neighbour
    min_bp = int(min_segment_mb * 1e6)
    merged = True
    while merged:
        merged = False
        new_segs = []
        i = 0
        while i < len(segments):
            s = segments[i]
            span = s[1] - s[0]
            if span < min_bp and len(segments) > 1:
                merged = True
                if i > 0 and new_segs:
                    new_segs[-1][1] = s[1]
                elif i + 1 < len(segments):
                    segments[i + 1][0] = s[0]
                else:
                    new_segs.append(s)
            else:
                new_segs.append(s)
            i += 1
        segments = new_segs

    # Compute per-segment pixel-class fractions from original classes array
    result = []
    for (bp_s, bp_e, cls) in segments:
        ps = int(bp_s / seq_len * n_pixels)
        pe = int(bp_e / seq_len * n_pixels)
        pe = min(pe, n_pixels)
        chunk = classes[ps:pe]
        cnt = len(chunk)
        if cnt == 0:
            fat = fgc = fnr = 0.0
        else:
            fat = np.sum(chunk == 0) / cnt
            fgc = np.sum(chunk == 1) / cnt
            fnr = np.sum(chunk == 2) / cnt
        result.append((bp_s, bp_e, cls, fat, fgc, fnr))

    return result


# ============================================================
# AGP + WINDOWMASKER
# ============================================================

def read_agp(agp_path):
    entries = defaultdict(list)
    with open(agp_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9 or cols[4] != 'W':
                continue
            entries[cols[0]].append((
                int(cols[1]) - 1, int(cols[2]), cols[5],
                int(cols[6]) - 1, int(cols[7]), cols[8]))
    return dict(entries)


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


def translate_wm_to_scaffold(agp_entries, wm_contig):
    intervals = []
    for (scf_s, scf_e, contig, ctg_s, ctg_e, orient) in agp_entries:
        if contig not in wm_contig:
            continue
        ctg_len = ctg_e - ctg_s
        for (wm_s, wm_e) in wm_contig[contig]:
            ws = max(wm_s, ctg_s) - ctg_s
            we = min(wm_e, ctg_e) - ctg_s
            if ws >= we:
                continue
            if orient == '-':
                ws, we = ctg_len - we, ctg_len - ws
            intervals.append((scf_s + ws, scf_s + we))
    return intervals


def build_wm_mask(intervals, length):
    mask = np.zeros(length, dtype=np.float32)
    for s, e in intervals:
        mask[max(0, s):min(e + 1, length)] = 1.0
    return mask


def wm_density_from_mask(mask, positions, window_size):
    hw = window_size // 2
    n = len(mask)
    return np.array([np.mean(mask[max(0, p - hw):min(n, p + hw)])
                     for p in positions])


# ============================================================
# SEQUENCE + GC
# ============================================================

def extract_scaffold(fasta_path, scaffold_name):
    seq = []
    in_target = False
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if in_target:
                    break
                in_target = (line[1:].split()[0] == scaffold_name)
                continue
            if in_target:
                seq.append(line.strip().upper())
    return ''.join(seq)


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


# ============================================================
# BIN-AVERAGED RENDERING
# ============================================================

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
# LOAD DATA
# ============================================================

print("Loading AGP ...")
agp = read_agp(AGP_FILE)

print("Loading WindowMasker intervals ...")
wm_contig = parse_windowmasker_intervals(WM_PLACED) \
    if os.path.isfile(WM_PLACED) else {}
print(f"  {len(wm_contig)} contigs")

print("Translating WM to scaffold coordinates ...")
scaffold_wm_masks = {}
for scf_name, scf_info in SCAFFOLDS.items():
    if scf_name not in agp:
        continue
    intervals = translate_wm_to_scaffold(agp[scf_name], wm_contig)
    mask = build_wm_mask(intervals, scf_info["length"])
    scaffold_wm_masks[scf_name] = mask
    print(f"  {scf_name}: masked={np.mean(mask):.2f}")


# ============================================================
# GENERATE
# ============================================================

def generate_scaffold_figure(scaffold_name, info):
    scf_len = info["length"]
    z_start = info["z_start"]
    z_end = info["z_end"]

    print(f"\n  {scaffold_name}: {scf_len / 1e6:.1f} Mb, "
          f"Z {z_start / 1e6:.0f}–{z_end / 1e6:.0f} Mb")

    full_seq = extract_scaffold(ASM_FASTA, scaffold_name)
    print(f"    {len(full_seq)} bp")

    wm_mask = scaffold_wm_masks.get(scaffold_name)

    pos, gc = sliding_window_gc(full_seq, WINDOW, STEP)
    wm = wm_density_from_mask(wm_mask, pos, WINDOW) \
        if wm_mask is not None else None
    print(f"    {len(pos)} windows")

    n_pixels = 1200
    img = render_banding_image(pos, gc, wm, scf_len, n_pixels)

    # ---- Pixel-based segment detection ----
    segments = pixel_segments(img, scf_len, n_pixels, min_segment_mb=1.0)
    tsv_path = os.path.join(OUT_DIR, f"z_banding_{scaffold_name}_segments.tsv")
    with open(tsv_path, "w") as fout:
        hdr = ["scaffold", "seg_start_Mb", "seg_end_Mb", "span_Mb",
               "rendered_class", "in_Z_region",
               "pix_frac_AT_rep", "pix_frac_GC_rep", "pix_frac_nonrep"]
        fout.write("\t".join(hdr) + "\n")
        for (bp_s, bp_e, cls, fat, fgc, fnr) in segments:
            span = (bp_e - bp_s) / 1e6
            mid = (bp_s + bp_e) / 2
            if mid >= z_start and mid < z_end:
                z_label = "Z"
            else:
                z_label = "non-Z"
            fout.write("\t".join([
                scaffold_name,
                f"{bp_s / 1e6:.2f}", f"{bp_e / 1e6:.2f}", f"{span:.2f}",
                CLASS_NAMES[cls], z_label,
                f"{fat:.4f}", f"{fgc:.4f}", f"{fnr:.4f}",
            ]) + "\n")
    print(f"    Segments: {len(segments)} → {tsv_path}")

    # Single strip figure
    fig_w = 170 / 25.4 * 0.30
    fig_h = 7.0

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    fig.patch.set_facecolor('white')

    ax.imshow(img.reshape(n_pixels, 1, 3), aspect='auto',
              extent=[0, 1, 1, 0], interpolation='nearest')

    # Outline
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.02, 0.0), 0.96, 1.0, boxstyle="round,pad=0.005",
        linewidth=0.6, edgecolor='black', facecolor='none', zorder=10))

    # Z boundary — BLUE
    z_s = z_start / scf_len
    z_e = z_end / scf_len
    for y in [z_s, z_e]:
        if 0 < y < 1:
            ax.axhline(y=y, color='#2060C0', linewidth=1.0,
                        linestyle='--', alpha=0.9, zorder=15)
    ax.annotate('', xy=(1.08, z_s), xytext=(1.08, z_e),
                arrowprops=dict(arrowstyle='|-|', color='#2060C0',
                                lw=1.0, mutation_scale=4),
                annotation_clip=False)
    ax.text(1.18, (z_s + z_e) / 2, 'Z', ha='left', va='center',
            fontsize=10, fontweight='bold', color='#2060C0', clip_on=False)

    ax.set_xlim(0, 1); ax.set_ylim(1, 0)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"{scaffold_name} ({scf_len / 1e6:.0f} Mb)",
                 fontsize=9, color='black', pad=6)

    # Legend — boxed, vertical stack, below strip
    ax_leg = fig.add_axes([0.02, 0.005, 0.78, 0.09])
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
        y = 0.75 - i * 0.30
        ec = '#cccccc' if col == COL_REP_POOR else 'none'
        ax_leg.add_patch(plt.Rectangle((0.05, y - 0.08), 0.08, 0.16,
                         facecolor=col, edgecolor=ec, linewidth=0.4))
        ax_leg.text(0.18, y, label, fontsize=5.5, va='center',
                    color='black')

    plt.subplots_adjust(bottom=0.13, top=0.95, left=0.05, right=0.82)

    # Save
    for fmt, kwargs in [
        ("png",  {"dpi": 300, "bbox_inches": "tight", "facecolor": "white"}),
        ("tiff", {"dpi": 300, "bbox_inches": "tight", "facecolor": "white",
                  "pil_kwargs": {"compression": "tiff_lzw"}}),
        ("pdf",  {"bbox_inches": "tight", "facecolor": "white"}),
    ]:
        out = os.path.join(OUT_DIR, f"z_banding_{scaffold_name}.{fmt}")
        fig.savefig(out, **kwargs)
        print(f"  Saved: {out}")
    plt.close(fig)


for scf_name, scf_info in SCAFFOLDS.items():
    generate_scaffold_figure(scf_name, scf_info)

print("\nDone.")
