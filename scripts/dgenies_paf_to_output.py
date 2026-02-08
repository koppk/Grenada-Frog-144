#!/usr/bin/env python3
"""
dgenies_paf_to_output.py - Generate D-Genies-equivalent output files from PAF files.

Replicates the downloadable outputs from the D-Genies web interface:
  1. Dot plot PNG image
  2. Association table (TSV)
  3. No-match queries list
  4. No-match targets list

Designed to handle large queries with many small/repetitive contigs that
cannot load in D-Genies' browser-based visualization.

Note: The PAF file must be pre-sorted by significance externally
(e.g. using awk/sort in bash). This script does NO internal sorting.

Usage:
    python3 dgenies_paf_to_output.py --paf sorted.paf --query-idx query.idx --target-idx target.idx -o output_dir/

    # With query contig reordering (gravity-based, to match target):
    python3 dgenies_paf_to_output.py --paf sorted.paf --query-idx query.idx --target-idx target.idx -o output_dir/ --sort

    # Limit lines read (for huge PAF files):
    python3 dgenies_paf_to_output.py --paf sorted.paf --query-idx query.idx --target-idx target.idx -o output_dir/ --max-lines 200000

Requires: Python 3.7+, matplotlib, numpy
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict, defaultdict
from math import sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import numpy as np


# ---------------------------------------------------------------------------
# Index loading / building
# ---------------------------------------------------------------------------

def load_index(idx_path):
    """Load a D-Genies .idx file.

    Returns:
        name (str): sample name (first line)
        order (list[str]): contig names in order
        contigs (dict[str, int]): contig_name -> length
        reversed_c (dict[str, bool]): contig_name -> is_reversed
        abs_start (dict[str, int]): contig_name -> absolute start position
        total_len (int): total length of all contigs
    """
    order = []
    contigs = OrderedDict()
    reversed_c = {}

    with open(idx_path) as fh:
        name = fh.readline().rstrip("\n")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            cname = parts[0]
            clen = int(parts[1])
            rev = bool(int(parts[2])) if len(parts) > 2 else False
            order.append(cname)
            contigs[cname] = clen
            reversed_c[cname] = rev

    # Compute absolute start positions (cumulative sum)
    abs_start = {}
    running = 0
    for cname in order:
        abs_start[cname] = running
        running += contigs[cname]
    total_len = running

    return name, order, contigs, reversed_c, abs_start, total_len


def build_index_from_paf(paf_path, which="query"):
    """Build an index from PAF columns when no .idx file is provided.

    Args:
        which: "query" to use columns 0,1 or "target" to use columns 5,6
    """
    contigs = OrderedDict()
    col_name = 0 if which == "query" else 5
    col_len = 1 if which == "query" else 6

    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            cname = parts[col_name]
            clen = int(parts[col_len])
            if cname not in contigs:
                contigs[cname] = clen
            else:
                contigs[cname] = max(contigs[cname], clen)

    order = list(contigs.keys())
    reversed_c = {c: False for c in order}
    abs_start = {}
    running = 0
    for c in order:
        abs_start[c] = running
        running += contigs[c]
    total_len = running
    name = which
    return name, order, contigs, reversed_c, abs_start, total_len


# ---------------------------------------------------------------------------
# Small-contig merging (###MIX### groups, matching D-Genies behavior)
# ---------------------------------------------------------------------------

MIN_CONTIG_FRAC = 0.002  # 0.2% of total length
MIN_CONSECUTIVE_SMALL = 5

def merge_small_contigs(order, contigs, abs_start, total_len):
    """Merge small consecutive contigs into ###MIX### groups.

    Returns new (order, contigs, abs_start) with merged groups.
    """
    threshold = total_len * MIN_CONTIG_FRAC
    new_order = []
    new_contigs = OrderedDict()
    new_abs_start = {}
    # Map original contig names to their (possibly merged) name
    contig_rename = {}

    i = 0
    mix_id = 0
    while i < len(order):
        c = order[i]
        if contigs[c] < threshold:
            # Collect consecutive small contigs
            group = []
            while i < len(order) and contigs[order[i]] < threshold:
                group.append(order[i])
                i += 1
            if len(group) >= MIN_CONSECUTIVE_SMALL:
                mix_name = f"###MIX###{mix_id}"
                mix_id += 1
                group_len = sum(contigs[g] for g in group)
                group_start = abs_start[group[0]]
                new_order.append(mix_name)
                new_contigs[mix_name] = group_len
                new_abs_start[mix_name] = group_start
                for g in group:
                    contig_rename[g] = mix_name
            else:
                # Not enough consecutive smalls -> keep individually
                for g in group:
                    new_order.append(g)
                    new_contigs[g] = contigs[g]
                    new_abs_start[g] = abs_start[g]
                    contig_rename[g] = g
        else:
            new_order.append(c)
            new_contigs[c] = contigs[c]
            new_abs_start[c] = abs_start[c]
            contig_rename[c] = c
            i += 1

    return new_order, new_contigs, new_abs_start, contig_rename


def rebuild_full_abs_start(merged_order, merged_abs_start, orig_order,
                           orig_contigs, contig_rename):
    """Rebuild abs_start for ALL original contig names from a merged/sorted order.

    After merging small contigs into ###MIX### groups and/or sorting,
    the abs_start dict only has group names (not individual member names).
    This function computes the correct absolute position for every original
    contig by determining its offset within its MIX group.

    For non-merged contigs, the position is taken directly from merged_abs_start.
    """
    # Build MIX group membership preserving original contig order within groups
    mix_members = defaultdict(list)
    for name in orig_order:
        renamed = contig_rename.get(name, name)
        if renamed.startswith("###MIX###"):
            mix_members[renamed].append(name)

    full = {}
    for name in merged_order:
        if name in mix_members:
            offset = merged_abs_start[name]
            for member in mix_members[name]:
                full[member] = offset
                offset += orig_contigs[member]
        else:
            full[name] = merged_abs_start[name]
    return full



# ---------------------------------------------------------------------------
# PAF parsing
# ---------------------------------------------------------------------------

def classify_identity(idy):
    """Classify identity into D-Genies identity class."""
    if idy < 0.25:
        return 0
    elif idy < 0.50:
        return 1
    elif idy < 0.75:
        return 2
    else:
        return 3


# D-Genies color scheme
IDY_COLORS = {
    3: "#094b09",  # dark green  (>= 75%)
    2: "#2ebd40",  # green       (50-75%)
    1: "#d5670b",  # orange      (25-50%)
    0: "#ffd84b",  # yellow      (< 25%)
}

IDY_LABELS = {
    3: "Identity >= 75%",
    2: "50% <= Identity < 75%",
    1: "25% <= Identity < 50%",
    0: "Identity < 25%",
}


def parse_paf(paf_path, q_abs_start, t_abs_start, q_contigs, t_contigs,
              q_contig_rename=None, t_contig_rename=None,
              q_abs_start_orig=None, t_abs_start_orig=None,
              max_lines=None):
    """Parse a PAF file and return match data.

    The PAF is expected to already be sorted by significance externally
    (e.g. using bash/awk/samtools). No internal sorting is performed.

    Args:
        q_abs_start, t_abs_start: merged abs_start dicts (for display grouping)
        q_abs_start_orig, t_abs_start_orig: original (pre-merge) abs_start
            dicts for correct coordinate computation. If not provided, falls
            back to q_abs_start/t_abs_start.

    Returns:
        matches: list of dicts with keys x1, x2, y1, y2, idy, idy_class,
                 q_name, t_name, length
        q_names_seen: set of query contig names with at least one match
        t_names_seen: set of target contig names with at least one match
        sampled: bool - True if max_lines was hit
    """
    if q_contig_rename is None:
        q_contig_rename = {}
    if t_contig_rename is None:
        t_contig_rename = {}

    # Use original abs_start for coordinate lookups (critical when contigs
    # are merged into ###MIX### groups - the merged dict only has group
    # starts, not individual contig positions)
    q_abs_lookup = q_abs_start_orig if q_abs_start_orig else q_abs_start
    t_abs_lookup = t_abs_start_orig if t_abs_start_orig else t_abs_start

    raw_lines = []
    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            raw_lines.append(line)

    sampled = False
    if max_lines and len(raw_lines) > max_lines:
        raw_lines = raw_lines[:max_lines]
        sampled = True

    matches = []
    q_names_seen = set()
    t_names_seen = set()

    for line in raw_lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue

        q_name_orig = parts[0]
        t_name_orig = parts[5]

        # Skip if contig not in original index
        if q_name_orig not in q_abs_lookup or t_name_orig not in t_abs_lookup:
            continue

        strand = 1 if parts[4] == "+" else -1
        q_start = int(parts[2])
        q_end = int(parts[3])
        t_start = int(parts[7])
        t_end = int(parts[8])
        n_matches = int(parts[9])
        block_len = int(parts[10])

        if block_len == 0:
            continue

        idy = n_matches / block_len

        # Compute absolute coordinates using ORIGINAL contig positions
        # (not merged MIX group positions)
        # x = target, y = query (D-Genies convention)
        q_abs = q_abs_lookup[q_name_orig]
        t_abs = t_abs_lookup[t_name_orig]

        if strand == 1:
            x1 = t_start + t_abs
            x2 = t_end + t_abs
            y1 = q_start + q_abs
            y2 = q_end + q_abs
        else:
            x1 = t_end + t_abs
            x2 = t_start + t_abs
            y1 = q_start + q_abs
            y2 = q_end + q_abs

        euclidean = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        matches.append({
            "x1": x1, "x2": x2,
            "y1": y1, "y2": y2,
            "idy": idy,
            "idy_class": classify_identity(idy),
            "q_name": q_name_orig,
            "t_name": t_name_orig,
            "length": euclidean,
            "strand": strand,
            "q_start": q_start,
            "q_end": q_end,
            "t_start": t_start,
            "t_end": t_end,
        })
        q_names_seen.add(q_name_orig)
        t_names_seen.add(t_name_orig)

    return matches, q_names_seen, t_names_seen, sampled


# ---------------------------------------------------------------------------
# Noise removal (matches D-Genies remove_noise)
# ---------------------------------------------------------------------------

def remove_noise(matches):
    """Remove small noisy matches using histogram-based threshold."""
    if not matches:
        return matches

    lengths = [m["length"] for m in matches]
    n = len(lengths)
    if n < 20:
        return matches

    n_bins = max(10, n // 10)
    hist, bin_edges = np.histogram(lengths, bins=n_bins)

    max_idx = np.argmax(hist)
    max_count = hist[max_idx]
    threshold_count = max_count / 50.0

    # Scan forward from max to find where count drops below threshold
    cutoff_idx = max_idx
    for i in range(max_idx + 1, len(hist)):
        if hist[i] < threshold_count:
            cutoff_idx = i
            break
    else:
        # Never dropped below -> no noise removal
        return matches

    noise_limit = bin_edges[cutoff_idx]
    filtered = [m for m in matches if m["length"] >= noise_limit]
    removed = len(matches) - len(filtered)
    if removed > 0:
        print(f"  Noise removal: removed {removed} matches below length threshold {noise_limit:.0f}")
    return filtered


# ---------------------------------------------------------------------------
# Contig sorting (sort query contigs to best match target order)
# ---------------------------------------------------------------------------

def sort_query_contigs(matches, q_order, q_contigs, q_abs_start, t_order,
                       t_contigs, t_abs_start, t_total_len,
                       t_abs_start_orig=None):
    """Sort query contigs to match target chromosome order using gravity.

    Args:
        t_abs_start_orig: original (pre-merge) target abs_start dict for
                          looking up individual contig positions. Falls back
                          to t_abs_start if not provided.

    Returns new (q_order, q_contigs, q_abs_start, reversed_contigs).
    """
    # Use original target abs_start for position lookups (handles merged contigs)
    t_abs_lookup = t_abs_start_orig if t_abs_start_orig else t_abs_start

    # Compute gravity for each (query_contig, target_contig) pair
    gravity = defaultdict(lambda: defaultdict(float))
    for m in matches:
        q_name = m["q_name"]
        t_name = m["t_name"]
        ln = m["length"]
        gravity[q_name][t_name] += (1 + ln) ** 2

    # For each query contig, find best target and compute gravity center
    contig_target = {}
    contig_gravity_center = {}

    for q_name in q_order:
        if q_name not in gravity:
            # No matches - put at end
            contig_target[q_name] = None
            contig_gravity_center[q_name] = t_total_len
            continue

        # Best target = highest gravity sum
        best_t = max(gravity[q_name], key=gravity[q_name].get)
        contig_target[q_name] = best_t

        # Compute gravity center: weighted median position on target
        weighted_pos = 0.0
        weight_sum = 0.0
        for m in matches:
            if m["q_name"] == q_name and m["t_name"] == best_t:
                t_offset = t_abs_lookup.get(best_t, 0)
                median_t = (m["t_start"] + m["t_end"]) / 2.0 + t_offset
                w = (1 + m["length"]) ** 2
                weighted_pos += median_t * w
                weight_sum += w

        if weight_sum > 0:
            contig_gravity_center[q_name] = weighted_pos / weight_sum
        else:
            contig_gravity_center[q_name] = t_total_len

    # Sort by gravity center
    new_q_order = sorted(q_order, key=lambda c: contig_gravity_center[c])

    # Check orientation for each contig
    reversed_contigs = {}
    for q_name in new_q_order:
        # Check if contig should be reversed by examining match slopes
        plus_weight = 0.0
        minus_weight = 0.0
        for m in matches:
            if m["q_name"] == q_name:
                w = m["length"] ** 2
                if m["strand"] == 1:
                    plus_weight += w
                else:
                    minus_weight += w
        reversed_contigs[q_name] = minus_weight > plus_weight

    # Rebuild abs_start with new order
    new_abs_start = {}
    running = 0
    for c in new_q_order:
        new_abs_start[c] = running
        running += q_contigs[c]

    return new_q_order, q_contigs, new_abs_start, reversed_contigs


# ---------------------------------------------------------------------------
# Output: Association table
# ---------------------------------------------------------------------------

def build_association_table(matches, q_order, q_contigs, t_contigs,
                            q_abs_start, t_abs_start, q_reversed):
    """Build query-to-target association table (matches D-Genies format).

    Returns list of rows: [Query, Target, Strand, Q-len, Q-start, Q-stop,
                           T-len, T-start, T-stop]
    """
    # Compute gravity per (query, target) to find best match
    gravity = defaultdict(lambda: defaultdict(float))
    q_on_t_coords = defaultdict(lambda: defaultdict(lambda: {"q_min": float("inf"),
                                                              "q_max": 0,
                                                              "t_min": float("inf"),
                                                              "t_max": 0}))

    for m in matches:
        q_name = m["q_name"]
        t_name = m["t_name"]
        ln = m["length"]
        gravity[q_name][t_name] += (1 + ln) ** 2
        coords = q_on_t_coords[q_name][t_name]
        coords["q_min"] = min(coords["q_min"], m["q_start"])
        coords["q_max"] = max(coords["q_max"], m["q_end"])
        coords["t_min"] = min(coords["t_min"], m["t_start"])
        coords["t_max"] = max(coords["t_max"], m["t_end"])

    rows = []
    for q_name in q_order:
        if q_name.startswith("###MIX###"):
            continue
        if q_name in gravity:
            best_t = max(gravity[q_name], key=gravity[q_name].get)
            coords = q_on_t_coords[q_name][best_t]
            strand = "-" if q_reversed.get(q_name, False) else "+"
            rows.append([
                q_name, best_t, strand,
                q_contigs[q_name],
                coords["q_min"], coords["q_max"],
                t_contigs[best_t],
                coords["t_min"], coords["t_max"],
            ])
        else:
            strand = "-" if q_reversed.get(q_name, False) else "+"
            rows.append([
                q_name, "None", strand,
                q_contigs[q_name],
                "na", "na", "na", "na", "na",
            ])
    return rows


def write_association_table(rows, out_path):
    """Write association table TSV file."""
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Query", "Target", "Strand", "Q-len", "Q-start",
                         "Q-stop", "T-len", "T-start", "T-stop"])
        for row in rows:
            writer.writerow(row)
    print(f"  Association table written: {out_path}")


# ---------------------------------------------------------------------------
# Output: No-match lists
# ---------------------------------------------------------------------------

def build_no_match_list(all_names, names_seen):
    """Return list of contig names with zero matches."""
    return [n for n in all_names if n not in names_seen and not n.startswith("###MIX###")]


def write_no_match_list(names, out_path, label):
    """Write no-match list to file."""
    with open(out_path, "w") as fh:
        for n in names:
            fh.write(n + "\n")
    print(f"  {label}: {len(names)} entries -> {out_path}")


# ---------------------------------------------------------------------------
# Output: Dot plot PNG
# ---------------------------------------------------------------------------

def _format_bp(val, total_len):
    """Format a base-pair value for axis tick labels."""
    if total_len >= 1e9:
        return f"{val/1e6:.0f} M"
    elif total_len >= 1e6:
        return f"{val/1e6:.1f} M"
    elif total_len >= 1e3:
        return f"{val/1e3:.0f} K"
    return str(int(val))


def draw_dotplot(matches, q_order, q_contigs, q_abs_start, q_total_len,
                 t_order, t_contigs, t_abs_start, t_total_len,
                 q_name, t_name, out_path, dpi=100, figsize=None,
                 sampled=False):
    """Draw an exact D-Genies-style dot plot and save as PNG.

    Replicates D-Genies rendering precisely:
      - 5000x5000 pixel output (50 inches at 100 dpi)
      - Internal 1000-unit coordinate space
      - 90% plot area, 5% axes on each side
      - Target contig names on TOP, Mbp scale on BOTTOM
      - Query contig names on RIGHT, Mbp scale on LEFT
      - Origin at bottom-left
      - Exact D-Genies line widths, colors, dash patterns
    """
    # D-Genies renders at 5000x5000 pixels.
    # We use 50x50 inches at 100 dpi = 5000x5000 px.
    if figsize is None:
        figsize = (50, 50)

    # D-Genies layout: 90% plot area, 5% margin on each side for axes
    # In the 50-inch figure, plot area = 45 inches, axes = 2.5 inches each side
    plot_frac = 0.90
    axis_frac = 0.05

    fig = plt.figure(figsize=figsize, facecolor="white")

    # Create main plot axes at the correct position (5% to 95% in both dims)
    ax = fig.add_axes([axis_frac, axis_frac, plot_frac, plot_frac])
    ax.set_facecolor("white")

    # D-Genies internal scale = 1000 units
    SCALE = 1000.0

    # --- D-Genies exact rendering parameters ---
    # Match line width: scale/400 = 2.5 in 1000-unit space
    # In figure coords: 2.5/1000 of the plot area width in inches, converted to points
    plot_inches = figsize[0] * plot_frac  # 45 inches
    plot_pixels = plot_inches * dpi       # 4500 pixels
    # 2.5 units out of 1000 = 0.25% of plot area
    lw_points = (2.5 / SCALE) * plot_inches * 72.0 / 1.0  # in points (72 pt/inch)
    # That gives ~8.1 points which is very thick. D-Genies SVG uses viewBox scaling
    # so the effective line width depends on the viewBox-to-pixel ratio.
    # SVG viewBox 0-1000 mapped to 90% of 5000px = 4500px.
    # stroke-width=2.5 in SVG coords means 2.5/1000 * 4500 = 11.25 pixels
    # In matplotlib points: 11.25 pixels / (dpi/72) = 11.25 / (100/72) = 8.1 points
    lw = 11.25 / (dpi / 72.0)

    # Break line width: scale/1500 = 0.667 in 1000-unit space
    # 0.667/1000 * 4500 = 3.0 pixels -> points
    break_lw = 3.0 / (dpi / 72.0)

    # Break line color and dash pattern
    break_color = "#7c7c7c"
    # D-Genies dash "3, 3" in 1000-unit space -> 3/1000 * 4500 = 13.5 pixels
    break_dash_px = 13.5 / (dpi / 72.0)  # in points
    break_dash = (break_dash_px, break_dash_px)

    # Axis background color
    axis_bg_color = "#f4f4f4"

    # Mix zone color
    mix_color = "#969696"

    # --- Draw axis background (the light gray border bands) ---
    # D-Genies draws trapezoid backgrounds. We approximate with filled rectangles
    # covering the axis areas.
    # Top axis background
    ax_top_bg = fig.add_axes([axis_frac, axis_frac + plot_frac, plot_frac, axis_frac])
    ax_top_bg.set_facecolor(axis_bg_color)
    ax_top_bg.set_xlim(0, 1); ax_top_bg.set_ylim(0, 1)
    ax_top_bg.set_xticks([]); ax_top_bg.set_yticks([])
    for spine in ax_top_bg.spines.values():
        spine.set_visible(False)

    # Right axis background
    ax_right_bg = fig.add_axes([axis_frac + plot_frac, axis_frac, axis_frac, plot_frac])
    ax_right_bg.set_facecolor(axis_bg_color)
    ax_right_bg.set_xlim(0, 1); ax_right_bg.set_ylim(0, 1)
    ax_right_bg.set_xticks([]); ax_right_bg.set_yticks([])
    for spine in ax_right_bg.spines.values():
        spine.set_visible(False)

    # Bottom axis background
    ax_bottom_bg = fig.add_axes([axis_frac, 0, plot_frac, axis_frac])
    ax_bottom_bg.set_facecolor(axis_bg_color)
    ax_bottom_bg.set_xlim(0, 1); ax_bottom_bg.set_ylim(0, 1)
    ax_bottom_bg.set_xticks([]); ax_bottom_bg.set_yticks([])
    for spine in ax_bottom_bg.spines.values():
        spine.set_visible(False)

    # Left axis background
    ax_left_bg = fig.add_axes([0, axis_frac, axis_frac, plot_frac])
    ax_left_bg.set_facecolor(axis_bg_color)
    ax_left_bg.set_xlim(0, 1); ax_left_bg.set_ylim(0, 1)
    ax_left_bg.set_xticks([]); ax_left_bg.set_yticks([])
    for spine in ax_left_bg.spines.values():
        spine.set_visible(False)

    # --- Draw alignment matches as line segments, grouped by identity class ---
    # D-Genies draws class 0 first (behind), then 1, 2, 3 (on top)
    # Within each class, sorted by identity ascending (low on bottom)
    for idy_class in [0, 1, 2, 3]:
        class_matches = [m for m in matches if m["idy_class"] == idy_class]
        if not class_matches:
            continue

        # Sort within class by identity ascending (lower drawn first)
        class_matches.sort(key=lambda m: m["idy"])

        segments = []
        for m in class_matches:
            # Scale to 0..SCALE (1000) coordinate space
            x1_s = m["x1"] / t_total_len * SCALE
            x2_s = m["x2"] / t_total_len * SCALE
            # Y-axis inverted: D-Genies does y = SCALE - (genomic_y / total * SCALE)
            y1_s = SCALE - (m["y1"] / q_total_len * SCALE)
            y2_s = SCALE - (m["y2"] / q_total_len * SCALE)

            segments.append([(x1_s, y1_s), (x2_s, y2_s)])

        lc = LineCollection(segments, colors=IDY_COLORS[idy_class],
                            linewidths=lw, capstyle="round",
                            zorder=2 + idy_class)
        ax.add_collection(lc)

    # --- Draw break lines (contig boundaries) ---
    # D-Genies: dashed "3,3", color #7c7c7c, width scale/1500
    # Drawn for all contigs except the last one
    for i, cname in enumerate(t_order):
        if i == len(t_order) - 1:
            break
        x = (t_abs_start[cname] + t_contigs[cname]) / t_total_len * SCALE
        ax.plot([x, x], [0, SCALE], color=break_color,
                linewidth=break_lw, linestyle=(0, break_dash),
                zorder=1)

    for i, cname in enumerate(q_order):
        if i == len(q_order) - 1:
            break
        # Y inverted
        y = SCALE - (q_abs_start[cname] + q_contigs[cname]) / q_total_len * SCALE
        ax.plot([0, SCALE], [y, y], color=break_color,
                linewidth=break_lw, linestyle=(0, break_dash),
                zorder=1)

    # --- Set up main plot area ---
    ax.set_xlim(0, SCALE)
    ax.set_ylim(0, SCALE)
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("black")

    # Remove default ticks from main plot
    ax.set_xticks([])
    ax.set_yticks([])

    # --- Bottom axis: Mbp scale (9 tick marks at 10%..90%) ---
    # D-Genies uses 9 ticks at positions 10%, 20%, ..., 90%
    ax_bottom = fig.add_axes([axis_frac, 0, plot_frac, axis_frac])
    ax_bottom.set_xlim(0, 500)
    ax_bottom.set_ylim(0, 20)
    ax_bottom.set_facecolor(axis_bg_color)
    for spine in ax_bottom.spines.values():
        spine.set_visible(False)

    tick_font_size = 6.5 * (figsize[0] / 12.0)  # scale font with figure size
    for i in range(1, 10):
        tick_x = 500.0 / 10 * i
        # Tick line from y=0 to y=5
        ax_bottom.plot([tick_x, tick_x], [0, 5], color="black", linewidth=0.5)
        # Label
        bp_val = t_total_len / 10 * i
        label = _format_bp(bp_val, t_total_len)
        ax_bottom.text(tick_x, 15, label, ha="center", va="center",
                       fontsize=tick_font_size, fontfamily="sans-serif")

    ax_bottom.set_xticks([])
    ax_bottom.set_yticks([])

    # --- Left axis: Mbp scale (9 tick marks, rotated) ---
    ax_left = fig.add_axes([0, axis_frac, axis_frac, plot_frac])
    ax_left.set_xlim(0, 20)
    ax_left.set_ylim(0, 500)
    ax_left.set_facecolor(axis_bg_color)
    for spine in ax_left.spines.values():
        spine.set_visible(False)

    for i in range(1, 10):
        # D-Genies left axis: ticks go from left edge, labels at tick positions
        # Y is inverted in the plot, so genomic position increases downward in SVG
        # but we want it to match: bottom = 0, top = max
        tick_y = 500.0 / 10 * i
        # Tick line from x=15 to x=20 (right edge)
        ax_left.plot([15, 20], [tick_y, tick_y], color="black", linewidth=0.5)
        # Label - D-Genies shows values increasing from bottom to top
        bp_val = q_total_len / 10 * i
        label = _format_bp(bp_val, q_total_len)
        ax_left.text(12, tick_y, label, ha="center", va="center",
                     fontsize=tick_font_size, fontfamily="sans-serif",
                     rotation=90)

    ax_left.set_xticks([])
    ax_left.set_yticks([])

    # --- Top axis: target contig/chromosome names ---
    ax_top = fig.add_axes([axis_frac, axis_frac + plot_frac, plot_frac, axis_frac])
    ax_top.set_xlim(0, 500)
    ax_top.set_ylim(0, 20)
    ax_top.set_facecolor(axis_bg_color)
    for spine in ax_top.spines.values():
        spine.set_visible(False)

    name_font_size = 6.0 * (figsize[0] / 12.0)
    title_font_size = 6.0 * (figsize[0] / 12.0)

    # Genome name (title) at top center, italic
    ax_top.text(250, 7.5, t_name, ha="center", va="center",
                fontsize=title_font_size, fontfamily="sans-serif",
                fontstyle="italic")

    # Contig names and separators
    max_labels = 80
    non_mix_t = [(cname, t_abs_start[cname], t_contigs[cname])
                 for cname in t_order if not cname.startswith("###MIX###")]
    mix_t = [(cname, t_abs_start[cname], t_contigs[cname])
             for cname in t_order if cname.startswith("###MIX###")]

    # Draw MIX zone rectangles
    for cname, start, length in mix_t:
        x0 = start / t_total_len * 500
        x1 = (start + length) / t_total_len * 500
        ax_top.fill_between([x0, x1], 12, 20, color=mix_color, linewidth=0)

    # Draw contig zone separators and names
    running_pos = 0
    for idx, cname in enumerate(t_order):
        start = t_abs_start[cname]
        length = t_contigs[cname]
        end = start + length

        # Zone separator (vertical line between zones)
        if idx > 0:
            sep_x = start / t_total_len * 500
            if not cname.startswith("###MIX###"):
                ax_top.plot([sep_x, sep_x], [12, 20], color="black", linewidth=0.5)

        # Contig name (skip MIX zones, skip if too many)
        if not cname.startswith("###MIX###") and len(non_mix_t) <= max_labels:
            center_x = (start + length / 2) / t_total_len * 500
            # Truncate long names
            display_name = cname
            zone_width = length / t_total_len * 500
            if len(display_name) > 3 and zone_width < len(display_name) * 1.5:
                max_chars = max(3, int(zone_width / 1.5))
                if len(display_name) > max_chars:
                    display_name = display_name[:max_chars-3] + "..."
            ax_top.text(center_x, 17, display_name, ha="center", va="center",
                        fontsize=name_font_size, fontfamily="sans-serif",
                        clip_on=True)

    # --- Right axis: query contig/chromosome names ---
    ax_right = fig.add_axes([axis_frac + plot_frac, axis_frac, axis_frac, plot_frac])
    ax_right.set_xlim(0, 20)
    ax_right.set_ylim(0, 500)
    ax_right.set_facecolor(axis_bg_color)
    for spine in ax_right.spines.values():
        spine.set_visible(False)

    # Genome name (title) at right center, italic, rotated 90 degrees
    ax_right.text(7.5, 250, q_name, ha="center", va="center",
                  fontsize=title_font_size, fontfamily="sans-serif",
                  fontstyle="italic", rotation=90)

    # Query contig names and separators
    non_mix_q = [(cname, q_abs_start[cname], q_contigs[cname])
                 for cname in q_order if not cname.startswith("###MIX###")]
    mix_q = [(cname, q_abs_start[cname], q_contigs[cname])
             for cname in q_order if cname.startswith("###MIX###")]

    # Draw MIX zone rectangles
    for cname, start, length in mix_q:
        # Y inverted to match plot
        y0 = (1.0 - (start + length) / q_total_len) * 500
        y1 = (1.0 - start / q_total_len) * 500
        ax_right.fill_between([0, 8], y0, y1, color=mix_color, linewidth=0)

    # Draw contig zone separators and names
    for idx, cname in enumerate(q_order):
        start = q_abs_start[cname]
        length = q_contigs[cname]

        # Zone separator (horizontal line between zones)
        if idx > 0:
            sep_y = (1.0 - start / q_total_len) * 500
            if not cname.startswith("###MIX###"):
                ax_right.plot([0, 8], [sep_y, sep_y], color="black", linewidth=0.5)

        # Contig name (skip MIX zones, skip if too many)
        if not cname.startswith("###MIX###") and len(non_mix_q) <= max_labels:
            center_y = (1.0 - (start + length / 2) / q_total_len) * 500
            display_name = cname
            zone_height = length / q_total_len * 500
            if len(display_name) > 3 and zone_height < len(display_name) * 1.5:
                max_chars = max(3, int(zone_height / 1.5))
                if len(display_name) > max_chars:
                    display_name = display_name[:max_chars-3] + "..."
            ax_right.text(17, center_y, display_name, ha="center", va="center",
                          fontsize=name_font_size, fontfamily="sans-serif",
                          rotation=90, clip_on=True)

    ax_right.set_xticks([])
    ax_right.set_yticks([])

    # --- Corner backgrounds (fill the 4 corners with axis bg color) ---
    for corner_pos in [(0, 0, axis_frac, axis_frac),
                       (0, axis_frac + plot_frac, axis_frac, axis_frac),
                       (axis_frac + plot_frac, 0, axis_frac, axis_frac),
                       (axis_frac + plot_frac, axis_frac + plot_frac, axis_frac, axis_frac)]:
        ax_corner = fig.add_axes(corner_pos)
        ax_corner.set_facecolor(axis_bg_color)
        ax_corner.set_xticks([]); ax_corner.set_yticks([])
        for spine in ax_corner.spines.values():
            spine.set_visible(False)

    # --- Save ---
    fig.savefig(out_path, dpi=dpi, facecolor="white", pad_inches=0)
    plt.close(fig)
    print(f"  Dot plot written: {out_path} ({int(figsize[0]*dpi)}x{int(figsize[1]*dpi)} px)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate D-Genies-equivalent output files from PAF + index files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--paf", required=True,
                        help="Input PAF file (e.g. map.paf from D-Genies)")
    parser.add_argument("--query-idx", default=None,
                        help="Query .idx file. If omitted, index is built from PAF.")
    parser.add_argument("--target-idx", default=None,
                        help="Target .idx file. If omitted, index is built from PAF.")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("--prefix", default=None,
                        help="Filename prefix for outputs (default: derived from PAF filename)")
    parser.add_argument("--max-lines", type=int, default=None,
                        help="Max PAF lines to process (default: all). "
                             "Lines are sorted by significance first, so the "
                             "most important matches are kept.")
    parser.add_argument("--remove-noise", action="store_true",
                        help="Remove small noisy matches using histogram-based threshold")
    parser.add_argument("--sort", action="store_true",
                        help="Sort query contigs to match target chromosome order")
    parser.add_argument("--no-merge-small", action="store_true",
                        help="Do NOT merge small contigs into ###MIX### groups")
    parser.add_argument("--dpi", type=int, default=100,
                        help="DPI for PNG output (default: 100, gives 5000x5000px)")
    parser.add_argument("--figsize", type=float, nargs=2, default=None,
                        metavar=("WIDTH", "HEIGHT"),
                        help="Figure size in inches (default: 50 50, gives 5000x5000px at 100dpi)")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip dot plot generation (only produce text files)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- Determine prefix ---
    if args.prefix:
        prefix = args.prefix
    else:
        prefix = os.path.splitext(os.path.basename(args.paf))[0]

    print(f"Processing PAF: {args.paf}")

    # --- Load or build indexes ---
    if args.query_idx:
        print(f"Loading query index: {args.query_idx}")
        q_name, q_order, q_contigs, q_reversed, q_abs_start, q_total = load_index(args.query_idx)
    else:
        print("Building query index from PAF...")
        q_name, q_order, q_contigs, q_reversed, q_abs_start, q_total = build_index_from_paf(args.paf, "query")
    print(f"  Query: {q_name} - {len(q_order)} contigs, total {q_total:,} bp")

    if args.target_idx:
        print(f"Loading target index: {args.target_idx}")
        t_name, t_order, t_contigs, t_reversed, t_abs_start, t_total = load_index(args.target_idx)
    else:
        print("Building target index from PAF...")
        t_name, t_order, t_contigs, t_reversed, t_abs_start, t_total = build_index_from_paf(args.paf, "target")
    print(f"  Target: {t_name} - {len(t_order)} contigs, total {t_total:,} bp")

    # --- Preserve original indexes ---
    q_order_orig = list(q_order)
    t_order_orig = list(t_order)
    q_abs_start_orig = dict(q_abs_start)
    t_abs_start_orig = dict(t_abs_start)

    # --- Parse PAF (first pass, using original indexes for sorting/stats) ---
    print("Parsing PAF file...")
    matches, q_seen, t_seen, sampled = parse_paf(
        args.paf, q_abs_start, t_abs_start, q_contigs, t_contigs,
        q_abs_start_orig=q_abs_start_orig, t_abs_start_orig=t_abs_start_orig,
        max_lines=args.max_lines)
    print(f"  Parsed {len(matches)} alignment matches")
    if sampled:
        print(f"  (sampled to {args.max_lines} lines)")

    # --- Optional noise removal ---
    if args.remove_noise:
        matches = remove_noise(matches)
        print(f"  After noise removal: {len(matches)} matches")

    # --- Optional sorting (BEFORE merging, so gravity works on real names) ---
    if args.sort:
        print("Sorting query contigs to match target order...")
        q_order, _, q_abs_start, q_reversed_new = sort_query_contigs(
            matches, q_order, q_contigs, q_abs_start,
            t_order, t_contigs, t_abs_start, t_total,
            t_abs_start_orig=t_abs_start_orig)
        q_reversed.update(q_reversed_new)
        q_abs_start_orig = dict(q_abs_start)
        q_order_orig = list(q_order)  # update orig to sorted order for merge

    # --- Merge small contigs (after sorting) ---
    q_contig_rename = {}
    t_contig_rename = {}
    if not args.no_merge_small:
        q_order, q_contigs_m, q_abs_start_merged, q_contig_rename = merge_small_contigs(
            q_order, q_contigs, q_abs_start, q_total)
        t_order, t_contigs_m, t_abs_start_merged, t_contig_rename = merge_small_contigs(
            t_order, t_contigs, t_abs_start, t_total)
        n_q_merged = sum(1 for v in q_contig_rename.values() if v.startswith("###MIX###"))
        n_t_merged = sum(1 for v in t_contig_rename.values() if v.startswith("###MIX###"))
        if n_q_merged > 0:
            print(f"  Merged {n_q_merged} small query contigs into MIX groups")
        if n_t_merged > 0:
            print(f"  Merged {n_t_merged} small target contigs into MIX groups")
        q_contigs_plot = q_contigs_m
        t_contigs_plot = t_contigs_m
        q_abs_start_plot = q_abs_start_merged
        t_abs_start_plot = t_abs_start_merged
    else:
        q_contigs_plot = q_contigs
        t_contigs_plot = t_contigs
        q_abs_start_plot = q_abs_start
        t_abs_start_plot = t_abs_start

    # --- Re-parse PAF if sorted (with new coordinates) ---
    if args.sort:
        print("Re-parsing PAF with sorted coordinates...")
        matches, q_seen, t_seen, sampled = parse_paf(
            args.paf, q_abs_start_plot, t_abs_start_plot,
            q_contigs_plot, t_contigs_plot,
            q_contig_rename, t_contig_rename,
            q_abs_start_orig=q_abs_start_orig, t_abs_start_orig=t_abs_start_orig,
            max_lines=args.max_lines)
        if args.remove_noise:
            matches = remove_noise(matches)
        print(f"  {len(matches)} matches after re-parse")

    # --- Generate outputs ---
    print("\nGenerating outputs...")

    # 1. Dot plot PNG
    if not args.no_plot:
        plot_path = os.path.join(args.outdir,
                                 f"map_{q_name}_to_{t_name}.png".replace(" ", "_"))
        draw_dotplot(matches, q_order, q_contigs_plot, q_abs_start_plot, q_total,
                     t_order, t_contigs_plot, t_abs_start_plot, t_total,
                     q_name, t_name, plot_path, dpi=args.dpi,
                     figsize=tuple(args.figsize) if args.figsize else None,
                     sampled=sampled)

    # 2. Association table (use original contigs, not merged)
    assoc_path = os.path.join(args.outdir,
                              f"{q_name}_{t_name}_assoc.tsv".replace(" ", "_"))
    assoc_rows = build_association_table(
        matches, list(q_contigs.keys()), q_contigs, t_contigs,
        q_abs_start_orig, t_abs_start_orig, q_reversed)
    write_association_table(assoc_rows, assoc_path)

    # 3. No-match queries
    no_q_path = os.path.join(args.outdir,
                             f"no_query_matches_{q_name}_to_{t_name}.txt".replace(" ", "_"))
    no_q = build_no_match_list(list(q_contigs.keys()), q_seen)
    write_no_match_list(no_q, no_q_path, "No-match queries")

    # 4. No-match targets
    no_t_path = os.path.join(args.outdir,
                             f"no_target_matches_{q_name}_to_{t_name}.txt".replace(" ", "_"))
    no_t = build_no_match_list(list(t_contigs.keys()), t_seen)
    write_no_match_list(no_t, no_t_path, "No-match targets")

    print("\nDone!")
    if sampled:
        print(f"NOTE: Output was sampled to {args.max_lines} most significant matches.")
    print(f"Output files in: {args.outdir}")


if __name__ == "__main__":
    main()
