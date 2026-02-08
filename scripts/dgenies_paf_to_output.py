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

Usage:
    python3 dgenies_paf_to_output.py --paf map.paf --query-idx query.idx --target-idx target.idx -o output_dir/

    # With sorting (reorder query contigs to match target):
    python3 dgenies_paf_to_output.py --paf map.paf --query-idx query.idx --target-idx target.idx -o output_dir/ --sort

    # Limit lines read (for huge PAF files):
    python3 dgenies_paf_to_output.py --paf map.paf --query-idx query.idx --target-idx target.idx -o output_dir/ --max-lines 200000

    # With noise removal:
    python3 dgenies_paf_to_output.py --paf map.paf --query-idx query.idx --target-idx target.idx -o output_dir/ --remove-noise

Requires: Python 3.7+, matplotlib, numpy
"""

import argparse
import csv
import os
import re
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


# ---------------------------------------------------------------------------
# PAF sorting by significance (matches D-Genies Sorter)
# ---------------------------------------------------------------------------

def compute_significance(parts):
    """Compute D-Genies significance score for a PAF line."""
    q_start = int(parts[2])
    q_end = int(parts[3])
    t_start = int(parts[7])
    t_end = int(parts[8])
    matches = int(parts[9])
    block_len = int(parts[10])
    if block_len == 0:
        return 0.0
    euclidean = sqrt((q_end - q_start) ** 2 + (t_end - t_start) ** 2)
    identity = matches / block_len
    return euclidean * identity


def sort_paf_lines(lines):
    """Sort PAF lines by descending significance score."""
    scored = []
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        score = compute_significance(parts)
        scored.append((score, line))
    scored.sort(key=lambda x: x[0], reverse=True)
    return [s[1] for s in scored]


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
              max_lines=None, presorted=False):
    """Parse a PAF file and return match data.

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

    raw_lines = []
    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            raw_lines.append(line)

    # Sort if not presorted
    if not presorted:
        raw_lines = sort_paf_lines(raw_lines)

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
        q_name = q_contig_rename.get(q_name_orig, q_name_orig)
        t_name = t_contig_rename.get(t_name_orig, t_name_orig)

        # Skip if contig not in index (shouldn't happen, but be safe)
        if q_name not in q_abs_start or t_name not in t_abs_start:
            # Try original names
            if q_name_orig not in q_abs_start or t_name_orig not in t_abs_start:
                continue
            q_name = q_name_orig
            t_name = t_name_orig

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

        # Compute absolute coordinates
        # x = target, y = query (D-Genies convention)
        q_abs = q_abs_start[q_name]
        t_abs = t_abs_start[t_name]

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

def draw_dotplot(matches, q_order, q_contigs, q_abs_start, q_total_len,
                 t_order, t_contigs, t_abs_start, t_total_len,
                 q_name, t_name, out_path, dpi=300, figsize=None,
                 sampled=False):
    """Draw a D-Genies-style dot plot and save as PNG.

    Target (reference) on X-axis, Query on Y-axis (inverted like D-Genies).
    """
    if figsize is None:
        figsize = (12, 12)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # --- Draw alignment matches as line segments, grouped by identity class ---
    for idy_class in [0, 1, 2, 3]:
        class_matches = [m for m in matches if m["idy_class"] == idy_class]
        if not class_matches:
            continue

        segments = []
        for m in class_matches:
            # Normalize to 0..1 scale
            x1_n = m["x1"] / t_total_len
            x2_n = m["x2"] / t_total_len
            # Invert Y axis (D-Genies convention: origin at top-left)
            y1_n = 1.0 - (m["y1"] / q_total_len)
            y2_n = 1.0 - (m["y2"] / q_total_len)
            segments.append([(x1_n, y1_n), (x2_n, y2_n)])

        lc = LineCollection(segments, colors=IDY_COLORS[idy_class],
                            linewidths=0.4, alpha=0.85,
                            label=IDY_LABELS[idy_class])
        ax.add_collection(lc)

    # --- Draw contig boundaries ---
    # Target (x-axis) boundaries
    for cname in t_order:
        x = (t_abs_start[cname] + t_contigs[cname]) / t_total_len
        if x < 1.0:
            is_mix = cname.startswith("###MIX###")
            ax.axvline(x=x, color="#888888" if is_mix else "#cccccc",
                       linewidth=0.3, linestyle="-")

    # Query (y-axis) boundaries
    for cname in q_order:
        y = 1.0 - (q_abs_start[cname] + q_contigs[cname]) / q_total_len
        if y > 0.0:
            is_mix = cname.startswith("###MIX###")
            ax.axhline(y=y, color="#888888" if is_mix else "#cccccc",
                       linewidth=0.3, linestyle="-")

    # --- Axis labels (show contig names if not too many) ---
    max_ticks = 50

    # Target x-axis ticks
    if len(t_order) <= max_ticks:
        xtick_pos = []
        xtick_labels = []
        for cname in t_order:
            if cname.startswith("###MIX###"):
                continue
            mid = (t_abs_start[cname] + t_contigs[cname] / 2) / t_total_len
            xtick_pos.append(mid)
            xtick_labels.append(cname)
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(xtick_labels, rotation=90, fontsize=5, ha="center")
    else:
        ax.set_xticks([])

    # Query y-axis ticks
    if len(q_order) <= max_ticks:
        ytick_pos = []
        ytick_labels = []
        for cname in q_order:
            if cname.startswith("###MIX###"):
                continue
            mid = 1.0 - (q_abs_start[cname] + q_contigs[cname] / 2) / q_total_len
            ytick_pos.append(mid)
            ytick_labels.append(cname)
        ax.set_yticks(ytick_pos)
        ax.set_yticklabels(ytick_labels, fontsize=5)
    else:
        ax.set_yticks([])

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(f"Target: {t_name}", fontsize=10)
    ax.set_ylabel(f"Query: {q_name}", fontsize=10)
    ax.set_aspect("equal")

    title = f"Dot plot: {q_name} vs {t_name}"
    if sampled:
        title += " (sampled)"
    ax.set_title(title, fontsize=12)

    # --- Legend ---
    legend_handles = []
    for idy_class in [3, 2, 1, 0]:
        class_matches = [m for m in matches if m["idy_class"] == idy_class]
        if class_matches:
            legend_handles.append(
                mpatches.Patch(color=IDY_COLORS[idy_class],
                               label=IDY_LABELS[idy_class]))
    if legend_handles:
        ax.legend(handles=legend_handles, loc="upper right", fontsize=6,
                  framealpha=0.8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Dot plot written: {out_path}")


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
    parser.add_argument("--dpi", type=int, default=300,
                        help="DPI for PNG output (default: 300)")
    parser.add_argument("--figsize", type=float, nargs=2, default=None,
                        metavar=("WIDTH", "HEIGHT"),
                        help="Figure size in inches (default: 12 12)")
    parser.add_argument("--presorted", action="store_true",
                        help="PAF is already sorted by significance (skip re-sorting)")
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

    # --- Preserve original indexes (before merging) for sort/association ---
    q_order_orig = list(q_order)
    q_abs_start_orig = dict(q_abs_start)
    t_abs_start_orig = dict(t_abs_start)

    # --- Merge small contigs ---
    q_contig_rename = {}
    t_contig_rename = {}
    if not args.no_merge_small:
        q_order, q_contigs_m, q_abs_start, q_contig_rename = merge_small_contigs(
            q_order, q_contigs, q_abs_start, q_total)
        t_order, t_contigs_m, t_abs_start, t_contig_rename = merge_small_contigs(
            t_order, t_contigs, t_abs_start, t_total)
        n_q_merged = sum(1 for v in q_contig_rename.values() if v.startswith("###MIX###"))
        n_t_merged = sum(1 for v in t_contig_rename.values() if v.startswith("###MIX###"))
        if n_q_merged > 0:
            print(f"  Merged {n_q_merged} small query contigs into MIX groups")
        if n_t_merged > 0:
            print(f"  Merged {n_t_merged} small target contigs into MIX groups")
        # Use merged contigs dict for plotting but keep original for association table
        q_contigs_plot = q_contigs_m
        t_contigs_plot = t_contigs_m
    else:
        q_contigs_plot = q_contigs
        t_contigs_plot = t_contigs

    # --- Parse PAF ---
    print(f"Parsing PAF file...")
    matches, q_seen, t_seen, sampled = parse_paf(
        args.paf, q_abs_start, t_abs_start, q_contigs_plot, t_contigs_plot,
        q_contig_rename, t_contig_rename,
        max_lines=args.max_lines, presorted=args.presorted)
    print(f"  Parsed {len(matches)} alignment matches")
    if sampled:
        print(f"  (sampled to {args.max_lines} lines)")

    # --- Optional noise removal ---
    if args.remove_noise:
        matches = remove_noise(matches)
        print(f"  After noise removal: {len(matches)} matches")

    # --- Optional sorting ---
    if args.sort:
        print("Sorting query contigs to match target order...")
        q_order, q_contigs_plot, q_abs_start, q_reversed_new = sort_query_contigs(
            matches, q_order, q_contigs_plot, q_abs_start,
            t_order, t_contigs_plot, t_abs_start, t_total,
            t_abs_start_orig=t_abs_start_orig)
        # Update q_reversed with sorting results
        q_reversed.update(q_reversed_new)
        # Re-parse with new absolute starts
        print("Re-parsing PAF with sorted coordinates...")
        matches, q_seen, t_seen, sampled = parse_paf(
            args.paf, q_abs_start, t_abs_start, q_contigs_plot, t_contigs_plot,
            q_contig_rename, t_contig_rename,
            max_lines=args.max_lines, presorted=args.presorted)
        if args.remove_noise:
            matches = remove_noise(matches)
        print(f"  {len(matches)} matches after re-parse")

    # --- Generate outputs ---
    print("\nGenerating outputs...")

    # 1. Dot plot PNG
    if not args.no_plot:
        plot_path = os.path.join(args.outdir,
                                 f"map_{q_name}_to_{t_name}.png".replace(" ", "_"))
        draw_dotplot(matches, q_order, q_contigs_plot, q_abs_start, q_total,
                     t_order, t_contigs_plot, t_abs_start, t_total,
                     q_name, t_name, plot_path, dpi=args.dpi,
                     figsize=tuple(args.figsize) if args.figsize else None,
                     sampled=sampled)

    # 2. Association table (use original contigs/abs_start, not merged)
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
