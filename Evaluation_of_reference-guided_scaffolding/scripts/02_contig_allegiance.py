#!/usr/bin/env python3
"""
02_contig_allegiance.py

Contig allegiance analysis at breakpoints within compound scaffolds.
Classifies each breakpoint as a predicted fusion or predicted scaffolding
artifact based on mixing score and transition sharpness of pre-scaffolding
contig alignments to D. ebraccatus.

Input (in input/):
    ragtag.scaffold.renamed.agp
    dgenies_pairwise/PriEup-Assembly_DenEbr/map_*.paf
    GCF_027789765.1_*.fna.fai

Output (in output/):
    allegiance_breakpoint_results.tsv
    allegiance_painting.tsv
    allegiance_full_scaffolds.{png,pdf,tiff}
    allegiance_breakpoint_zooms.{png,pdf,tiff}

Dependencies: Python 3.8+, matplotlib

Author: Kopp K, Pristimantis euphronides genome project
"""


import sys
import csv
import os
from collections import defaultdict

# ============================================================
# Default breakpoints from tier (ii) (four-species synteny)
# ============================================================
DEFAULT_BREAKPOINTS = [
    ("scaffold_1", 120_000_000, "De6", "De5"),
    ("scaffold_3",  30_000_000, "De8", "De4"),
    ("scaffold_3",  95_000_000, "De4", "De15"),
    ("scaffold_4",  65_000_000, "De8", "De11"),
    ("scaffold_6",  70_000_000, "De12", "De13"),
]

COMPOUND_SCAFFOLDS = {"scaffold_1", "scaffold_3", "scaffold_4", "scaffold_6"}


# ============================================================
# Parsing
# ============================================================

def parse_agp(agp_path):
    contig_to_scaffold = {}
    scaffold_contigs = defaultdict(list)
    with open(agp_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[4] != 'W':
                continue
            scaffold = parts[0]
            scaf_start = int(parts[1])
            scaf_end = int(parts[2])
            contig = parts[5]
            orientation = parts[8]
            contig_to_scaffold[contig] = (scaffold, scaf_start, scaf_end, orientation)
            scaffold_contigs[scaffold].append((contig, scaf_start, scaf_end, orientation))
    for scaf in scaffold_contigs:
        scaffold_contigs[scaf].sort(key=lambda x: x[1])
    return contig_to_scaffold, dict(scaffold_contigs)


def parse_fai(fai_path, n_chr=15):
    sizes = []
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sizes.append((parts[0], int(parts[1])))
    sizes.sort(key=lambda x: x[1], reverse=True)
    chr_list = sizes[:n_chr]
    chr_map = {name: f"De{i}" for i, (name, size) in enumerate(chr_list, 1)}
    return chr_map, chr_list


def parse_paf(paf_path, min_align_len=5000, min_mapq=10):
    contig_alignments = defaultdict(list)
    total = 0
    passed = 0
    with open(paf_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            total += 1
            query_name = parts[0]
            query_len = int(parts[1])
            query_start = int(parts[2])
            query_end = int(parts[3])
            strand = parts[4]
            target_name = parts[5]
            target_len = int(parts[6])
            target_start = int(parts[7])
            target_end = int(parts[8])
            n_matches = int(parts[9])
            align_block_len = int(parts[10])
            mapq = int(parts[11])
            align_len = query_end - query_start
            if align_len < min_align_len or mapq < min_mapq:
                continue
            passed += 1
            contig_alignments[query_name].append({
                'target': target_name, 'target_len': target_len,
                'target_start': target_start, 'target_end': target_end,
                'query_start': query_start, 'query_end': query_end,
                'query_len': query_len, 'align_len': align_len,
                'mapq': mapq, 'strand': strand,
                'n_matches': n_matches, 'block_len': align_block_len,
            })
    return dict(contig_alignments), total, passed


# ============================================================
# Allegiance determination (unchanged)
# ============================================================

def determine_contig_allegiance(alignments, chr_map):
    chr_bases = defaultdict(int)
    for a in alignments:
        de_label = chr_map.get(a['target'])
        if de_label:
            chr_bases[de_label] += a['align_len']
    if not chr_bases:
        return None, {}, False, set()
    total = sum(chr_bases.values())
    primary = max(chr_bases, key=chr_bases.get)
    significant = {ch for ch, bases in chr_bases.items() if bases >= total * 0.15}
    is_split = len(significant) >= 2
    return primary, dict(chr_bases), is_split, significant


def paint_scaffold(scaffold_contig_list, contig_alignments, chr_map):
    painted = []
    for contig, scaf_start, scaf_end, orientation in scaffold_contig_list:
        aligns = contig_alignments.get(contig, [])
        primary, allegiances, is_split, split_chrs = \
            determine_contig_allegiance(aligns, chr_map)
        painted.append({
            'contig': contig,
            'scaf_start': scaf_start, 'scaf_end': scaf_end,
            'scaf_midpoint': (scaf_start + scaf_end) // 2,
            'contig_len': scaf_end - scaf_start + 1,
            'orientation': orientation,
            'primary_chr': primary, 'allegiances': allegiances,
            'is_split': is_split, 'split_chrs': split_chrs,
            'n_alignments': len(aligns),
            'total_aligned_bases': sum(a['align_len'] for a in aligns),
        })
    return painted


# ============================================================
# Transition-sharpness analysis
# ============================================================

def find_transition_point(painted, expected_left, expected_right):
    """
    Walk through mapped contigs and find where the dominant allegiance
    switches from expected_left to expected_right.

    Returns (transition_pos, transition_width_Mb) or (None, None).
    Transition width = distance from last expected_left-dominated window
    to first expected_right-dominated window.
    """
    mapped = [p for p in painted if p['primary_chr'] is not None]
    if len(mapped) < 10:
        return None, None

    win = 10
    # Record where each window's dominant allegiance is
    window_results = []
    for i in range(len(mapped) - win + 1):
        window = mapped[i:i + win]
        counts = defaultdict(int)
        for c in window:
            counts[c['primary_chr']] += 1
        dominant = max(counts, key=counts.get)
        dominant_frac = counts[dominant] / win
        mid_pos = window[win // 2]['scaf_midpoint']
        left_count = counts.get(expected_left, 0)
        right_count = counts.get(expected_right, 0)
        window_results.append({
            'mid_pos': mid_pos,
            'dominant': dominant,
            'dominant_frac': dominant_frac,
            'left_count': left_count,
            'right_count': right_count,
        })

    # Find the switch: last window dominated by expected_left,
    # first window dominated by expected_right
    last_left_pos = None
    first_right_pos = None
    transition_pos = None

    for wr in window_results:
        if wr['dominant'] == expected_left:
            last_left_pos = wr['mid_pos']

    for wr in window_results:
        if wr['dominant'] == expected_right and last_left_pos is not None:
            if wr['mid_pos'] > last_left_pos:
                first_right_pos = wr['mid_pos']
                break

    if last_left_pos is not None and first_right_pos is not None:
        transition_pos = (last_left_pos + first_right_pos) // 2
        transition_width = (first_right_pos - last_left_pos) / 1e6
    else:
        transition_pos = None
        transition_width = None

    return transition_pos, transition_width


def analyze_breakpoint_zone(painted, breakpoint_pos, expected_left,
                                expected_right, zone_size=5_000_000):
    """
    Zone analysis focused on transition sharpness.
    """
    zone_left = breakpoint_pos - zone_size
    zone_right = breakpoint_pos + zone_size

    # All mapped contigs in zone
    zone_contigs = [p for p in painted
                    if p['scaf_end'] >= zone_left and p['scaf_start'] <= zone_right
                    and p['primary_chr'] is not None]

    # Split contigs (still reported, just not used for classification)
    split_in_zone = [p for p in zone_contigs if p['is_split']]
    relevant_splits = [p for p in split_in_zone
                       if expected_left in p['split_chrs']
                       and expected_right in p['split_chrs']]

    # ----- METRIC 1: Zone mixing score -----
    # How well-represented are BOTH expected chromosomes in the zone?
    n_left = sum(1 for c in zone_contigs if c['primary_chr'] == expected_left)
    n_right = sum(1 for c in zone_contigs if c['primary_chr'] == expected_right)
    n_expected_total = n_left + n_right
    n_zone_total = len(zone_contigs)

    if n_expected_total > 0:
        minority_count = min(n_left, n_right)
        mixing_score = minority_count / n_expected_total  # 0 = pure, 0.5 = perfect mix
    else:
        mixing_score = 0.0

    # ----- METRIC 2: Both-present in zone? -----
    both_present = n_left > 0 and n_right > 0

    # ----- METRIC 3: Find actual transition point and width -----
    transition_pos, transition_width = find_transition_point(
        painted, expected_left, expected_right)

    offset_Mb = abs(transition_pos - breakpoint_pos) / 1e6 if transition_pos else None

    # ----- METRIC 4: Wider zone scan (±15 Mb) for allegiance profile -----
    wide_left = breakpoint_pos - 15_000_000
    wide_right = breakpoint_pos + 15_000_000
    wide_contigs = [p for p in painted
                    if p['scaf_end'] >= wide_left and p['scaf_start'] <= wide_right
                    and p['primary_chr'] is not None]
    wide_n_left = sum(1 for c in wide_contigs if c['primary_chr'] == expected_left)
    wide_n_right = sum(1 for c in wide_contigs if c['primary_chr'] == expected_right)
    wide_total = wide_n_left + wide_n_right
    wide_mixing = min(wide_n_left, wide_n_right) / wide_total if wide_total > 0 else 0.0

    # Full allegiance counts
    zone_counts = defaultdict(int)
    for c in zone_contigs:
        zone_counts[c['primary_chr']] += 1

    wide_counts = defaultdict(int)
    for c in wide_contigs:
        wide_counts[c['primary_chr']] += 1

    return {
        'n_contigs_in_zone': n_zone_total,
        'n_left_in_zone': n_left,
        'n_right_in_zone': n_right,
        'mixing_score_zone': mixing_score,
        'both_present_zone': both_present,
        'n_contigs_wide': len(wide_contigs),
        'n_left_wide': wide_n_left,
        'n_right_wide': wide_n_right,
        'mixing_score_wide': wide_mixing,
        'both_present_wide': wide_n_left > 0 and wide_n_right > 0,
        'transition_pos': transition_pos,
        'transition_width_Mb': transition_width,
        'offset_Mb': offset_Mb,
        'zone_allegiance_counts': dict(zone_counts),
        'wide_allegiance_counts': dict(wide_counts),
        'n_split_in_zone': len(split_in_zone),
        'n_relevant_splits': len(relevant_splits),
        'relevant_split_names': [p['contig'] for p in relevant_splits],
    }


def classify_breakpoint(zr, expected_left, expected_right):
    """
    Classifier using transition sharpness.

    Key principle: a predicted fusion has a GRADUAL, DIFFUSE transition where
    both D. ebraccatus chromosomes intermingle. A scaffolding artifact has
    a SHARP boundary.

    Scoring:
      +2  both chromosomes well-represented in ±5 Mb zone (mixing >= 0.2)
      +1  both chromosomes present in ±5 Mb zone (mixing 0.05-0.2)
      +1  transition width > 5 Mb (gradual/diffuse)
      +1  mixing score in wide zone >= 0.3
      -2  only one chromosome in ±5 Mb zone
      -1  transition width < 2 Mb (sharp)
      -1  zone dominated by unexpected chromosome (transition far from expected)
    """
    reasons = []
    score = 0

    # --- Zone mixing (±5 Mb) ---
    ms = zr['mixing_score_zone']
    n_left = zr['n_left_in_zone']
    n_right = zr['n_right_in_zone']

    if ms >= 0.2:
        score += 2
        reasons.append(
            f"Good mixing in +-5Mb zone: {expected_left}={n_left}, "
            f"{expected_right}={n_right} (score={ms:.2f})"
        )
    elif zr['both_present_zone']:
        score += 1
        reasons.append(
            f"Both present in +-5Mb zone: {expected_left}={n_left}, "
            f"{expected_right}={n_right} (score={ms:.2f})"
        )
    else:
        score -= 2
        dominant = expected_left if n_left > 0 else (expected_right if n_right > 0 else "neither")
        reasons.append(
            f"Only {dominant} in +-5Mb zone: {expected_left}={n_left}, "
            f"{expected_right}={n_right}"
        )

    # --- Check if zone is dominated by a THIRD chromosome (breakpoint position off) ---
    n_other = zr['n_contigs_in_zone'] - n_left - n_right
    if n_other > n_left + n_right and zr['n_contigs_in_zone'] > 5:
        score -= 1
        reasons.append(
            f"Zone dominated by other chromosomes ({n_other} other vs "
            f"{n_left + n_right} expected)"
        )

    # --- Wide zone mixing (±15 Mb) for context ---
    wms = zr['mixing_score_wide']
    if wms >= 0.3:
        score += 1
        reasons.append(
            f"Good mixing in +-15Mb: {expected_left}={zr['n_left_wide']}, "
            f"{expected_right}={zr['n_right_wide']} (score={wms:.2f})"
        )
    elif not zr['both_present_wide']:
        score -= 1
        reasons.append(
            f"Only one chr in +-15Mb: {expected_left}={zr['n_left_wide']}, "
            f"{expected_right}={zr['n_right_wide']}"
        )

    # --- Transition width ---
    tw = zr['transition_width_Mb']
    if tw is not None:
        if tw > 5:
            score += 1
            reasons.append(f"Diffuse transition: {tw:.1f} Mb wide")
        elif tw < 2:
            score -= 1
            reasons.append(f"Sharp transition: {tw:.1f} Mb wide")
        else:
            reasons.append(f"Moderate transition: {tw:.1f} Mb wide")
    else:
        reasons.append("No clear transition detected between expected pair")

    # --- Transition offset ---
    if zr['offset_Mb'] is not None:
        reasons.append(f"Transition at {zr['transition_pos']/1e6:.1f} Mb "
                       f"(offset {zr['offset_Mb']:.1f} Mb from expected)")

    # --- Split contigs (bonus, not primary) ---
    if zr['n_relevant_splits'] > 0:
        score += 1
        reasons.append(f"+{zr['n_relevant_splits']} split-allegiance contigs (bonus)")

    # --- Classify ---
    if score >= 1:
        cls = "PREDICTED_FUSION"
    elif score <= -2:
        cls = "PREDICTED_SCAFFOLDING_ARTIFACT"
    else:
        cls = "UNCLASSIFIED"

    return cls, score, "; ".join(reasons)


# ============================================================
# Plotting
# ============================================================

def make_plot(painted_scaffolds, breakpoints, zone_results, classifications,
              scores, chr_map, outdir):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("  [WARN] matplotlib not available, skipping plots.", file=sys.stderr)
        return

    all_de = sorted(set(chr_map.values()), key=lambda x: int(x[2:]))
    cmap = matplotlib.colormaps.get_cmap('tab20').resampled(max(len(all_de), 1))
    de_colors = {label: cmap(i) for i, label in enumerate(all_de)}
    de_colors[None] = (0.85, 0.85, 0.85, 1.0)

    cls_colors = {
        'PREDICTED_FUSION': '#2ca02c',
        'UNCLASSIFIED': '#ffbb78',
        'PREDICTED_SCAFFOLDING_ARTIFACT': '#d62728',
    }

    # ---- Full scaffold paintings ----
    unique_scaffolds = list(dict.fromkeys(bp[0] for bp in breakpoints))

    fig, axes = plt.subplots(len(unique_scaffolds), 1,
                              figsize=(6.69, 3.0 * len(unique_scaffolds)),
                              squeeze=False)

    for row, scaffold in enumerate(unique_scaffolds):
        ax = axes[row, 0]
        painted = painted_scaffolds.get(scaffold, [])
        if not painted:
            ax.text(0.5, 0.5, f"No data for {scaffold}", transform=ax.transAxes)
            continue
        scaffold_size = max(p['scaf_end'] for p in painted)

        for p in painted:
            color = de_colors.get(p['primary_chr'], de_colors[None])
            w = (p['scaf_end'] - p['scaf_start']) / 1e6
            ax.barh(0, w, left=p['scaf_start'] / 1e6,
                    height=0.6, color=color, edgecolor='none', alpha=0.8)
            if p['is_split']:
                ax.barh(0, w, left=p['scaf_start'] / 1e6,
                        height=0.6, color='none', edgecolor='black',
                        linewidth=1.2)

        for bp in breakpoints:
            if bp[0] != scaffold:
                continue
            bp_key = f"{bp[0]}:{bp[1]}"
            cls = classifications.get(bp_key, "UNCLASSIFIED")
            sc = scores.get(bp_key, 0)
            color = cls_colors.get(cls, 'gray')
            ax.axvline(x=bp[1] / 1e6, color=color, linewidth=2.5,
                        linestyle='--', zorder=10)
            ax.text(bp[1] / 1e6, 0.75,
                    f"  {bp[2]}|{bp[3]}\n  {cls.replace('_', ' ')} ({sc:+d})",
                    fontsize=8, fontweight='bold', color=color,
                    ha='center', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                              edgecolor=color, alpha=0.85))

        ax.set_xlim(0, scaffold_size / 1e6)
        ax.set_ylim(-0.5, 1.4)
        ax.set_xlabel('Position (Mb)', fontsize=11)
        ax.set_yticks([0])
        ax.set_yticklabels(['Contig allegiance'], fontsize=9)
        ax.set_title(
            f'{scaffold} ({scaffold_size/1e6:.0f} Mb)',
            fontsize=11, loc='left')
        ax.grid(axis='x', alpha=0.3)

    relevant_de = set()
    for bp in breakpoints:
        relevant_de.update([bp[2], bp[3]])
    for scaf in unique_scaffolds:
        for p in painted_scaffolds.get(scaf, []):
            if p['primary_chr']:
                relevant_de.add(p['primary_chr'])

    patches = [mpatches.Patch(color=de_colors.get(de, 'gray'), label=de)
               for de in sorted(relevant_de, key=lambda x: int(x[2:]))
               if de in de_colors]
    patches.append(mpatches.Patch(color=de_colors[None], label='No alignment'))
    patches.append(mpatches.Patch(facecolor='none', edgecolor='black',
                                   linewidth=1.2, label='Split-allegiance'))
    fig.legend(handles=patches, loc='lower center',
               ncol=min(10, len(patches)), fontsize=9,
               bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    plt.savefig(os.path.join(outdir, 'allegiance_full_scaffolds.png'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'allegiance_full_scaffolds.pdf'),
                bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'allegiance_full_scaffolds.tiff'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Full scaffold plot saved (.png, .pdf, .tiff)", file=sys.stderr)

    # ---- Zoomed breakpoint plots ----
    fig2, axes2 = plt.subplots(len(breakpoints), 1,
                                figsize=(6.69, 3.5 * len(breakpoints)),
                                squeeze=False)

    for idx, bp in enumerate(breakpoints):
        ax = axes2[idx, 0]
        scaffold, bp_pos, left_anc, right_anc = bp
        painted = painted_scaffolds.get(scaffold, [])
        scaffold_size = max(p['scaf_end'] for p in painted) if painted else 0
        zr = zone_results[idx]

        vr_radius = 20_000_000
        vl = max(0, bp_pos - vr_radius)
        vr_end = min(scaffold_size, bp_pos + vr_radius)

        visible = [p for p in painted
                   if p['scaf_end'] >= vl and p['scaf_start'] <= vr_end]

        for p in visible:
            color = de_colors.get(p['primary_chr'], de_colors[None])
            x0 = max(p['scaf_start'], vl) / 1e6
            x1 = min(p['scaf_end'], vr_end) / 1e6
            ax.barh(0, x1 - x0, left=x0, height=0.6, color=color,
                    edgecolor='none', alpha=0.8)
            if p['is_split']:
                ax.barh(0, x1 - x0, left=x0, height=0.6, color='none',
                        edgecolor='black', linewidth=1.5)

        bp_key = f"{scaffold}:{bp_pos}"
        cls = classifications.get(bp_key, "UNCLASSIFIED")
        sc = scores.get(bp_key, 0)
        color = cls_colors.get(cls, 'gray')

        # Expected breakpoint
        ax.axvline(x=bp_pos / 1e6, color=color, linewidth=3,
                    linestyle='-', zorder=10, label='Expected BP')

        # Detected transition
        if zr.get('transition_pos'):
            ax.axvline(x=zr['transition_pos'] / 1e6, color='purple',
                        linewidth=2, linestyle=':', zorder=9,
                        label='Detected transition')

        # Zone shading (±5 Mb)
        z5l = max(vl, bp_pos - 5_000_000) / 1e6
        z5r = min(vr_end, bp_pos + 5_000_000) / 1e6
        ax.axvspan(z5l, z5r, alpha=0.08, color=color)

        # Ancestral labels
        ax.text((vl / 1e6 + bp_pos / 1e6) / 2, 0.55, left_anc,
                fontsize=16, fontweight='bold', color='blue', alpha=0.4,
                ha='center')
        ax.text((bp_pos / 1e6 + vr_end / 1e6) / 2, 0.55, right_anc,
                fontsize=16, fontweight='bold', color='red', alpha=0.4,
                ha='center')

        # Info box
        ms5 = zr.get('mixing_score_zone', 0)
        ms15 = zr.get('mixing_score_wide', 0)
        tw = zr.get('transition_width_Mb')
        tw_str = f"{tw:.1f} Mb" if tw is not None else "N/A"
        info = (f"Classification: {cls.replace('_', ' ')} (score={sc:+d})\n"
                f"Zone mixing (+-5Mb): {ms5:.2f}\n"
                f"Wide mixing (+-15Mb): {ms15:.2f}\n"
                f"Transition width: {tw_str}\n"
                f"Zone: {left_anc}={zr.get('n_left_in_zone', 0)}, "
                f"{right_anc}={zr.get('n_right_in_zone', 0)}")
        ax.text(0.02, 0.95, info, transform=ax.transAxes, fontsize=9,
                va='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        ax.set_xlim(vl / 1e6, vr_end / 1e6)
        ax.set_ylim(-0.5, 1.3)
        ax.set_xlabel('Position (Mb)')
        ax.set_yticks([])
        ax.legend(fontsize=8, loc='upper right')
        ax.set_title(
            f'{scaffold} @ {bp_pos/1e6:.0f} Mb - {left_anc} | {right_anc}',
            fontsize=11, loc='left')

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'allegiance_breakpoint_zooms.png'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'allegiance_breakpoint_zooms.pdf'),
                bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'allegiance_breakpoint_zooms.tiff'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Zoomed breakpoint plot saved (.png, .pdf, .tiff)", file=sys.stderr)


# ============================================================
# Main
# ============================================================

# ── Hardcoded input paths ──────────────────────────────────────────
BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

AGP = os.path.join(BASEDIR, 'input', 'ragtag.scaffold.renamed.agp')
PAF = os.path.join(BASEDIR,
    'input', 'dgenies_pairwise', 'PriEup-Assembly_DenEbr', 'map_final_medaka_polished_assembly_consensus_to_GCF_027789765.1_aDenEbr1.pat_genomic.paf')
DENEBR_FAI = os.path.join(BASEDIR,
    'input', 'GCF_027789765.1_aDenEbr1.pat_genomic.fna.fai')
OUTDIR = os.path.join(BASEDIR, 'output')
MIN_ALIGN = 5000
MIN_MAPQ = 10
ZONE_SIZE = 5_000_000


def main():

    os.makedirs(OUTDIR, exist_ok=True)

    breakpoints = DEFAULT_BREAKPOINTS

    print("=" * 70, file=sys.stderr)
    print("CONTIG ALLEGIANCE ANALYSIS", file=sys.stderr)
    print("=" * 70, file=sys.stderr)

    # ---- Load inputs ----
    print(f"\n  Loading AGP: {AGP}", file=sys.stderr)
    contig_to_scaffold, scaffold_contigs = parse_agp(AGP)
    print(f"    Contigs: {len(contig_to_scaffold)}  |  "
          f"Scaffolds: {len(scaffold_contigs)}", file=sys.stderr)

    print(f"\n  Loading DenEbr .fai: {DENEBR_FAI}", file=sys.stderr)
    chr_map, chr_sizes = parse_fai(DENEBR_FAI, n_chr=15)
    print(f"    D. ebraccatus chromosomes:", file=sys.stderr)
    for name, size in chr_sizes:
        print(f"      {chr_map[name]:>5s} = {name:20s}  ({size/1e6:.1f} Mb)",
              file=sys.stderr)

    print(f"\n  Loading PAF: {PAF}", file=sys.stderr)
    print(f"    Filters: align_len >= {MIN_ALIGN}, mapq >= {MIN_MAPQ}",
          file=sys.stderr)
    contig_alignments, total_aln, passed_aln = parse_paf(
        PAF, min_align_len=MIN_ALIGN, min_mapq=MIN_MAPQ)
    print(f"    Total alignments: {total_aln:,}", file=sys.stderr)
    print(f"    Passed filters:   {passed_aln:,}", file=sys.stderr)
    print(f"    Contigs with alignments: {len(contig_alignments):,}",
          file=sys.stderr)

    # Name match
    agp_names = set(contig_to_scaffold.keys())
    paf_names = set(contig_alignments.keys())
    overlap = agp_names & paf_names
    print(f"\n  Name matching: {len(overlap):,} / {len(agp_names):,} AGP contigs "
          f"have PAF alignments ({len(overlap)/len(agp_names)*100:.0f}%)",
          file=sys.stderr)
    if len(overlap) == 0:
        print("  [FATAL] No names match!", file=sys.stderr)
        sys.exit(1)

    # ---- Paint compound scaffolds ----
    print(f"\n  Painting compound scaffolds...", file=sys.stderr)
    scaffolds_to_paint = COMPOUND_SCAFFOLDS | set(bp[0] for bp in breakpoints)
    painted_scaffolds = {}

    for scaffold in sorted(scaffolds_to_paint):
        if scaffold not in scaffold_contigs:
            print(f"    [WARN] {scaffold} not in AGP", file=sys.stderr)
            continue
        painted = paint_scaffold(
            scaffold_contigs[scaffold], contig_alignments, chr_map)
        painted_scaffolds[scaffold] = painted

        n_mapped = sum(1 for p in painted if p['primary_chr'] is not None)
        n_split = sum(1 for p in painted if p['is_split'])
        n_total = len(painted)

        allg_bases = defaultdict(int)
        for p in painted:
            if p['primary_chr']:
                allg_bases[p['primary_chr']] += p['contig_len']
        total_mb = sum(allg_bases.values())

        print(f"\n    {scaffold}: {n_total} contigs, "
              f"{n_mapped} mapped ({n_mapped/n_total*100:.0f}%), "
              f"{n_split} split", file=sys.stderr)
        for de in sorted(allg_bases, key=lambda x: allg_bases[x], reverse=True):
            pct = allg_bases[de] / total_mb * 100 if total_mb else 0
            print(f"      {de:>5s}: {allg_bases[de]/1e6:>7.1f} Mb ({pct:.1f}%)",
                  file=sys.stderr)

    # ---- Analyze breakpoints ----
    print(f"\n{'='*70}", file=sys.stderr)
    print("BREAKPOINT ANALYSIS", file=sys.stderr)
    print(f"{'='*70}", file=sys.stderr)

    zone_results = []
    classifications = {}
    scores_dict = {}

    for bp in breakpoints:
        scaffold, bp_pos, left_anc, right_anc = bp
        painted = painted_scaffolds.get(scaffold, [])

        if not painted:
            print(f"\n  [SKIP] {scaffold}", file=sys.stderr)
            zone_results.append({})
            continue

        zr = analyze_breakpoint_zone(
            painted, bp_pos, left_anc, right_anc, ZONE_SIZE)
        cls, score, reasons = classify_breakpoint(zr, left_anc, right_anc)

        bp_key = f"{scaffold}:{bp_pos}"
        classifications[bp_key] = cls
        scores_dict[bp_key] = score
        zone_results.append(zr)

        print(f"\n  {'='*60}", file=sys.stderr)
        print(f"  {scaffold} @ {bp_pos/1e6:.0f} Mb - {left_anc} | {right_anc}",
              file=sys.stderr)
        print(f"  {'='*60}", file=sys.stderr)
        print(f"    Zone (+-5 Mb): {zr['n_contigs_in_zone']} contigs | "
              f"{left_anc}={zr['n_left_in_zone']}, "
              f"{right_anc}={zr['n_right_in_zone']} | "
              f"mixing={zr['mixing_score_zone']:.3f}",
              file=sys.stderr)
        print(f"    Wide (+-15 Mb): {zr['n_contigs_wide']} contigs | "
              f"{left_anc}={zr['n_left_wide']}, "
              f"{right_anc}={zr['n_right_wide']} | "
              f"mixing={zr['mixing_score_wide']:.3f}",
              file=sys.stderr)
        if zr['transition_pos']:
            print(f"    Transition: {zr['transition_pos']/1e6:.1f} Mb "
                  f"(width={zr['transition_width_Mb']:.1f} Mb, "
                  f"offset={zr['offset_Mb']:.1f} Mb)",
                  file=sys.stderr)
        else:
            print(f"    Transition: not detected for {left_anc}->{right_anc}",
                  file=sys.stderr)
        print(f"    Zone allegiance: {dict(zr['zone_allegiance_counts'])}",
              file=sys.stderr)
        print(f"    Split contigs ({left_anc}/{right_anc}): "
              f"{zr['n_relevant_splits']}", file=sys.stderr)
        print(f"    -> CLASSIFICATION: {cls} (score={score:+d})", file=sys.stderr)
        print(f"       {reasons}", file=sys.stderr)

    # ---- Summary ----
    print(f"\n{'='*70}", file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print(f"{'='*70}", file=sys.stderr)

    markers = {
        'PREDICTED_FUSION': ' FUSION ',
        'UNCLASSIFIED': '  ???   ',
        'PREDICTED_SCAFFOLDING_ARTIFACT': 'ARTIFACT',
    }

    n_real = n_artifact = n_ambiguous = 0

    for bp, zr in zip(breakpoints, zone_results):
        if not zr:
            continue
        scaffold, bp_pos, left_anc, right_anc = bp
        bp_key = f"{scaffold}:{bp_pos}"
        cls = classifications.get(bp_key, "?")
        sc = scores_dict.get(bp_key, 0)
        marker = markers.get(cls, '???')
        ms = zr.get('mixing_score_zone', 0)
        tw = zr.get('transition_width_Mb')
        tw_str = f"{tw:.1f}Mb" if tw is not None else "N/A"

        print(f"  [{marker}]  {scaffold:12s} @ {bp_pos/1e6:>5.0f} Mb  "
              f"({left_anc:>4s}|{right_anc:<5s})  "
              f"mix={ms:.2f}  width={tw_str}  score={sc:+d}",
              file=sys.stderr)

        if 'REAL' in cls:
            n_real += 1
        elif 'ARTIFACT' in cls:
            n_artifact += 1
        else:
            n_ambiguous += 1

    print(f"\n  Real: {n_real}  |  Artifact: {n_artifact}  |  "
          f"Ambiguous: {n_ambiguous}", file=sys.stderr)
    print(f"  Expected: 2 real + 3 artifacts (n=16 - n=13 = 3 extra chr)",
          file=sys.stderr)

    print(f"\n  KARYOTYPE IMPLICATIONS:", file=sys.stderr)
    for bp in breakpoints:
        scaffold, bp_pos, left_anc, right_anc = bp
        bp_key = f"{scaffold}:{bp_pos}"
        cls = classifications.get(bp_key, "?")
        if 'ARTIFACT' in cls:
            print(f"    {scaffold} splits at {bp_pos/1e6:.0f} Mb -> "
                  f"{left_anc} and {right_anc} are SEPARATE in P. euphronides",
                  file=sys.stderr)
        elif 'REAL' in cls:
            print(f"    {scaffold} INTACT -> "
                  f"{left_anc}+{right_anc} fused in both species",
                  file=sys.stderr)
        else:
            print(f"    {scaffold} @ {bp_pos/1e6:.0f} Mb: inconclusive",
                  file=sys.stderr)

    # ---- Write results TSV ----
    results_path = os.path.join(OUTDIR, 'allegiance_breakpoint_results.tsv')
    with open(results_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'scaffold', 'breakpoint_pos_bp', 'left_ancestor', 'right_ancestor',
            'classification', 'score',
            'mixing_zone_5Mb', 'mixing_wide_15Mb',
            'n_left_zone', 'n_right_zone', 'both_present_zone',
            'n_left_wide', 'n_right_wide',
            'transition_pos_bp', 'transition_width_Mb', 'offset_Mb',
            'n_split_relevant', 'zone_allegiance', 'reasons',
        ])
        for bp, zr in zip(breakpoints, zone_results):
            if not zr:
                continue
            scaffold, bp_pos, left_anc, right_anc = bp
            bp_key = f"{scaffold}:{bp_pos}"
            cls = classifications.get(bp_key, "?")
            sc = scores_dict.get(bp_key, 0)
            _, _, reasons = classify_breakpoint(zr, left_anc, right_anc)
            writer.writerow([
                scaffold, bp_pos, left_anc, right_anc, cls, sc,
                f"{zr['mixing_score_zone']:.3f}",
                f"{zr['mixing_score_wide']:.3f}",
                zr['n_left_in_zone'], zr['n_right_in_zone'],
                zr['both_present_zone'],
                zr['n_left_wide'], zr['n_right_wide'],
                zr.get('transition_pos', ''),
                f"{zr['transition_width_Mb']:.1f}" if zr.get('transition_width_Mb') is not None else '',
                f"{zr['offset_Mb']:.1f}" if zr.get('offset_Mb') is not None else '',
                zr['n_relevant_splits'],
                str(zr['zone_allegiance_counts']),
                reasons,
            ])
    print(f"\n  [OK] Results: {results_path}", file=sys.stderr)

    # ---- Write painting TSV ----
    painting_path = os.path.join(OUTDIR, 'allegiance_painting.tsv')
    with open(painting_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'scaffold', 'scaf_start', 'scaf_end', 'contig', 'contig_len',
            'primary_chr', 'is_split', 'split_chrs',
            'n_alignments', 'total_aligned_bases', 'allegiances'])
        for scaffold in sorted(painted_scaffolds):
            for p in painted_scaffolds[scaffold]:
                writer.writerow([
                    scaffold, p['scaf_start'], p['scaf_end'], p['contig'],
                    p['contig_len'], p['primary_chr'] or 'none',
                    p['is_split'],
                    ','.join(sorted(p['split_chrs'])) if p['split_chrs'] else '',
                    p['n_alignments'], p['total_aligned_bases'],
                    str(p['allegiances']) if p['allegiances'] else '',
                ])
    print(f"  [OK] Painting: {painting_path}", file=sys.stderr)

    # ---- Plots ----
    make_plot(painted_scaffolds, breakpoints, zone_results, classifications,
              scores_dict, chr_map, OUTDIR)

    print(f"\nDone.", file=sys.stderr)


if __name__ == '__main__':
    main()
