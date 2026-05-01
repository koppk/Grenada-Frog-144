#!/usr/bin/env python3
"""
03_karyotype_schematic.py

Generates the karyotype reconstruction figure from tier (ii) and
tier (iii) analysis outputs plus Z-candidate region boundaries.

Input:
    input/Pristimantis_euphronides.genome.fasta.fai
    input/z_candidate_regions.tsv
    output/Table_compound_chromosomes.tsv
    output/Table_karyotype_reconstruction.tsv
    output/allegiance_breakpoint_results.tsv

Output (in output/):
    Fig_karyotype_schematic.{png,pdf,tiff}

Dependencies: Python 3.8+, matplotlib

Author: Kopp K, Pristimantis euphronides genome project
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


# ── Aesthetics ─────────────────────────────────────────────────────
DE_COLOURS = {
    'De1':  '#7986CB', 'De2':  '#5C6BC0', 'De3':  '#42A5F5',
    'De4':  '#66BB6A', 'De5':  '#4CAF50', 'De6':  '#E07B54',
    'De7':  '#26A69A', 'De8':  '#FF7043', 'De9':  '#EF5350',
    'De10': '#EC407A', 'De11': '#8D6E63', 'De12': '#BDBDBD',
    'De13': '#FDD835', 'De14': '#78909C', 'De15': '#AB47BC',
}
FL = 7.5; FU = 7; FB = 7; FZ = 7; FO = 7; FG = 7


# ── Data loading ───────────────────────────────────────────────────

def load_fai(path):
    """Read .fai -> {scaffold_num: size_bp} for scaffold_1..scaffold_13."""
    sizes = {}
    with open(path) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            name = cols[0]
            if name.startswith('scaffold_'):
                try:
                    num = int(name.split('_')[1])
                    if 1 <= num <= 13:
                        sizes[num] = int(cols[1])
                except (ValueError, IndexError):
                    pass
    if len(sizes) != 13:
        print('WARNING: expected 13 scaffolds in .fai, found {}'.format(
            len(sizes)), file=sys.stderr)
    return sizes


def load_compound_and_karyotype(compound_path, karyotype_path):
    """Build ancestry map: {scaffold_num: [(de_unit, start_mb, end_mb), ...]}
    
    Compound scaffolds: read spatial order from Table_compound_chromosomes.tsv
    and infer segment boundaries from breakpoints.
    Simple scaffolds: read from Table_karyotype_reconstruction.tsv.
    """
    # First pass: identify compound scaffolds and their units in order
    compound = {}  # {scaffold_num: [de_unit, ...]} in row order
    compound_aligned = {}  # {scaffold_num: {de_unit: aligned_mb}}
    with open(compound_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            scaf = row['PriEup_scaffold'].strip()
            num = int(scaf.replace('scaffold_', ''))
            de = row['DenEbr_unit'].strip()
            aligned = float(row['aligned_Mb'].strip())
            compound.setdefault(num, []).append(de)
            compound_aligned.setdefault(num, {})[de] = aligned
    compound_scaffolds = set(compound.keys())

    # Simple scaffolds from karyotype reconstruction
    simple = {}  # {scaffold_num: de_unit}
    with open(karyotype_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            de = row['DenEbr_unit'].strip()
            scaffolds_str = row['PriEup_scaffolds'].strip()
            # Parse entries like "scaffold_2 (42 Mb)"
            for part in scaffolds_str.split(';'):
                part = part.strip()
                if not part:
                    continue
                scaf_name = part.split('(')[0].strip()
                num = int(scaf_name.replace('scaffold_', ''))
                if num not in compound_scaffolds:
                    simple[num] = de

    return compound, compound_aligned, simple


def load_allegiance(path):
    """Read allegiance_breakpoint_results.tsv ->
    breakpoints: {scaffold_num: [(pos_mb, classification), ...]}
    spatial:     {scaffold_num: [(pos_mb, left_unit, right_unit), ...]}
    
    Classification mapping (harmonised):
        PREDICTED_SCAFFOLDING_ARTIFACT -> ART  (Predicted scaffolding artifact)
        PREDICTED_FUSION -> FUSE  (Predicted fusion)
    """
    bps = {}
    spatial = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            scaf = row['scaffold'].strip()
            num = int(scaf.replace('scaffold_', ''))
            pos_bp = int(row['breakpoint_pos_bp'].strip())
            pos_mb = pos_bp / 1e6
            cls_raw = row['classification'].strip().upper()
            if 'ARTIFACT' in cls_raw:
                cls = 'ART'
            elif 'FUSION' in cls_raw:
                cls = 'FUSE'
            else:
                print('WARNING: unrecognised classification: {}'.format(
                    cls_raw), file=sys.stderr)
                continue
            left = row['left_ancestor'].strip()
            right = row['right_ancestor'].strip()
            bps.setdefault(num, []).append((pos_mb, cls))
            spatial.setdefault(num, []).append((pos_mb, left, right))
    for num in bps:
        bps[num].sort(key=lambda x: x[0])
    for num in spatial:
        spatial[num].sort(key=lambda x: x[0])
    return bps, spatial


def load_z_candidates(path):
    """Read z_candidate_regions.tsv -> {scaffold_num: (start_mb, end_mb)}"""
    zr = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            scaf = row['scaffold'].strip()
            num = int(scaf.replace('scaffold_', ''))
            s = float(row['start_Mb'])
            e = float(row['end_Mb'])
            zr[num] = (s, e)
    return zr


# ── Build panel data ───────────────────────────────────────────────

def build_ancestry(sizes, compound, compound_aligned, simple, breakpoints, spatial):
    """Merge compound and simple scaffolds into full ancestry map.
    
    For compound scaffolds, segment boundaries and unit ORDER are derived
    from breakpoint positions and left/right allegiance labels.
    For simple scaffolds, the entire scaffold is one unit.
    """
    ancestry = {}
    for n in range(1, 14):
        size_mb = sizes[n] / 1e6
        if n in spatial:
            # Determine spatial order from breakpoint left/right labels
            bp_info = spatial[n]  # [(pos_mb, left_unit, right_unit), ...]
            # Build ordered unit list from breakpoints
            # First unit = leftmost breakpoint's left_unit
            # Between breakpoints = right of previous = left of next
            # Last unit = rightmost breakpoint's right_unit
            boundaries = [0.0]
            units_ordered = [bp_info[0][1]]  # left of first breakpoint
            for pos, left, right in bp_info:
                boundaries.append(pos)
                units_ordered.append(right)
            boundaries.append(size_mb)
            
            segments = []
            for i, de in enumerate(units_ordered):
                segments.append((de, boundaries[i], boundaries[i + 1]))
            ancestry[n] = segments
        elif n in simple:
            de = simple[n]
            ancestry[n] = [(de, 0.0, size_mb)]
        else:
            print('WARNING: scaffold_{} has no ancestry assignment'.format(n),
                  file=sys.stderr)
            ancestry[n] = [('?', 0.0, size_mb)]
    return ancestry


def build_panel_a(sizes, ancestry, breakpoints, z_cands):
    """Build items for panel A (deposited 13-scaffold assembly)."""
    items = []
    for n in range(1, 14):
        size_mb = sizes[n] / 1e6
        units = []
        for de, s, e in ancestry.get(n, []):
            col = DE_COLOURS.get(de, '#999999')
            units.append((de, e - s, col))
        compound = n in breakpoints
        bp_list = []
        for pos, cls in breakpoints.get(n, []):
            bp_list.append(('~{:.0f}'.format(pos), cls))
        zr = None
        if n in z_cands:
            zs, ze = z_cands[n]
            zr = (zs, ze, '{:.0f}\u2013{:.0f}'.format(zs, ze))
        has_art = any(c == 'ART' for _, c in breakpoints.get(n, []))
        has_fuse = any(c == 'FUSE' for _, c in breakpoints.get(n, []))
        items.append(dict(
            label='s{}'.format(n), size=size_mb, compound=compound,
            is_split=has_art, is_fusion=has_fuse,
            units=units, breakpoints=bp_list, z_region=zr))
    items.append(dict(is_w=True))
    return items


def build_panel_b(sizes, ancestry, breakpoints, z_cands):
    """Build items for panel B (inferred 16-unit karyotype).
    Split at artifact breakpoints; keep predicted fusions intact."""
    raw_units = []

    for n in range(1, 14):
        size_mb = sizes[n] / 1e6
        anc_list = ancestry.get(n, [])
        art_positions = [pos for pos, cls in breakpoints.get(n, [])
                         if cls == 'ART']
        fuse_positions = [pos for pos, cls in breakpoints.get(n, [])
                          if cls == 'FUSE']

        if not art_positions:
            has_fuse = len(fuse_positions) > 0
            if has_fuse:
                sub = []
                for de, s, e in anc_list:
                    col = DE_COLOURS.get(de, '#999999')
                    sub.append((de, e - s, col))
                raw_units.append(dict(
                    size=size_mb, src='s{}'.format(n), is_split=False,
                    is_fusion=True, sub=sub, zr=z_cands.get(n)))
            else:
                de = anc_list[0][0] if anc_list else '?'
                col = DE_COLOURS.get(de, '#999999')
                raw_units.append(dict(
                    de=de, col=col, size=size_mb, src='s{}'.format(n),
                    is_split=False, is_fusion=False, sub=None,
                    zr=z_cands.get(n)))
        else:
            cuts = sorted([0.0] + art_positions + [size_mb])
            for i in range(len(cuts) - 1):
                seg_start = cuts[i]
                seg_end = cuts[i + 1]
                seg_size = seg_end - seg_start
                seg_anc = []
                for de, s, e in anc_list:
                    if s < seg_end and e > seg_start:
                        cs = max(s, seg_start)
                        ce = min(e, seg_end)
                        col = DE_COLOURS.get(de, '#999999')
                        seg_anc.append((de, ce - cs, col))
                seg_fuse = any(seg_start < fp < seg_end
                               for fp in fuse_positions)
                zr = None
                if n in z_cands:
                    zs, ze = z_cands[n]
                    if zs < seg_end and ze > seg_start:
                        zr = (max(zs - seg_start, 0),
                              min(ze - seg_start, seg_size))
                if seg_fuse and len(seg_anc) > 1:
                    raw_units.append(dict(
                        size=seg_size, src='s{}'.format(n), is_split=True,
                        is_fusion=True, sub=seg_anc, zr=None))
                else:
                    de = seg_anc[0][0] if seg_anc else '?'
                    col = DE_COLOURS.get(de, '#999999')
                    raw_units.append(dict(
                        de=de, col=col, size=seg_size, src='s{}'.format(n),
                        is_split=True, is_fusion=False, sub=None, zr=zr))

    # Sort by descending size, assign PrChr 1-N
    raw_sorted = sorted(raw_units, key=lambda x: -x['size'])
    items = []
    for i, u in enumerate(raw_sorted):
        nm = 'PrChr {}'.format(i + 1)
        if u.get('sub'):
            units = u['sub']
        else:
            units = [(u.get('de', '?'), u['size'],
                      u.get('col', '#999999'))]
        zr = None
        if u.get('zr') and isinstance(u['zr'], tuple) and len(u['zr']) == 2:
            zs, ze = u['zr']
            zr = (zs, ze, '{:.0f}\u2013{:.0f}'.format(zs, ze))
        items.append(dict(
            label=nm, size=u['size'], is_split=u['is_split'],
            is_fusion=u['is_fusion'], units=units, z_region=zr,
            origin=u['src'], breakpoints=[]))
    items.append(dict(is_w=True))
    return items


# ── Drawing functions ──────────────────────────────────────────────

def draw_z(ax, x0, y, bh, sc, zr):
    zs, ze, zl = zr
    zx0 = x0 + zs * sc
    zx1 = x0 + ze * sc
    ax.add_patch(plt.Rectangle(
        (zx0, y - bh / 2), zx1 - zx0, bh,
        facecolor='none', edgecolor='black', linewidth=0.5,
        hatch='//////', alpha=0.7, zorder=4))
    by = y + bh / 2 + 0.06
    ax.plot([zx0, zx0, zx1, zx1], [by - .03, by, by, by - .03],
            color='black', lw=1.0, zorder=5)
    ax.text((zx0 + zx1) / 2, by + 0.03,
            'Z-cand. ({})'.format(zl),
            ha='center', va='bottom', fontsize=FZ,
            color='black', fontweight='bold', zorder=5)


def draw_bars(ax, items, sc, bh, sp, x0, show_origin=False):
    for i, item in enumerate(items):
        y = (len(items) - 1 - i) * sp

        if item.get('is_w'):
            bar_w = 0.35
            bar_x = x0 + 0.05
            ax.add_patch(plt.Rectangle(
                (bar_x, y - bh / 2), bar_w, bh,
                facecolor='#E8E0F0', edgecolor='#7B1FA2', lw=1.0,
                linestyle=':', hatch='..'))
            ax.text(bar_x - .01, y,
                    'W-cand.',
                    ha='right', va='center', fontsize=FL,
                    color='#7B1FA2', style='italic')
            continue

        # Colour blocks
        cx = x0
        for un, us, uc in item['units']:
            w = us * sc
            ax.add_patch(plt.Rectangle(
                (cx, y - bh / 2), w, bh,
                facecolor=uc, edgecolor='none', alpha=0.85))
            if w > 0.04:
                tc = 'white' if uc not in ['#FDD835', '#BDBDBD'] else 'black'
                ax.text(cx + w / 2, y, un, ha='center', va='center',
                        fontsize=FU, fontweight='bold', color=tc)
            cx += w

        if item.get('z_region'):
            draw_z(ax, x0, y, bh, sc, item['z_region'])

        for bl, bc in item.get('breakpoints', []):
            bm = float(bl.replace('~', ''))
            bx = x0 + bm * sc
            c = '#D32F2F' if bc == 'ART' else '#1B5E20'
            ls = '--' if bc == 'ART' else '-'
            ax.plot([bx, bx], [y - bh / 2 - .06, y + bh / 2 + .06],
                    color=c, lw=1.8, ls=ls, zorder=5)
            sym = '\u2702' if bc == 'ART' else '\u25CF'
            ax.text(bx, y + bh / 2 + .10,
                    '{} {}'.format(sym, bl),
                    ha='center', va='bottom', fontsize=FB,
                    color=c, fontweight='bold')

        lbl_color = 'black'
        lbl_weight = 'normal'
        if item.get('is_fusion'):
            lbl_color = '#1B5E20'; lbl_weight = 'bold'
        elif item.get('compound') or item.get('is_split'):
            lbl_color = '#D32F2F'; lbl_weight = 'bold'
        ax.text(x0 - .01, y,
                '{}  ({:.1f})'.format(item['label'], item['size']),
                ha='right', va='center', fontsize=FL,
                fontweight=lbl_weight, color=lbl_color)

        if show_origin and item.get('origin'):
            ax.text(x0 + item['size'] * sc + .008, y,
                    '\u2190 {}'.format(item['origin']),
                    ha='left', va='center', fontsize=FO,
                    color='gray', style='italic')


# ── Main ───────────────────────────────────────────────────────────

# ── Hardcoded input paths ──────────────────────────────────────────
# All paths relative to package root

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FAI = os.path.join(BASEDIR,
    'input', 'Pristimantis_euphronides.genome.fasta.fai')

COMPOUND_TSV = os.path.join(BASEDIR,
    'output', 'Table_compound_chromosomes.tsv')

KARYOTYPE_TSV = os.path.join(BASEDIR,
    'output', 'Table_karyotype_reconstruction.tsv')

ALLEGIANCE_TSV = os.path.join(BASEDIR,
    'output', 'allegiance_breakpoint_results.tsv')

Z_CANDIDATES = os.path.join(BASEDIR,
    'input', 'z_candidate_regions.tsv')

OUTDIR = os.path.join(BASEDIR, 'output')
DPI = 300


def parse_args():
    p = argparse.ArgumentParser(
        description='Generate karyotype reconstruction figure from '
                    'Tier (ii) + 2 output files.')
    p.add_argument('--outdir', default=OUTDIR,
                   help='Output directory (default: {})'.format(OUTDIR))
    p.add_argument('--dpi', type=int, default=DPI,
                   help='PNG/TIFF resolution (default: {})'.format(DPI))
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load data from actual analysis outputs
    sizes = load_fai(FAI)
    compound, compound_aligned, simple = load_compound_and_karyotype(
        COMPOUND_TSV, KARYOTYPE_TSV)
    breakpoints, spatial = load_allegiance(ALLEGIANCE_TSV)
    z_cands = load_z_candidates(Z_CANDIDATES)

    # Build ancestry from compound + simple + breakpoints + spatial order
    ancestry = build_ancestry(sizes, compound, compound_aligned,
                              simple, breakpoints, spatial)

    # Validate
    for n in range(1, 14):
        if n not in sizes:
            print('ERROR: scaffold_{} missing from .fai'.format(n),
                  file=sys.stderr)
            sys.exit(1)
        if n not in ancestry:
            print('ERROR: scaffold_{} missing from ancestry'.format(n),
                  file=sys.stderr)
            sys.exit(1)

    # Build items
    items_a = build_panel_a(sizes, ancestry, breakpoints, z_cands)
    items_b = build_panel_b(sizes, ancestry, breakpoints, z_cands)

    # Layout
    mx = max(sizes[n] / 1e6 for n in range(1, 14))
    sc = 0.82 / mx
    bh = 0.40
    sp = 0.95
    x0 = 0.12
    n_a = len(items_a)
    n_b = len(items_b)

    fig = plt.figure(figsize=(7.5, 10.5))
    h_a = n_a * sp + 0.5
    h_b = n_b * sp + 0.5
    total = h_a + h_b + 1.0
    frac_a = h_a / total
    frac_lg = 0.5 / total

    # Panel A
    ax_a = fig.add_axes([0.03, 1.0 - frac_a + 0.01, 0.96, frac_a - 0.02])
    ax_a.set_xlim(-0.02, 1.02)
    ax_a.set_ylim(-0.5, n_a * sp + 0.3)
    ax_a.axis('off')
    ax_a.text(-0.02, n_a * sp + 0.25, 'A', fontsize=14,
              fontweight='bold', va='top', ha='left')
    draw_bars(ax_a, items_a, sc, bh, sp, x0,
              show_origin=False)

    # Legend — matches figure: "Predicted scaffolding artifact" / "Predicted fusion"
    lg_y = 1.0 - frac_a - frac_lg
    ax_lg = fig.add_axes([0.03, lg_y, 0.96, frac_lg])
    ax_lg.set_xlim(0, 1); ax_lg.set_ylim(0, 1); ax_lg.axis('off')
    handles = [
        Line2D([0], [0], marker='|', color='#D32F2F', markersize=10,
               markeredgewidth=2, linestyle='None',
               label='Predicted scaffolding artifact'),
        Line2D([0], [0], marker='|', color='#1B5E20', markersize=10,
               markeredgewidth=2, linestyle='None',
               label='Predicted fusion'),
        mpatches.Patch(facecolor='white', edgecolor='black', lw=1.0,
                       hatch='//////', label='Z-candidate'),
    ]
    ax_lg.legend(handles=handles, loc='center', fontsize=FG,
                 frameon=True, ncol=3, handlelength=0.8,
                 columnspacing=1.0)

    # Panel B
    ax_b = fig.add_axes([0.03, 0.01, 0.96, lg_y - 0.02])
    ax_b.set_xlim(-0.02, 1.02)
    ax_b.set_ylim(-0.5, n_b * sp + 0.3)
    ax_b.axis('off')
    ax_b.text(-0.02, n_b * sp + 0.25, 'B', fontsize=14,
              fontweight='bold', va='top', ha='left')
    draw_bars(ax_b, items_b, sc, bh, sp, x0,
              show_origin=True)

    for ext in ['png', 'pdf']:
        out = os.path.join(args.outdir,
                           'Fig_karyotype_schematic.{}'.format(ext))
        kw = {'dpi': args.dpi} if ext == 'png' else {}
        plt.savefig(out, facecolor='white', **kw)
        print('Saved: {}'.format(out))
    plt.close()

    try:
        from PIL import Image
        png = os.path.join(args.outdir, 'Fig_karyotype_schematic.png')
        tif = os.path.join(args.outdir, 'Fig_karyotype_schematic.tiff')
        Image.open(png).save(tif, dpi=(300, 300), compression='tiff_lzw')
        print('Saved: {}'.format(tif))
    except ImportError:
        print('Pillow not available; TIFF not generated.', file=sys.stderr)


if __name__ == '__main__':
    main()
