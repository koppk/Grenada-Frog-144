#!/usr/bin/env python3
"""
w_ideogram_schmid.py

Supplementary Figure: ZW chromosome properties of Pristimantis euphronides
===========================================================================
Reconstructed from Schmid et al. (2002), Cytogenet Genome Res 97:81-94.

Panel A: Grayscale C-banding ideogram (standard textbook convention).
Panel B: Fluorochrome response per chromosomal region, using representative
         emission colors. Schematic — zone-level only.

Z left, W right, following ZW convention and Schmid's presentation.

Author: Kopp K., P. euphronides genome project
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

OUTDIR = "ZW-Schmid_2002-ideogram-Figure"

# ---------------------------------------------------------------------------
# Zone model
# ---------------------------------------------------------------------------
W_ZONES = [
    (0.000, 0.040, 'tip_p',       'tip'),
    (0.040, 0.200, 'proximal_p',  'p-arm'),
    (0.200, 0.235, 'centromere',   'CEN'),
    (0.235, 0.600, 'proximal_q',  'prox. q'),
    (0.600, 0.960, 'distal_q',    'dist. q'),
    (0.960, 1.000, 'tip_q',       'tip'),
]

Z_ZONES = [
    (0.000, 0.100, 'z_cen',  'CEN'),
    (0.100, 0.850, 'z_eu',   'eu'),
    (0.850, 1.000, 'z_tel',  'tel'),
]

# ---------------------------------------------------------------------------
# Panel A: C-banding grayscale
# ---------------------------------------------------------------------------
ZONE_GRAY = {
    'tip_p':       '#F0F0F0',
    'proximal_p':  '#505050',
    'centromere':  '#1A1A1A',
    'proximal_q':  '#404040',
    'distal_q':    '#606060',
    'tip_q':       '#F0F0F0',
    'z_cen':       '#707070',
    'z_eu':        '#F0F0F0',
    'z_tel':       '#505050',
}

# ---------------------------------------------------------------------------
# Panel B: Fluorochrome colors
# ---------------------------------------------------------------------------
Q_BRIGHT = (0.75, 0.92, 0.12)
Q_DIM    = (0.18, 0.25, 0.04)
Q_MOD    = (0.40, 0.50, 0.06)

H_BRIGHT = (0.40, 0.82, 0.90)
H_DIM    = (0.10, 0.25, 0.28)
H_MOD    = (0.20, 0.45, 0.50)

M_BRIGHT = (0.92, 0.78, 0.12)
M_DIM    = (0.22, 0.18, 0.04)
M_MOD    = (0.50, 0.42, 0.06)

D_BRIGHT = (0.25, 0.40, 0.92)
D_DIM    = (0.06, 0.09, 0.25)
D_MOD    = (0.12, 0.18, 0.45)

FLUORO_COLOR_MAP = {
    'Q': {'bright': Q_BRIGHT, 'mod': Q_MOD, 'dim': Q_DIM},
    'H': {'bright': H_BRIGHT, 'mod': H_MOD, 'dim': H_DIM},
    'M': {'bright': M_BRIGHT, 'mod': M_MOD, 'dim': M_DIM},
    'D': {'bright': D_BRIGHT, 'mod': D_MOD, 'dim': D_DIM},
}

FLUORO_KEYS = ['Q', 'H', 'M', 'D']

FLUORO_W = {
    'tip_p':        ('dim',    'dim',    'dim',    'dim'),
    'proximal_p':   ('bright', 'bright', 'dim',    'dim'),
    'centromere':   ('dim',    'bright', 'dim',    'bright'),
    'proximal_q':   ('mixed',  'bright', 'mixed',  'dim'),
    'distal_q':     ('mixed',  'mod',    'mixed',  'dim'),
    'tip_q':        ('dim',    'dim',    'dim',    'dim'),
}

FLUORO_Z = {
    'z_cen':  ('dim', 'dim', 'dim', 'dim'),
    'z_eu':   ('dim', 'dim', 'dim', 'dim'),
    'z_tel':  ('dim', 'dim', 'bright', 'dim'),
}


# ---------------------------------------------------------------------------
# Drawing functions
# ---------------------------------------------------------------------------

def draw_zone_chromosome(ax, zones, x_center, width, y_bottom, y_top,
                          color_map, label_below, label_fs=10):
    total_h = y_top - y_bottom
    for (s, e, zid, _) in zones:
        ys = y_top - s * total_h
        ye = y_top - e * total_h
        color = color_map.get(zid, '#CCCCCC')
        ax.add_patch(plt.Rectangle((x_center - width/2, ye), width, ys - ye,
                                    facecolor=color, edgecolor='none'))
    ax.add_patch(FancyBboxPatch(
        (x_center - width/2, y_bottom), width, total_h,
        boxstyle=f"round,pad={width*0.15}",
        linewidth=0.8, edgecolor='black', facecolor='none', zorder=10))
    for (s, e, zid, _) in zones:
        if 'centromere' in zid or zid == 'z_cen':
            cy = y_top - ((s + e) / 2) * total_h
            nw = width * 0.15
            ax.plot([x_center - width/2 - 0.005, x_center - width/2 + nw],
                    [cy, cy], color='black', lw=1.0, zorder=11)
            ax.plot([x_center + width/2 + 0.005, x_center + width/2 - nw],
                    [cy, cy], color='black', lw=1.0, zorder=11)
    ax.text(x_center, y_bottom - 0.03, label_below, fontsize=label_fs,
            ha='center', va='top', fontweight='bold')


def draw_fluoro_zone_chr(ax, zones, response_map, x_center, width,
                          y_bottom, y_top, fluoro_key):
    total_h = y_top - y_bottom
    cmap = FLUORO_COLOR_MAP[fluoro_key]
    fi = FLUORO_KEYS.index(fluoro_key)

    for (s, e, zid, _) in zones:
        ys = y_top - s * total_h
        ye = y_top - e * total_h
        zh = ys - ye
        resp = response_map.get(zid, ('dim',)*4)[fi]

        if resp == 'mixed':
            bright_c = cmap['bright']
            dim_c = cmap['dim']
            n_stripes = max(4, int(zh / 0.012))
            for si in range(n_stripes):
                sy = ye + (si / n_stripes) * zh
                sh = zh / n_stripes
                c = bright_c if si % 2 == 0 else dim_c
                ax.add_patch(plt.Rectangle((x_center - width/2, sy), width, sh,
                                            facecolor=c, edgecolor='none'))
        else:
            color = cmap.get(resp, cmap['dim'])
            ax.add_patch(plt.Rectangle((x_center - width/2, ye), width, zh,
                                        facecolor=color, edgecolor='none'))

    ax.add_patch(FancyBboxPatch(
        (x_center - width/2, y_bottom), width, total_h,
        boxstyle=f"round,pad={width*0.12}",
        linewidth=0.5, edgecolor='#555555', facecolor='none', zorder=10))
    for (s, e, zid, _) in zones:
        if 'centromere' in zid or zid == 'z_cen':
            cy = y_top - ((s + e) / 2) * total_h
            nw = width * 0.12
            ax.plot([x_center - width/2 - 0.003, x_center - width/2 + nw],
                    [cy, cy], color='#AAAAAA', lw=0.6, zorder=11)
            ax.plot([x_center + width/2 + 0.003, x_center + width/2 - nw],
                    [cy, cy], color='#AAAAAA', lw=0.6, zorder=11)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def make_figure():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
    })

    FIG_W_MM = 170
    FIG_H_MM = 165   # slightly taller for legend room
    fig = plt.figure(figsize=(FIG_W_MM/25.4, FIG_H_MM/25.4))

    ax_a = fig.add_axes([0.02, 0.18, 0.34, 0.76])
    ax_b = fig.add_axes([0.40, 0.18, 0.58, 0.76])

    w_y_top = 0.90
    w_y_bot = 0.09
    total_w = w_y_top - w_y_bot

    z_frac = 1.4 / 5.5
    z_h = z_frac * total_w
    z_y_top = w_y_top
    z_y_bot = z_y_top - z_h

    # =================== Panel A ===================
    ax_a.text(-0.02, 1.01, 'A', fontsize=13, fontweight='bold',
              transform=ax_a.transAxes, va='bottom', ha='left')

    z_x_a = 0.22;  z_w_a = 0.06
    w_x_a = 0.62;  w_w_a = 0.13

    draw_zone_chromosome(ax_a, Z_ZONES, z_x_a, z_w_a, z_y_bot, z_y_top,
                          ZONE_GRAY, 'Z', label_fs=10)
    draw_zone_chromosome(ax_a, W_ZONES, w_x_a, w_w_a, w_y_bot, w_y_top,
                          ZONE_GRAY, 'W', label_fs=10)

    ax_a.text(z_x_a, z_y_top + 0.03, '1.4 \u00b5m', fontsize=6.5,
              ha='center', color='#333333')
    ax_a.text(w_x_a, w_y_top + 0.03, '5.5 \u00b5m', fontsize=6.5,
              ha='center', color='#333333')

    cen_y = w_y_top - 0.2175 * total_w
    ax_a.text(w_x_a - w_w_a/2 - 0.07, (w_y_top + cen_y)/2, 'p',
              fontsize=9, ha='center', va='center', style='italic', color='#555555')
    ax_a.text(w_x_a - w_w_a/2 - 0.07, (cen_y + w_y_bot)/2, 'q',
              fontsize=9, ha='center', va='center', style='italic', color='#555555')

    um_per_unit = 5.5 / total_w
    bar_h = 1.0 / um_per_unit
    bar_y = w_y_bot + 0.01
    bar_x = 0.92
    ax_a.plot([bar_x, bar_x], [bar_y, bar_y + bar_h],
              color='black', linewidth=1.5, solid_capstyle='butt')
    ax_a.text(bar_x + 0.02, bar_y + bar_h/2, '1 \u00b5m', fontsize=6,
              va='center', ha='left')

    annot_x0 = w_x_a + w_w_a/2 + 0.015
    annotations = [
        (0.120, 'Q\u207a band (i)'),
        (0.2175, 'H\u207a\u207a, DAPI\u207a'),
        (0.418, 'Q\u207a (ii, iii) alt.\nwith M\u207a'),
        (0.780, 'Q\u207a (iv); weaker\nH, DA/D'),
    ]
    for (frac, text) in annotations:
        y = w_y_top - frac * total_w
        ax_a.annotate(text, xy=(annot_x0, y), xytext=(0.99, y),
                      fontsize=5, va='center', ha='right', color='#333333',
                      arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=0.5))

    ax_a.set_xlim(0, 1);  ax_a.set_ylim(0, 1)
    ax_a.set_xticks([]);  ax_a.set_yticks([])
    ax_a.spines[:].set_visible(False)

    # =================== Panel B ===================
    ax_b.text(-0.03, 1.01, 'B', fontsize=13, fontweight='bold',
              transform=ax_b.transAxes, va='bottom', ha='left')
    ax_b.set_facecolor('black')

    # More generous spacing: each column gets 0.24 of axes width
    # Within each column: Z at offset 0.0, W at offset 0.08
    col_spacing = 0.24
    z_w_b = 0.025
    w_w_b = 0.055

    for fi in range(4):
        fkey = FLUORO_KEYS[fi]
        x_base = 0.05 + fi * col_spacing

        z_x = x_base
        w_x = x_base + 0.065

        draw_fluoro_zone_chr(ax_b, Z_ZONES, FLUORO_Z,
                              z_x, z_w_b, z_y_bot, z_y_top, fkey)
        draw_fluoro_zone_chr(ax_b, W_ZONES, FLUORO_W,
                              w_x, w_w_b, w_y_bot, w_y_top, fkey)

        ax_b.text(z_x, z_y_bot - 0.035, 'Z', fontsize=8, ha='center',
                  fontweight='bold', color='white')
        ax_b.text(w_x, w_y_bot - 0.035, 'W', fontsize=9, ha='center',
                  fontweight='bold', color='white')

        # Header — consistent single-letter abbreviation
        header_x = (z_x + w_x) / 2
        ax_b.text(header_x, 0.97, FLUORO_KEYS[fi],
                  fontsize=9, ha='center', va='bottom', fontweight='bold',
                  color='white')

    ax_b.set_xlim(-0.02, 1.0)
    ax_b.set_ylim(0, 1.02)
    ax_b.set_xticks([]);  ax_b.set_yticks([])
    ax_b.spines[:].set_visible(False)

    # =================== Legend A ===================
    ax_la = fig.add_axes([0.02, 0.02, 0.34, 0.12])
    ax_la.set_xlim(0, 1);  ax_la.set_ylim(0, 1)
    ax_la.axis('off')

    a_items = [
        ('#1A1A1A', 'hetero-\nchromatic'),
        ('#606060', 'het. (less\nintense)'),
        ('#F0F0F0', 'euchro-\nmatic'),
    ]
    for i, (color, label) in enumerate(a_items):
        x = 0.02 + i * 0.34
        ax_la.add_patch(plt.Rectangle((x, 0.35), 0.06, 0.40,
                         facecolor=color, edgecolor='black', linewidth=0.4))
        ax_la.text(x + 0.10, 0.55, label, fontsize=5, va='center',
                   linespacing=1.1)

    # =================== Legend B ===================
    # Align legend columns directly under Panel B fluorochrome columns.
    # Panel B header_x positions (in ax_b axes coords):
    #   Q: 0.0825, H: 0.3225, M: 0.5625, D: 0.8025
    # ax_b and ax_lb share the same figure x range, so same axes coords align.

    ax_lb = fig.add_axes([0.40, 0.02, 0.58, 0.12])
    ax_lb.set_xlim(-0.02, 1.0)   # match ax_b xlim
    ax_lb.set_ylim(0, 1)
    ax_lb.axis('off')

    fluoro_legend = [
        ('Q', Q_BRIGHT, Q_DIM, True),   # has mixed
        ('H', H_BRIGHT, H_DIM, False),
        ('M', M_BRIGHT, M_DIM, True),   # has mixed
        ('D', D_BRIGHT, D_DIM, False),
    ]

    # Column centers — shifted slightly right to align visually under the
    # Z/W pair (W is wider, so visual center sits right of arithmetic midpoint)
    header_xs = [0.0975, 0.3375, 0.5775, 0.8175]

    box_w = 0.028
    box_h = 0.30
    gap = 0.006
    pair_w = 2 * box_w + gap  # total width of bright+dim pair

    row1_y = 0.55

    for i, (label, bright, dim, has_mixed) in enumerate(fluoro_legend):
        cx = header_xs[i]
        x0 = cx - pair_w / 2  # left edge of bright box

        # Bright box
        ax_lb.add_patch(plt.Rectangle((x0, row1_y), box_w, box_h,
                         facecolor=bright, edgecolor='#444444', linewidth=0.4))
        # Dim box
        ax_lb.add_patch(plt.Rectangle((x0 + box_w + gap, row1_y), box_w, box_h,
                         facecolor=dim, edgecolor='#444444', linewidth=0.4))
        # Label
        ax_lb.text(x0 + pair_w + 0.012, row1_y + box_h / 2,
                   f'{label}: bright\n     / dim', fontsize=4.5, va='center',
                   color='black', linespacing=1.1)

        # Mixed stripe box (row 2) — only for Q and M
        if has_mixed:
            row2_y = 0.10
            mx0 = cx - box_w / 2  # center single box under column
            n_s = 4
            stripe_h = box_h / n_s
            for si in range(n_s):
                sy = row2_y + si * stripe_h
                c = bright if si % 2 == 0 else dim
                ax_lb.add_patch(plt.Rectangle((mx0, sy), box_w, stripe_h,
                                 facecolor=c, edgecolor='none'))
            ax_lb.add_patch(plt.Rectangle((mx0, row2_y), box_w, box_h,
                             facecolor='none', edgecolor='#444444', linewidth=0.4))
            ax_lb.text(mx0 + box_w + 0.012, row2_y + box_h / 2,
                       'mixed (alternating\nbright/dim in zone)',
                       fontsize=4.5, va='center', color='black',
                       linespacing=1.1)

    # =================== Save ===================
    os.makedirs(OUTDIR, exist_ok=True)
    DPI = 300
    for ext in ['png', 'pdf']:
        out = os.path.join(OUTDIR, f'W_chromosome_ideogram_Schmid2002.{ext}')
        fig.savefig(out, dpi=DPI, bbox_inches='tight', pad_inches=0.03,
                    facecolor='white')
        print(f"Saved: {out}")

    out_tiff = os.path.join(OUTDIR, 'W_chromosome_ideogram_Schmid2002.tiff')
    fig.savefig(out_tiff, dpi=DPI, bbox_inches='tight', pad_inches=0.03,
                facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
    print(f"Saved: {out_tiff}")

    from PIL import Image
    out_png = os.path.join(OUTDIR, 'W_chromosome_ideogram_Schmid2002.png')
    img = Image.open(out_png)
    w_px, h_px = img.size
    print(f"\nDimensions: {w_px}x{h_px} px = {w_px/DPI*25.4:.1f}x{h_px/DPI*25.4:.1f} mm at {DPI} DPI")
    print(f"TIFF: {os.path.getsize(out_tiff)/1024/1024:.2f} MB")
    plt.close()


if __name__ == '__main__':
    make_figure()
