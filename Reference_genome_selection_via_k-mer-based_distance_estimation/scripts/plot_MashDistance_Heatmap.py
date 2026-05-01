# Author: Kopp K, Pristimantis euphronides genome project
#
# Generate a lower-triangle heatmap of Mash distances from mash triangle output.
#
# Input:  ref_triangle.tab  (mash triangle output, relaxed Phylip format)
#         label_map.tsv     (two columns: filename_substring → display_label,
#                            listed in desired sort order)
# Output: <input_basename>_heatmap.tiff, <input_basename>_heatmap.png (300 dpi)
#
# Usage:  python3 plot_MashDistance_Heatmap.py ref_triangle.tab label_map.tsv
#
# Requirements: Python 3.x, numpy, matplotlib, seaborn

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

# --- Parse arguments ---
input_file = sys.argv[1] if len(sys.argv) > 1 else "ref_triangle.tab"
label_file = sys.argv[2] if len(sys.argv) > 2 else "label_map.tsv"

# --- Parse relaxed Phylip lower-triangle format ---
# Handles both with and without the leading count line
with open(input_file) as f:
    lines = [line.rstrip('\n') for line in f if line.strip()]

# Detect whether first line is the count or the first label
try:
    n = int(lines[0].strip())
    data_start = 1
except ValueError:
    # No count header; infer n from number of lines
    n = len(lines)
    data_start = 0

raw_labels = []
matrix = np.zeros((n, n))

for i in range(n):
    parts = lines[data_start + i].split('\t')
    raw_labels.append(parts[0])
    for j, val in enumerate(parts[1:]):
        if val:
            matrix[i][j] = float(val)
            matrix[j][i] = float(val)

# --- Read label map (substring -> display label, in sort order) ---
substrings = []
display_labels = []
with open(label_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        key, label = line.split('\t', 1)
        substrings.append(key)
        display_labels.append(label)

# --- Match raw labels to display labels and determine sort order ---
raw_to_display = {}
raw_to_order = {}
for raw in raw_labels:
    matched = False
    for idx, sub in enumerate(substrings):
        if sub in raw:
            raw_to_display[raw] = display_labels[idx]
            raw_to_order[raw] = idx
            matched = True
            break
    if not matched:
        print(f"WARNING: no match for '{raw}' in label map")
        raw_to_display[raw] = raw
        raw_to_order[raw] = len(substrings) + len(raw_to_order)

# --- Reorder matrix by label map order ---
sorted_indices = sorted(range(n), key=lambda i: raw_to_order[raw_labels[i]])
matrix = matrix[np.ix_(sorted_indices, sorted_indices)]
row_labels = [raw_to_display[raw_labels[i]] for i in sorted_indices]
col_labels = ['|'.join(l.split('|')[1:]) for l in row_labels]

# --- Mask upper triangle ---
mask_upper = np.triu(np.ones_like(matrix, dtype=bool))

# --- Colour map: dark green to yellow to dark red ---
green_red_cmap = LinearSegmentedColormap.from_list(
    "green_red", ["#006400", "#FFFF00", "#e60000"]
)

# --- Plot ---
size = max(8, n * 0.49)
fig, ax = plt.subplots(figsize=(size, size), dpi=300)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.3)

sns.heatmap(
    matrix,
    mask=mask_upper,
    cmap=green_red_cmap,
    linewidths=0.3,
    linecolor='gray',
    ax=ax,
    square=True,
    xticklabels=col_labels,
    yticklabels=row_labels,
    cbar_ax=cax,
    cbar_kws={'label': 'Mash Distance'}
)

ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=6)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
ax.tick_params(axis='x', direction='out')
ax.tick_params(axis='y', direction='out')
cax.tick_params(direction='out', labelsize=6)

plt.subplots_adjust(left=0.38, right=0.91, top=0.97, bottom=0.08)

# --- Export ---
basename = input_file.rsplit('.', 1)[0]
for fmt in ['tiff', 'png']:
    outfile = f"{basename}_heatmap.{fmt}"
    plt.savefig(outfile, dpi=300, bbox_inches='tight', format=fmt)
    print(f"Heatmap saved to: {outfile}")
plt.close()
