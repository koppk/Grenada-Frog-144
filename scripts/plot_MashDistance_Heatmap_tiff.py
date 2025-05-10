"""
Generate a 24x24 lower triangle heatmap of Mash distances for Hyloidea genomes.

Requirements:
- Python 3.x
- pandas, numpy, matplotlib, seaborn, odfpy 
- might have to install with: pip install pandas numpy matplotlib seaborn odfpy
- Input: .ods file with 24x24 numeric matrix (lower triangle) and row/column labels in column 0
- Output: high-resolution TIFF suitable for LibreOffice and BMC Genomics supplementary material

Author: [Your Name]
Date: [YYYY-MM-DD]
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# -----------------------------
# Step 1: Load the input matrix
# -----------------------------
# Input file: ODS with labels in column 0, and a 24x24 numeric matrix (lower triangle) from columns 1–24
input_file = "ref_triangle_sorted_lower_triangle.ods"  # Make sure this file is in the same directory

df = pd.read_excel(input_file, engine="odf", header=None)

# Extract the 24 labels from column 0 (rows 1–24)
labels_24 = df.iloc[1:25, 0].tolist()

# Extract the 24x24 matrix (rows 1–24, columns 1–24)
matrix = df.iloc[1:25, 1:25]
matrix = matrix.apply(pd.to_numeric, errors='coerce')

# Set index and column labels:
# - Row labels: full label (Family|Species|AssemblyName)
# - Column labels: trimmed to Species|AssemblyName
row_labels = labels_24
col_labels = ['|'.join(label.split('|')[1:]) for label in labels_24]

matrix.index = row_labels
matrix.columns = col_labels

# -----------------------------------------------------
# Step 2: Create a mask to display only the lower half
# -----------------------------------------------------
mask_upper = np.triu(np.ones_like(matrix, dtype=bool))

# -----------------------------------
# Step 3: Define the color map scheme
# -----------------------------------
green_red_cmap = LinearSegmentedColormap.from_list(
    "green_red", ["#006400", "#FFFF00", "#e60000"]
)

# --------------------------------------------
# Step 4: Plot and configure the heatmap layout
# --------------------------------------------
fig, ax = plt.subplots(figsize=(11.69, 11.69), dpi=300)  # A4 portrait square

# Create external axis for the colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.3)

# Plot the heatmap
sns.heatmap(
    matrix,
    mask=mask_upper,
    cmap=green_red_cmap,
    linewidths=0.3,
    linecolor='gray',
    ax=ax,
    square=True,
    xticklabels=True,
    yticklabels=True,
    cbar_ax=cax,
    cbar_kws={'label': 'Mash Distance'}
)

# Format the tick labels
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=6)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
ax.tick_params(axis='x', direction='out')
ax.tick_params(axis='y', direction='out')
cax.tick_params(direction='out', labelsize=6)

# Adjust layout to avoid clipping labels
plt.subplots_adjust(left=0.38, right=0.91, top=0.97, bottom=0.08)

# ---------------------------------
# Step 5: Export to high-quality TIFF
# ---------------------------------
output_file = "MashDistance_Heatmap_24x24_A4_300dpi_FINAL.tiff"
plt.savefig(output_file, dpi=300, bbox_inches='tight', format='tiff')
plt.close()

print(f"✅ Heatmap successfully saved to: {output_file}")

