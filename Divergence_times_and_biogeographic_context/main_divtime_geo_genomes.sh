#!/bin/bash
# main_divtime_geo_genomes.sh
#
# Reference genome selection and biogeographic context analysis for P. euphronides.
#
# Identifies the phylogenetically closest chromosome-level reference genome for
# scaffolding based on divergence time (tMRCA), quantifies the data gap (species
# closer than best reference but lacking genomes), and adds geographic context.
# Author: Kopp K, Pristimantis euphronides genome project

set -euo pipefail

FOCAL="Pristimantis_euphronides"
OUTDIR="output"

# clean output directory
if [ -d "$OUTDIR" ]; then
    rm -rf "$OUTDIR"
fi
mkdir -p "$OUTDIR"

# Step 1: compute tMRCA and topological distance from focal to all tree tips,
#         flag species with chromosome-level genomes
python scripts/build_tree_metrics_chromlevel.py \
  --tree input/Supplementary_File_S3_time_tree.tre \
  --chromosome_level_list input/Anura_chromosome-level_genomes_unique.tsv \
  --focal "$FOCAL" \
  --out "$OUTDIR/tmrca_topo_chromlevel.tsv"

# Step 2: add biogeographic context - ISO country codes from AmphibiaWeb,
#         minimum distance to Grenada (haversine), Caribbean/mainland flags
python scripts/add_isocc_and_distances.py \
  --tmrca_tsv "$OUTDIR/tmrca_topo_chromlevel.tsv" \
  --amphibia_names input/amphib_names.txt \
  --countries_csv input/countries.csv \
  --out "$OUTDIR/tmrca_topo_chromlevel_with_geo.tsv"

# Step 3: generate summary - rank chr-level genomes, report phylogenetic gap
python scripts/create_reference_selection_summary.py \
  --input "$OUTDIR/tmrca_topo_chromlevel_with_geo.tsv" \
  --focal "$FOCAL" \
  --outdir "$OUTDIR"

# Step 4: refine biogeographic context for the 12 closest mainland Pristimantis
#         congeners - extract per-species type-locality coordinates either from
#         the original taxonomic description (where given as DMS or decimal) or
#         from manual Google Maps lookups recorded in the input file, then
#         compute haversine distances to Pointe Salines (the southwesternmost
#         point of Grenada and the closest single point of the island to the
#         South American mainland)
python scripts/compute_typelocality_distances.py \
  --asw_with_googlecoords input/12_Pristimantis_amphibiansoftheworld_WithGoogleCoords.txt \
  --out "$OUTDIR/12_Pristimantis_typelocality_distances.tsv"

echo "Done. Results in $OUTDIR/"
