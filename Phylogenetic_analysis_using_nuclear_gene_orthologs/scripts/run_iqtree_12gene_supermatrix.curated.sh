#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# run_iqtree_12gene_supermatrix.curated.sh
# Run IQ-TREE on the curated 12-gene supermatrix with gene-partitioned model
# selection (MFP+MERGE). 1000 ultrafast bootstrap, 1000 SH-aLRT, --bnni.
#
# Input:  supermatrix_12genes/supermatrix_12genes.curated.fasta
#         supermatrix_12genes/partition.curated.txt
# Output: supermatrix_12genes/iqtree_12genes_curated.*

set -euo pipefail

WORKDIR="supermatrix_12genes"
ALIGNMENT="${WORKDIR}/supermatrix_12genes.curated.fasta"
PARTITION="${WORKDIR}/partition.curated.txt"
PREFIX="${WORKDIR}/iqtree_12genes_curated"

iqtree3 -s "$ALIGNMENT" \
        -p "$PARTITION" \
        -m MFP+MERGE \
        -B 1000 \
        --alrt 1000 \
        --bnni \
        -T AUTO \
        -pre "$PREFIX"

echo "Done: ${PREFIX}.treefile"
