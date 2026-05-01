#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# run_iqtree_5gene_supermatrix.curated.sh
# Run IQ-TREE on the curated 5-gene supermatrix with gene-partitioned model
# selection (MFP+MERGE). 1000 ultrafast bootstrap, 1000 SH-aLRT, --bnni.
#
# Input:  supermatrix_5genes/supermatrix.5genes.curated.fasta
#         supermatrix_5genes/partition.curated.txt
# Output: supermatrix_5genes/iqtree_5genes_curated.*

set -euo pipefail

WORKDIR="supermatrix_5genes"
ALIGNMENT="${WORKDIR}/supermatrix.5genes.curated.fasta"
PARTITION="${WORKDIR}/partition.curated.txt"
PREFIX="${WORKDIR}/iqtree_5genes_curated"

iqtree3 -s "$ALIGNMENT" \
        -p "$PARTITION" \
        -m MFP+MERGE \
        -B 1000 \
        --alrt 1000 \
        --bnni \
        -T AUTO \
        -pre "$PREFIX"

echo "Done: ${PREFIX}.treefile"
