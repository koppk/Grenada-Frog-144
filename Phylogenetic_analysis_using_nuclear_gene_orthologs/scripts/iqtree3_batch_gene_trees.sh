#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# iqtree3_batch_gene_trees.sh
# Run IQ-TREE v2.2.0.3 on each of the 12 individual gene alignments.
# Model selection via ModelFinder (MFP), 1000 ultrafast bootstrap replicates,
# 1000 SH-aLRT replicates, with --bnni for bootstrap tree refinement.
#
# Input:  refined/<GENE>_mafft.fasta
# Output: iqtree_single_genes/<GENE>.treefile, .iqtree, .contree, .ufboot, etc.

set -euo pipefail

INPUT_DIR="refined"
OUT_DIR="iqtree_single_genes"

mkdir -p "$OUT_DIR"

for aln in "$INPUT_DIR"/*_mafft.fasta; do
    gene=$(basename "$aln" _mafft.fasta)
    out_prefix="${OUT_DIR}/${gene}"

    echo "$gene: running IQ-TREE"
    iqtree3 -s "$aln" \
            -m MFP \
            -B 1000 \
            --alrt 1000 \
            --bnni \
            -T AUTO \
            -pre "$out_prefix"
    echo "$gene: done"
done
