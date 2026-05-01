#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# full_alignment_pipeline.sh
# Per-gene pipeline: MAFFT v7.475 initial alignment, terminal gap trimming,
# MAFFT high-accuracy re-alignment (--globalpair --maxiterate 1000).
#
# Input:  concatenated/<GENE>_PriEup_Hyloidea_RanTemp.fasta
# Output: aligned/<GENE>_aligned.fasta
#         trimmed/<GENE>_trimmed.fasta
#         refined/<GENE>_mafft.fasta

set -euo pipefail

INPUT_DIR="concatenated"
ALIGN_DIR="aligned"
TRIM_DIR="trimmed"
FINAL_DIR="refined"
BLUNT_TRIM_SCRIPT="./blunt_trim_alignment.sh"

mkdir -p "$ALIGN_DIR" "$TRIM_DIR" "$FINAL_DIR"

for fasta in "$INPUT_DIR"/*.fasta; do
    gene=$(basename "$fasta" .fasta)
    aln_out="${ALIGN_DIR}/${gene}_aligned.fasta"
    trim_out="${TRIM_DIR}/${gene}_trimmed.fasta"
    final_out="${FINAL_DIR}/${gene}_mafft.fasta"

    echo "$gene: initial MAFFT alignment"
    mafft --auto "$fasta" > "$aln_out"

    echo "$gene: trimming terminal gaps"
    "$BLUNT_TRIM_SCRIPT" "$aln_out" "$trim_out"

    echo "$gene: high-accuracy MAFFT realignment"
    mafft --maxiterate 1000 --globalpair "$trim_out" > "$final_out"

    echo "$gene: done"
done
