#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# final_refined_mafft.sh
# Final MAFFT re-alignment of curated alignments after Jalview editing.
# Uses high-accuracy settings (--globalpair --maxiterate 1000).
#
# Input:  trimmed/<GENE>_trimmed.fasta
# Output: refined/<GENE>_mafft.fasta

set -euo pipefail

INPUT_DIR="trimmed"
OUTPUT_DIR="refined"

mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*_trimmed.fasta; do
    gene=$(basename "$file" _trimmed.fasta)
    output="${OUTPUT_DIR}/${gene}_mafft.fasta"

    echo "$gene: high-accuracy MAFFT alignment"
    mafft --globalpair --maxiterate 1000 --thread -1 "$file" > "$output"
    echo "$gene: done"
done
