#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# build_supermatrix_clean.sh
# Construct 12-gene concatenated supermatrix with partition file from all
# 12 retained gene alignments. Missing genes are filled with gap characters.
#
# Input:  renamed/<GENE>_renamed.fasta
# Output: supermatrix/supermatrix.fasta
#         supermatrix/partition.txt
#         supermatrix/all_taxa.txt

set -euo pipefail

INPUT_DIR="renamed"
OUTPUT_DIR="supermatrix"

mkdir -p "$OUTPUT_DIR"

SUPERMATRIX="${OUTPUT_DIR}/supermatrix.fasta"
PARTITIONS="${OUTPUT_DIR}/partition.txt"
TAXA_LIST="${OUTPUT_DIR}/all_taxa.txt"

> "$SUPERMATRIX"
> "$PARTITIONS"
> "$TAXA_LIST"

for file in "$INPUT_DIR"/*.fasta; do
    grep "^>" "$file" | sed 's/^>//; s/ *|.*//'
done | sort -u > "$TAXA_LIST"

declare -A SUPERSEQ
while IFS= read -r taxon; do
    SUPERSEQ["$taxon"]=""
done < "$TAXA_LIST"

TOTAL_START=1

for file in "$INPUT_DIR"/*.fasta; do
    gene=$(basename "$file" _renamed.fasta)
    echo "$gene: adding to supermatrix"

    declare -A GENESEQ
    ALIGN_LEN=0

    awk '
    BEGIN { RS = ">"; FS = "\n" }
    NR > 1 {
        header = $1
        gsub(/ .*/, "", header)
        gsub(/ *\|.*/, "", header)
        seq = ""
        for (i = 2; i <= NF; i++) seq = seq $i
        print header "\t" seq
    }' "$file" > "$OUTPUT_DIR/${gene}_temp.tsv"

    while IFS=$'\t' read -r name seq; do
        GENESEQ["$name"]="$seq"
        ALIGN_LEN=${#seq}
    done < "$OUTPUT_DIR/${gene}_temp.tsv"

    for taxon in "${!SUPERSEQ[@]}"; do
        if [[ -n "${GENESEQ[$taxon]:-}" ]]; then
            SUPERSEQ["$taxon"]+="${GENESEQ[$taxon]}"
        else
            SUPERSEQ["$taxon"]+=$(printf '%*s' "$ALIGN_LEN" | tr ' ' '-')
        fi
    done

    END=$((TOTAL_START + ALIGN_LEN - 1))
    echo "DNA, $gene = $TOTAL_START-$END" >> "$PARTITIONS"
    TOTAL_START=$((END + 1))

    unset GENESEQ
    rm -f "$OUTPUT_DIR/${gene}_temp.tsv"
done

while IFS= read -r taxon; do
    echo ">$taxon" >> "$SUPERMATRIX"
    echo "${SUPERSEQ[$taxon]}" >> "$SUPERMATRIX"
done < "$TAXA_LIST"

echo "Supermatrix: $SUPERMATRIX"
echo "Partition: $PARTITIONS"
