#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# build_5gene_supermatrix.sh
# Construct 5-gene concatenated supermatrix (CXCR4, POMC, SIAH1, SLC8A3, TYR)
# with partition file from selected.families.gene.set.
#
# Input:  acc_taxid_metadata/selected.families.gene.set
#         refined/<GENE>_mafft.renamed.fasta
# Output: supermatrix_5genes/supermatrix.fasta
#         supermatrix_5genes/partition.txt

set -euo pipefail

TABLE="./acc_taxid_metadata/selected.families.gene.set"
ALIGN_DIR="./refined"
OUT_DIR="./supermatrix_5genes"
OUT_FASTA="${OUT_DIR}/supermatrix.fasta"
PARTITION_FILE="${OUT_DIR}/partition.txt"
GENES=("SLC8A3" "POMC" "TYR" "SIAH1" "CXCR4")

mkdir -p "$OUT_DIR"
> "$OUT_FASTA"
> "$PARTITION_FILE"

IFS=$'\t' read -r -a HEADERS < <(head -n 1 "$TABLE")
declare -A COL_IDX
for i in "${!HEADERS[@]}"; do
    COL_IDX[${HEADERS[$i]}]=$((i + 1))
done

total_start=1
first_species=""

tail -n +2 "$TABLE" | while IFS=$'\t' read -r -a FIELDS; do
    species="${FIELDS[0]}"
    species_id="${species// /_}"
    concat_seq=""

    [[ -z "$first_species" ]] && first_species="$species_id"

    for gene in "${GENES[@]}"; do
        acc="${FIELDS[${COL_IDX[$gene]}-1]}"
        fasta_file="${ALIGN_DIR}/${gene}_mafft.renamed.fasta"

        if [[ ! -f "$fasta_file" ]]; then
            echo "Missing alignment: $fasta_file" >&2
            exit 1
        fi

        seq=$(awk -v acc=">$acc" '
            $0 ~ "^"acc"$" {getline; print; found=1; next}
            found && $0 ~ "^>" {exit}
            found {printf "%s", $0}
        ' "$fasta_file")

        if [[ -z "$seq" ]]; then
            echo "Missing sequence: $acc in $gene" >&2
            exit 1
        fi

        seq_len=${#seq}
        concat_seq+="$seq"

        if [[ "$species_id" == "$first_species" ]]; then
            start=$total_start
            end=$((start + seq_len - 1))
            echo "DNA, $gene = $start-$end" >> "$PARTITION_FILE"
            total_start=$((end + 1))
        fi
    done

    echo ">$species_id" >> "$OUT_FASTA"
    echo "$concat_seq" >> "$OUT_FASTA"
done

echo "Supermatrix: $OUT_FASTA"
echo "Partition: $PARTITION_FILE"
