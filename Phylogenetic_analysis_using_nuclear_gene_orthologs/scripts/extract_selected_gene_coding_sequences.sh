#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# extract_selected_gene_coding_sequences.sh
# Query BESTMODELS-FINAL.bedDB column 26 for target gene names (with synonym
# mapping to Portik et al. 2023 nomenclature), extract CDS from
# BESTMODELS-FINAL.CDS.fa using SeqKit, rename headers to >PriEup_GENENAME.
#
# Input:  BESTMODELS-FINAL.bedDB, BESTMODELS-FINAL.CDS.fa
# Output: PriEup_gene_FASTA/PriEup_<GENE>.fasta

set -euo pipefail

OUTDIR="PriEup_gene_FASTA"
mkdir -p "$OUTDIR"

declare -A gene_map=(
    [NCX1]="SLC8A1"
    [CMYC]="MYC"
    [CRYBA]="CRYBA1"
    [MYH]="MYH6"
    [NT3]="NTF3"
    [SIA]="SIAH1"
)

genes=(
    CXCR4 RAG1 TYR NCX1 SLC8A3 BDNF BMP2 CMYC CRYBA MC1R
    MYH NT3 POMC RAG2 RHO SIA TNS3 ZFX
)

beddb="BESTMODELS-FINAL.bedDB"
cds_fa="BESTMODELS-FINAL.CDS.fa"

for gene in "${genes[@]}"; do
    mapped_name=${gene_map[$gene]:-}
    query=${mapped_name:-$gene}

    hannoID=$(awk -F'\t' -v g="$query" '$26 == g {print $4}' "$beddb")

    if [[ -z "$hannoID" ]]; then
        echo "$gene ($query): not found in $beddb"
        continue
    fi

    if [[ -n "$mapped_name" ]]; then
        filename="${OUTDIR}/PriEup_${mapped_name}_${gene}.fasta"
        header=">PriEup_${mapped_name}_${gene}"
    else
        filename="${OUTDIR}/PriEup_${gene}.fasta"
        header=">PriEup_${gene}"
    fi

    seqkit grep -r -p "$hannoID" "$cds_fa" > "$filename"
    sed -Ei "s/^>.*/$header/" "$filename"
    echo "$gene: extracted to $filename"
done
