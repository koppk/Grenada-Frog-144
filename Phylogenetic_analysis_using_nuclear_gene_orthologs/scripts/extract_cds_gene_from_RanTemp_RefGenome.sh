#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# extract_cds_gene_from_RanTemp_RefGenome.sh
# Extract longest CDS isoform per gene from the Rana temporaria reference
# genome (GCF_905171775.1) using SeqKit. R. temporaria serves as the
# non-Hyloidea outgroup.
#
# Input:  GCF_905171775.1_aRanTem1.1_cds_from_genomic.fna.gz
#         gene_synonyms_ForRanTemp.txt
# Output: cds_genes_from_CDS_fasta_longest/<GENE>_Rana_temporaria_CDS.fasta

set -euo pipefail

CDS_FASTA="GCF_905171775.1_aRanTem1.1_cds_from_genomic.fna.gz"
GENE_LIST="gene_synonyms_ForRanTemp.txt"
OUTDIR="cds_genes_from_CDS_fasta_longest"

mkdir -p "$OUTDIR"

while IFS= read -r line || [[ -n "$line" ]]; do
    IFS=',' read -ra SYNONYMS <<< "$line"
    GENE="${SYNONYMS[0]}"

    REGEX=""
    for synonym in "${SYNONYMS[@]}"; do
        [[ -z "$REGEX" ]] && REGEX="\\[gene=${synonym}\\]" || REGEX="${REGEX}|\\[gene=${synonym}\\]"
    done

    TMP_MATCHES=$(mktemp)
    seqkit grep -r -n -p "$REGEX" "$CDS_FASTA" > "$TMP_MATCHES"

    NUM_MATCHES=$(grep -c "^>" "$TMP_MATCHES" || true)
    if [[ "$NUM_MATCHES" -eq 0 ]]; then
        echo "$GENE: no CDS found"
        rm -f "$TMP_MATCHES"
        continue
    fi

    seqkit sort -l -r "$TMP_MATCHES" | seqkit head -n 1 > "$OUTDIR/${GENE}_Rana_temporaria_CDS.fasta"
    echo "$GENE: longest CDS saved ($NUM_MATCHES isoforms)"
    rm -f "$TMP_MATCHES"

done < "$GENE_LIST"
