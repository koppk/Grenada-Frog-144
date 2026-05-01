#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# concat_gene_alignments.sh
# Concatenate P. euphronides, Hyloidea, and R. temporaria sequences per gene
# into multi-species FASTA files for alignment.
#
# Input:  PriEup_gene_FASTA/, Hyloidea_sequences/, cds_genes_from_CDS_fasta_longest/
# Output: concatenated/<GENE>_PriEup_Hyloidea_RanTemp.fasta

set -euo pipefail

GENES="ZFX SLC8A1 MYC CRYBA1 MYH6 NTF3 SIAH1 CXCR4 RAG1 TYR SLC8A3 BDNF BMP2 MC1R POMC RAG2 RHO TNS3"
OUTDIR="concatenated"

mkdir -p "$OUTDIR"

for GENE in $GENES; do
    FILE1="./PriEup_gene_FASTA/PriEup_${GENE}.fasta"
    FILE2="./Hyloidea_sequences/${GENE}_Hyloidea.fasta"
    FILE3="./cds_genes_from_CDS_fasta_longest/${GENE}_Rana_temporaria_CDS.fasta"
    OUTFILE="${OUTDIR}/${GENE}_PriEup_Hyloidea_RanTemp.fasta"

    MISSING=""
    [[ ! -f "$FILE1" ]] && MISSING="$MISSING $FILE1"
    [[ ! -f "$FILE2" ]] && MISSING="$MISSING $FILE2"
    [[ ! -f "$FILE3" ]] && MISSING="$MISSING $FILE3"

    if [[ -n "$MISSING" ]]; then
        echo "$GENE: missing$MISSING"
        continue
    fi

    cat "$FILE1" "$FILE2" "$FILE3" > "$OUTFILE"
    echo "$GENE: concatenated"
done
