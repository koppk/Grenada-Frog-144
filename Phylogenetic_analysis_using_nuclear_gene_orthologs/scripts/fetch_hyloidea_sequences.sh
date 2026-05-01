#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# fetch_hyloidea_sequences.sh
# Download Hyloidea homologs from NCBI GenBank for each target gene using
# Entrez queries with gene-specific synonyms.
#
# Input:  gene_synonyms.txt (gene names with comma-separated synonyms)
# Output: Hyloidea_sequences/<GENE>_Hyloidea.fasta
#         Hyloidea_sequences/<GENE>_uids.txt

set -euo pipefail

OUTDIR="Hyloidea_sequences"
LOGFILE="gene_fetch_log.txt"

mkdir -p "$OUTDIR"
> "$LOGFILE"

FILTERS='AND (cds[Feature key] OR gene[Title] OR mRNA[Title]) NOT genomic[Title] NOT chromosome[Title] NOT "complete genome"[Title]'

exec 3< gene_synonyms.txt

while IFS= read -r line <&3 || [[ -n "$line" ]]; do
    IFS=',' read -ra SYNONYMS <<< "$line"
    GENE_NAME="${SYNONYMS[0]}"

    GENE_QUERY=""
    for synonym in "${SYNONYMS[@]}"; do
        [[ -z "$GENE_QUERY" ]] && GENE_QUERY="${synonym}[Gene]" || GENE_QUERY="${GENE_QUERY} OR ${synonym}[Gene]"
    done

    FULL_QUERY="(${GENE_QUERY}) AND Hyloidea[Organism] $FILTERS"
    echo "Query: $FULL_QUERY" >> "$LOGFILE"

    UID_FILE="${OUTDIR}/${GENE_NAME}_uids.txt"
    esearch -db nucleotide -query "$FULL_QUERY" | efetch -format uid > "$UID_FILE"

    NUM_UIDS=$(grep -c . "$UID_FILE" || true)
    if [[ "$NUM_UIDS" -eq 0 ]]; then
        echo "$GENE_NAME: no UIDs found, skipping" | tee -a "$LOGFILE"
        continue
    fi

    echo "$GENE_NAME: $NUM_UIDS UIDs" | tee -a "$LOGFILE"

    FASTA_OUT="${OUTDIR}/${GENE_NAME}_Hyloidea.fasta"
    > "$FASTA_OUT"

    while IFS= read -r uid || [[ -n "$uid" ]]; do
        TMP_FASTA=$(mktemp)
        if efetch -db nucleotide -id "$uid" -format fasta > "$TMP_FASTA" 2>> "$LOGFILE"; then
            if grep -q "^>" "$TMP_FASTA"; then
                cat "$TMP_FASTA" >> "$FASTA_OUT"
            else
                echo "$GENE_NAME: UID $uid returned no FASTA header" >> "$LOGFILE"
            fi
        else
            echo "$GENE_NAME: efetch failed for UID $uid" >> "$LOGFILE"
        fi
        rm -f "$TMP_FASTA"
        sleep 0.3
    done < "$UID_FILE"

    SEQ_COUNT=$(grep -c "^>" "$FASTA_OUT" || true)
    echo "$GENE_NAME: $SEQ_COUNT sequences saved" | tee -a "$LOGFILE"

done

exec 3<&-
