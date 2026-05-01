#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# get_taxids.sh
# Retrieve NCBI TaxIDs for GenBank accessions using NCBI Entrez efetch.
# P. euphronides accessions (PriEup_*) are assigned TaxID 448649 directly.
#
# Usage: get_taxids.sh <path/to/Acc_<GENE>.txt>
# Input:  Accession list (one per line)
# Output: <input_basename>.taxid.tsv (Accession, TaxID)

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <path/to/Acc_<GENE>.txt>"
    exit 1
fi

input="$1"

if [[ ! -f "$input" ]]; then
    echo "Error: file not found: $input"
    exit 1
fi

output="${input%.txt}.taxid.tsv"

echo -e "Accession\tTaxID" > "$output"

while read -r ACC; do
    ACC="${ACC//$'\r'/}"
    ACC="${ACC%%[[:space:]]}"

    echo -n -e "${ACC}\t" >> "$output"

    if [[ $ACC == PriEup_* ]]; then
        echo "448649" >> "$output"
    else
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=xml" \
            | grep -m1 TSeq_taxid \
            | cut -d '>' -f2 \
            | cut -d '<' -f1 \
            | tr -d "\n" \
            >> "$output"
        echo >> "$output"
    fi
done < "$input"

echo "Wrote $output"
