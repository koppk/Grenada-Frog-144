#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# make_table.sh
# Fetch TaxID, scientific name, genus, and family for GenBank accessions
# using NCBI Entrez utilities. Alternative to get_taxids.sh + taxonkit
# workflow for cases where taxonkit is unavailable.
#
# Input:  accession.txt (one accession per line)
# Output: stdout (TSV: Accession, TaxID, ScientificName, Genus, Family)

set -euo pipefail

INPUT="accession.txt"

echo -e "Accession\tTaxID\tScientificName\tGenus\tFamily"

while read -r acc; do
    echo "Processing $acc" >&2

    TAXID=$(esearch -db nucleotide -query "$acc" | \
            elink -target taxonomy | \
            efetch -format xml | \
            xtract -pattern Taxon -element TaxId)

    if [[ -z "$TAXID" ]]; then
        echo -e "$acc\tNA\tNA\tNA\tNA"
        continue
    fi

    XML=$(esearch -db taxonomy -query "$TAXID" | efetch -format xml)

    SCI_NAME=$(echo "$XML" | xtract -pattern Taxon -element ScientificName)

    GENUS=$(echo "$XML" | \
        xtract -pattern Taxon -block LineageEx -if Rank -equals genus -element ScientificName)

    FAMILY=$(echo "$XML" | \
        xtract -pattern Taxon -block LineageEx -if Rank -equals family -element ScientificName)

    echo -e "$acc\t$TAXID\t$SCI_NAME\t$GENUS\t$FAMILY"
done < "$INPUT"
