#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# create_taxid_family_genus_species.sh
# Assign family, genus, and species from TaxIDs using taxonkit v0.19.0
# with NCBI Taxonomy database (taxdump).
#
# Input:  Acc_taxid.tsv (Accession, TaxID)
# Output: Acc_taxid_family_genus_species.tsv

set -euo pipefail

INPUT="Acc_taxid.tsv"
OUTPUT="Acc_taxid_family_genus_species.tsv"

tail -n+2 "$INPUT" \
    | cut -f2 \
    | taxonkit lineage -n -r \
    | taxonkit reformat -f '{f};{g};{s}' \
    | csvtk -H -t cut -f 1,5 \
    | csvtk -H -t sep -f 2 -s ';' -R \
    | csvtk -H -t add-header -n "taxid,family,genus,species" \
    > "$OUTPUT"

echo "Wrote $OUTPUT"
