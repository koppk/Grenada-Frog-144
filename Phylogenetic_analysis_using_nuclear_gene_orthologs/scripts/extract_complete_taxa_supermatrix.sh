#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# extract_complete_taxa_supermatrix.sh
# Filter supermatrix to species with sequences available for all 12 genes
# (12-gene supermatrix taxon selection).
#
# Input:  supermatrix/supermatrix.fasta, complete_taxa.txt
# Output: supermatrix/supermatrix_complete_taxa.fasta

set -euo pipefail

SUPERMATRIX="supermatrix/supermatrix.fasta"
TAXON_LIST="complete_taxa.txt"
OUT="supermatrix/supermatrix_complete_taxa.fasta"

awk '
    NR == FNR { keep[$1] = 1; next }
    /^>/ {
        split($0, parts, /[| ]/)
        taxon = substr(parts[1], 2)
        printing = (taxon in keep)
    }
    printing
' "$TAXON_LIST" "$SUPERMATRIX" > "$OUT"

echo "Filtered to $(grep -c "^>" "$OUT") taxa: $OUT"
