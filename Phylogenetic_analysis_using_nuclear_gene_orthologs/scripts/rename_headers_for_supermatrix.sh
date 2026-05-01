#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# rename_headers_for_supermatrix.sh
# Standardise FASTA headers to Species_name__ACCESSION__|_Genus_|_Family format
# across all 12 retained gene alignments for supermatrix construction.
#
# Input:  refined/<GENE>_mafft.fasta
# Output: renamed/<GENE>_renamed.fasta

set -euo pipefail

INPUT_DIR="refined"
OUTPUT_DIR="renamed"

mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*_mafft.fasta; do
    gene=$(basename "$file" _mafft.fasta)
    outfile="$OUTPUT_DIR/${gene}_renamed.fasta"

    awk '
    BEGIN { OFS = "" }
    /^>/ {
        header = substr($0, 2)
        if (header ~ /PriEup/) {
            print ">Pristimantis_euphronides | PriEup"
            next
        }
        if (header ~ /lcl\|NC_053499.1/) {
            print ">Rana_temporaria | NC_053499.1"
            next
        }
        acc = gensub(/^([^\/ >]+).*/, "\\1", "g", header)
        species = gensub(/^.* ([A-Z][a-z]+ [a-z0-9._-]+).*/, "\\1", "g", header)
        gsub(" ", "_", species)
        print ">"species" | "acc
        next
    }
    { print }
    ' "$file" > "$outfile"

    echo "$gene: headers renamed"
done
