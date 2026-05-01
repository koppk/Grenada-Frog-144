#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Postprocess RagTag scaffold output:
#   1. Build renaming map (NC_ → scaffold_1–13, NW_ → scaffold_14–31)
#   2. Rename AGP
#   3. Rename and filter FASTA (remove contigs < 200 bp)
#   4. Extract contig length files for boxplot (Figure SR3)
#
# Input:  ragtag.scaffold.agp, ragtag.scaffold.fasta (RagTag output directory)
# Output: scaffold_renaming_map.tsv, ragtag.scaffold.renamed.agp,
#         Pristimantis_euphronides_genome.fasta, all_contigs.len, top13_contigs.len
#
# Usage:  bash postprocess_ragtag.sh <ragtag_output_dir> <output_dir>
#
# Requirements: awk, sed, seqkit

set -euo pipefail

RAGTAG_DIR="${1:?Usage: bash postprocess_ragtag.sh <ragtag_output_dir> <output_dir>}"
OUT_DIR="${2:?Usage: bash postprocess_ragtag.sh <ragtag_output_dir> <output_dir>}"

AGP="${RAGTAG_DIR}/ragtag.scaffold.agp"
FASTA="${RAGTAG_DIR}/ragtag.scaffold.fasta"

mkdir -p "${OUT_DIR}"

# --- Step 1: Build renaming map ---
# Extract unique scaffold names (NC_ and NW_), compute total span per scaffold,
# sort NC_ and NW_ separately by size (descending), assign scaffold numbers.

echo "Step 1: Building renaming map..."

# Get scaffold name and total span (max coordinate in column 3)
awk -F'\t' '$1 ~ /_RagTag$/ {
    if ($3 > max[$1]) max[$1] = $3
}
END {
    for (s in max) print s, max[s]
}' "${AGP}" | sort -k2,2nr > "${OUT_DIR}/scaffold_sizes_unsorted.tmp"

# NC_ scaffolds → scaffold_1 to scaffold_N (by E. coqui accession number)
# NW_ scaffolds → scaffold_(N+1) onwards (by accession number)
awk '
    /^NC_/ { nc[NR] = $0 }
    /^NW_/ { nw[NR] = $0 }
    END {
        for (i in nc) { split(nc[i], a, " "); printf "%s\t%s\n", a[1], a[2] }
    }
' "${OUT_DIR}/scaffold_sizes_unsorted.tmp" | sort > "${OUT_DIR}/scaffold_renaming_map_nc.tmp"

awk '
    /^NW_/ { nw[NR] = $0 }
    END {
        for (i in nw) { split(nw[i], a, " "); printf "%s\t%s\n", a[1], a[2] }
    }
' "${OUT_DIR}/scaffold_sizes_unsorted.tmp" | sort > "${OUT_DIR}/scaffold_renaming_map_nw.tmp"

NC_COUNT=$(wc -l < "${OUT_DIR}/scaffold_renaming_map_nc.tmp")

# Renumber NC_ scaffolds 1..N (E. coqui chromosome order)
awk -F'\t' '{printf "%s\tscaffold_%d\t%s\n", $1, NR, $2}' \
    "${OUT_DIR}/scaffold_renaming_map_nc.tmp" > "${OUT_DIR}/scaffold_renaming_map.tsv"

# Renumber NW_ scaffolds (N+1)..M (accession order)
awk -F'\t' -v start="${NC_COUNT}" '{printf "%s\tscaffold_%d\t%s\n", $1, NR+start, $2}' \
    "${OUT_DIR}/scaffold_renaming_map_nw.tmp" >> "${OUT_DIR}/scaffold_renaming_map.tsv"

rm "${OUT_DIR}/scaffold_sizes_unsorted.tmp" \
   "${OUT_DIR}/scaffold_renaming_map_nc.tmp" \
   "${OUT_DIR}/scaffold_renaming_map_nw.tmp"

TOTAL_SCAFFOLDS=$(wc -l < "${OUT_DIR}/scaffold_renaming_map.tsv")
echo "  ${NC_COUNT} NC_ scaffolds (scaffold_1–${NC_COUNT})"
echo "  $((TOTAL_SCAFFOLDS - NC_COUNT)) NW_ scaffolds (scaffold_$((NC_COUNT+1))–scaffold_${TOTAL_SCAFFOLDS})"

# --- Step 2: Rename AGP ---
echo "Step 2: Renaming AGP..."

# Build sed commands from renaming map
awk -F'\t' '{print "s/^" $1 "\\t/" $2 "\\t/"}' "${OUT_DIR}/scaffold_renaming_map.tsv" \
    > "${OUT_DIR}/rename_sed.tmp"

sed -f "${OUT_DIR}/rename_sed.tmp" "${AGP}" > "${OUT_DIR}/ragtag.scaffold.renamed.agp"
rm "${OUT_DIR}/rename_sed.tmp"

echo "  Written: ragtag.scaffold.renamed.agp"

# --- Step 3: Rename and filter FASTA ---
echo "Step 3: Renaming FASTA and removing contigs < 200 bp..."

# Build seqkit replace pattern from renaming map
awk -F'\t' '{gsub(/\./, "\\.", $1); print $1 "\t" $2}' "${OUT_DIR}/scaffold_renaming_map.tsv" \
    > "${OUT_DIR}/rename_seqkit.tmp"

seqkit replace -p '(.+)' -r '{kv}' -k "${OUT_DIR}/rename_seqkit.tmp" --keep-key "${FASTA}" \
    | seqkit seq -m 200 | gzip > "${OUT_DIR}/Pristimantis_euphronides_genome.fasta.gz"

rm "${OUT_DIR}/rename_seqkit.tmp"

echo "  Written: Pristimantis_euphronides_genome.fasta.gz"

# --- Step 4: Extract contig length files for boxplot ---
echo "Step 4: Extracting contig length files from renamed AGP..."

# all_contigs.len — all W-lines
awk -F'\t' '$5=="W" {print $6, $8}' "${OUT_DIR}/ragtag.scaffold.renamed.agp" \
    > "${OUT_DIR}/all_contigs.len"

# top13_contigs.len — contigs in scaffold_1 through scaffold_13
awk -F'\t' '$5=="W" && $1 ~ /^scaffold_(1[0-3]?|[1-9])$/ {print $6, $8}' \
    "${OUT_DIR}/ragtag.scaffold.renamed.agp" > "${OUT_DIR}/top13_contigs.len"

ALL_N=$(wc -l < "${OUT_DIR}/all_contigs.len")
TOP13_N=$(wc -l < "${OUT_DIR}/top13_contigs.len")
echo "  all_contigs.len: ${ALL_N} contigs"
echo "  top13_contigs.len: ${TOP13_N} contigs"

echo "Done."
