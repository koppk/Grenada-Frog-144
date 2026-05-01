#!/usr/bin/env bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Extract contig length and placement information from the RagTag AGP file
# and split contigs into three groups based on their placement:
#   (1) contigs placed in scaffolds 1-13
#   (2) contigs placed in scaffolds 14-31
#   (3) unplaced contigs
#
# Contigs <200 bp are excluded to match the deposited genome FASTA.
#
# Input:  Workflow_1/output/ragtag.scaffold.renamed.agp
# Output: Workflow_1/output/scaffolds_1-13.len
#         Workflow_1/output/scaffolds_14-31.len
#         Workflow_1/output/unplaced_contigs.len
#
# Each .len file is two-column: <contig_name> <length_bp>
#
# Usage:  bash split_contig_length_and_placement.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$(cd "${SCRIPT_DIR}/../output" && pwd)"
AGP="${OUT_DIR}/ragtag.scaffold.renamed.agp"

if [[ ! -f "${AGP}" ]]; then
    echo "ERROR: AGP file not found: ${AGP}" >&2
    exit 1
fi

echo "Reading AGP: ${AGP}" >&2
echo "Writing to:  ${OUT_DIR}" >&2

# Group 1: scaffolds 1-13
awk -F'\t' '$5=="W" && $1 ~ /^scaffold_([1-9]|1[0-3])$/ && ($8 - $7 + 1) >= 200 \
            {print $6, $8 - $7 + 1}' \
    "${AGP}" > "${OUT_DIR}/scaffolds_1-13.len"

# Group 2: scaffolds 14-31
awk -F'\t' '$5=="W" && $1 ~ /^scaffold_(1[4-9]|2[0-9]|3[01])$/ && ($8 - $7 + 1) >= 200 \
            {print $6, $8 - $7 + 1}' \
    "${AGP}" > "${OUT_DIR}/scaffolds_14-31.len"

# Group 3: unplaced
awk -F'\t' '$5=="W" && $1 !~ /^scaffold_([1-9]|1[0-9]|2[0-9]|3[01])$/ && ($8 - $7 + 1) >= 200 \
            {print $6, $8 - $7 + 1}' \
    "${AGP}" > "${OUT_DIR}/unplaced_contigs.len"

# Sanity check
n1=$(wc -l < "${OUT_DIR}/scaffolds_1-13.len")
n2=$(wc -l < "${OUT_DIR}/scaffolds_14-31.len")
n3=$(wc -l < "${OUT_DIR}/unplaced_contigs.len")
total=$(( n1 + n2 + n3 ))

echo >&2
echo "scaffolds_1-13.len      ${n1} contigs" >&2
echo "scaffolds_14-31.len     ${n2} contigs" >&2
echo "unplaced_contigs.len    ${n3} contigs" >&2
echo "                        ----- " >&2
echo "total                   ${total} contigs" >&2
