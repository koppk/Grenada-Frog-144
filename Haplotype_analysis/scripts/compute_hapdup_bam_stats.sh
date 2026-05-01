#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Compute read mapping, haplotagging, and structural variant statistics
# from HapDup intermediate files.
# Reproduces numbers in Additional file 5, Section
# "Haplotype Phasing and Dual Assembly".
#
# Input:  <hapdup_run_dir> containing:
#           lr_mapping.bam                        (read mapping to assembly)
#           hapdup/margin/MARGIN_PHASED.haplotagged.bam  (haplotagged reads)
#           hapdup/structural/inversions.bed       (detected inversions)
#           hapdup/structural/breakpoints_all.csv  (all breakpoints)
#           hapdup/structural/breakpoints_balanced.csv
#
# Output: stdout (redirect to file for documentation)
#
# Usage: ./compute_hapdup_bam_stats.sh <hapdup_run_dir> [threads]
#        Default: 16 threads
#
# Note: lr_mapping.bam and MARGIN_PHASED.haplotagged.bam are not deposited
# on Zenodo due to size. They can be regenerated using run_minimap2_HapDup_prep.sh
# and run_docker_HapDup.sh from the deposited reads and assembly.

DIR=${1:-.}

LR_BAM="${DIR}/lr_mapping.bam"
HTAG_BAM="${DIR}/hapdup/margin/MARGIN_PHASED.haplotagged.bam"
INV_BED="${DIR}/hapdup/structural/inversions.bed"
BP_ALL="${DIR}/hapdup/structural/breakpoints_all.csv"
BP_BAL="${DIR}/hapdup/structural/breakpoints_balanced.csv"

THREADS=${2:-16}

# Check inputs
for f in "$LR_BAM" "$HTAG_BAM" "$INV_BED" "$BP_ALL" "$BP_BAL"; do
    if [ ! -f "$f" ]; then
        echo "Error: file not found: $f" >&2
        exit 1
    fi
done

echo "=== Read mapping statistics (lr_mapping.bam) ==="
samtools flagstat -@ ${THREADS} "$LR_BAM"
echo ""
echo "Average read length:"
samtools stats -@ ${THREADS} "$LR_BAM" | grep "^SN" | grep "average length"

echo ""
echo "=== Haplotagging statistics (MARGIN_PHASED.haplotagged.bam) ==="
TOTAL=$(samtools view -@ ${THREADS} -c "$HTAG_BAM")
samtools view -@ ${THREADS} "$HTAG_BAM" | awk -v total="$TOTAL" '{
    hp_found = 0
    for (i = 12; i <= NF; i++) {
        if ($i == "HP:i:1") { h1++; hp_found = 1; break }
        if ($i == "HP:i:2") { h2++; hp_found = 1; break }
    }
    if (!hp_found) h0++
} END {
    printf "Total reads in haplotagged BAM: %d\n", total
    printf "Haplotype 1 (HP:i:1): %d\n", h1
    printf "Haplotype 2 (HP:i:2): %d\n", h2
    printf "Unphased (no HP tag): %d\n", h0
    printf "Total haplotagged (HP:i:1 + HP:i:2): %d\n", h1 + h2
}'

echo ""
echo "=== Structural variant detection (HapDup step 5) ==="
INV_COUNT=$(awk 'NF>0 && !/^#/' "$INV_BED" | wc -l)
BP_ALL_COUNT=$(awk 'NF>0 && !/^#/' "$BP_ALL" | wc -l)
BP_BAL_COUNT=$(awk 'NF>0 && !/^#/' "$BP_BAL" | wc -l)
echo "Inversions detected: ${INV_COUNT}"
echo "Breakpoints (all): ${BP_ALL_COUNT}"
echo "Breakpoints (balanced): ${BP_BAL_COUNT}"
