#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Compute HapDup phasing statistics from PEPPER and Margin output.
# Reproduces numbers in Additional file 5, Section
# "Haplotype Phasing and Dual Assembly".
#
# Input:  hapdup_intermediates/ containing:
#           PEPPER_VARIANT_FULL.vcf.gz
#           MARGIN_PHASED.phased.vcf.gz
#           MARGIN_PHASED.phaseset.bed
# Output: stdout

DIR=${1:-.}

echo "=== Variant counts ==="
TOTAL=$(zcat "${DIR}/PEPPER_VARIANT_FULL.vcf.gz" | grep -c -v "^#")
PHASED=$(zcat "${DIR}/MARGIN_PHASED.phased.vcf.gz" | grep -v "^#" | grep -c "|")
echo "Total variants (PEPPER): ${TOTAL}"
echo "Phased variants (genotype separator '|'): ${PHASED}"
awk "BEGIN {printf \"Phased fraction: %.1f%%\n\", ${PHASED}/${TOTAL}*100}"

echo ""
echo "=== Phase block statistics ==="
awk -F'\t' '{
    len = $3 - $2
    sum += len
    n++
    lens[n] = len
}
END {
    printf "Number of phase blocks: %d\n", n
    printf "Total phased span (bp): %d\n", sum
    printf "Average block length (bp): %.0f\n", sum / n

    asort(lens)
    cumsum = 0
    for (i = 1; i <= n; i++) {
        cumsum += lens[i]
        if (cumsum >= sum / 2) {
            printf "Block N50 (bp): %d\n", lens[i]
            break
        }
    }
}' "${DIR}/MARGIN_PHASED.phaseset.bed"
