#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Compute haplotype alignment statistics from D-Genies output files.
# Reproduces Table HR2 in Additional file 5.
#
# Input:  <comparison_dir> containing:
#           map_*.paf.gz           (D-Genies PAF alignment)
#           *_assoc.tsv            (D-Genies association table)
# Output: stdout
#
# Metrics reported:
#   - Haplotype 2 bases aligned (bp and % of total), from assoc.tsv
#   - Number of alignment blocks, from PAF
#   - Total matching base pairs, from PAF column 10
#   - Mean alignment block length (bp), from PAF column 10
#   - Weighted mean per-base divergence (%), from PAF dv tag

DIR=${1:-.}

PAF=$(ls "${DIR}"/map_*.paf.gz 2>/dev/null | head -1)
ASSOC=$(ls "${DIR}"/*_assoc.tsv 2>/dev/null | head -1)

if [ -z "$PAF" ] || [ -z "$ASSOC" ]; then
    echo "Error: PAF.gz and assoc.tsv not found in ${DIR}" >&2
    exit 1
fi

echo "=== Input files ==="
echo "PAF: $(basename $PAF)"
echo "Assoc: $(basename $ASSOC)"
echo ""

echo "=== Haplotype 2 (query) alignment span from assoc.tsv ==="
awk -F'\t' '
    NR > 1 && $2 != "None" { span += $6 - $5; total += $4 }
    NR > 1 && $2 == "None" { total += $4 }
    END {
        printf "Aligned span: %d bp\n", span
        printf "Total haplotype 2 length: %d bp\n", total
        printf "Fraction aligned: %.1f%%\n", span / total * 100
    }
' "$ASSOC"

echo ""
echo "=== Alignment block statistics from PAF ==="
zcat "$PAF" | awk -F'\t' '{
    blocks++
    match_bp += $10
    for (i = 13; i <= NF; i++) {
        if ($i ~ /^dv:/) {
            split($i, a, ":")
            dv_sum += a[3] * $11
            dv_total += $11
        }
    }
}
END {
    printf "Number of alignment blocks: %d\n", blocks
    printf "Total matching base pairs: %d bp\n", match_bp
    printf "Mean alignment block length: %.0f bp\n", match_bp / blocks
    printf "Weighted mean per-base divergence: %.2f%%\n", (dv_sum / dv_total) * 100
}'
