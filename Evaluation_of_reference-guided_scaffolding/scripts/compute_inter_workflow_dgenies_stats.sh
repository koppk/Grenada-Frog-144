#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Compute inter-workflow concordance statistics from D-Genies v1.5.0
# output files.
#
# Reproduces the inter-workflow alignment statistics reported in
# Additional file 3:
#   - "Assembly quality assessment" (pre-scaffolding comparison)
#   - "Evaluation of reference-guided scaffolding" (post-scaffolding)
#
# Input:  <comparison_dir> containing:
#           map_*.paf or map_*.paf.gz  (D-Genies PAF alignment)
#           *_assoc.tsv               (D-Genies association table)
#           *.tar.gz                  (D-Genies backup archive with
#                                      query.idx and target.idx)
#           no_query_matches_*        (unmatched query list)
#           no_target_matches_*       (unmatched target list)
#
# Output: summary to stdout and <comparison_dir>/stats_summary.txt
#
# Usage:
#   bash compute_inter_workflow_dgenies_stats.sh <comparison_dir>

set -euo pipefail

DIR="${1:-.}"
OUTFILE="${DIR}/stats_summary.txt"
exec > >(tee "$OUTFILE") 2>&1

# --- Locate files ---

PAF_GZ=$(ls "${DIR}"/map_*.paf.gz 2>/dev/null | head -1 || true)
PAF_PLAIN=$(ls "${DIR}"/map_*.paf 2>/dev/null | head -1 || true)
ASSOC=$(ls "${DIR}"/*_assoc.tsv 2>/dev/null | head -1 || true)
BACKUP=$(ls "${DIR}"/*.tar.gz 2>/dev/null | head -1 || true)
NO_QUERY=$(ls "${DIR}"/no_query_matches_* 2>/dev/null | head -1 || true)
NO_TARGET=$(ls "${DIR}"/no_target_matches_* 2>/dev/null | head -1 || true)

if [ -n "$PAF_GZ" ]; then
    PAF="$PAF_GZ"
    CAT="zcat"
elif [ -n "$PAF_PLAIN" ]; then
    PAF="$PAF_PLAIN"
    CAT="cat"
else
    echo "Error: no map_*.paf or map_*.paf.gz found in ${DIR}" >&2
    exit 1
fi

[ -z "$ASSOC" ] && { echo "Error: no *_assoc.tsv found in ${DIR}" >&2; exit 1; }
[ -z "$BACKUP" ] && { echo "Error: no *.tar.gz backup found in ${DIR}" >&2; exit 1; }

# --- Extract index files from backup archive ---

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

tar xzf "$BACKUP" -C "$TMPDIR" query.idx target.idx 2>/dev/null

QUERY_IDX=$(find "$TMPDIR" -name query.idx | head -1)
TARGET_IDX=$(find "$TMPDIR" -name target.idx | head -1)

[ -z "$QUERY_IDX" ] && { echo "Error: query.idx not found in backup archive" >&2; exit 1; }
[ -z "$TARGET_IDX" ] && { echo "Error: target.idx not found in backup archive" >&2; exit 1; }

# --- Input summary ---

echo "=== Input files ==="
echo "PAF:        $(basename "$PAF")"
echo "Assoc:      $(basename "$ASSOC")"
echo "Backup:     $(basename "$BACKUP")"
[ -n "$NO_QUERY" ] && echo "No-query:   $(basename "$NO_QUERY")"
[ -n "$NO_TARGET" ] && echo "No-target:  $(basename "$NO_TARGET")"
echo ""

# --- Assembly sizes from index files ---

echo "=== Assembly sizes from index files ==="
awk 'NR==1 {name=$0} NR>1 {count++; sum+=$2} END {
    printf "Query:  %s — %d sequences, %d bp (%.1f Mb)\n", name, count, sum, sum/1e6
}' "$QUERY_IDX"

awk 'NR==1 {name=$0} NR>1 {count++; sum+=$2} END {
    printf "Target: %s — %d sequences, %d bp (%.1f Mb)\n", name, count, sum, sum/1e6
}' "$TARGET_IDX"
echo ""

# --- Query statistics from association table ---

echo "=== Query statistics from association table ==="
awk -F'\t' '
    NR > 1 && $2 != "None" { matched++; matched_bp += $4 }
    NR > 1 && $2 == "None" { unmatched++; unmatched_bp += $4 }
    END {
        total = matched + unmatched
        total_bp = matched_bp + unmatched_bp
        printf "Total query sequences: %d (%d bp, %.1f Mb)\n", total, total_bp, total_bp/1e6
        printf "Matched:   %d of %d (%.1f%%)\n", matched, total, matched/total*100
        printf "Matched bp: %d of %d (%.1f Mb of %.1f Mb)\n", matched_bp, total_bp, matched_bp/1e6, total_bp/1e6
        if (total_bp > 0) printf "Fraction:  %.2f%%\n", matched_bp/total_bp*100
        printf "Unmatched: %d (%d bp, %.3f Mb)\n", unmatched, unmatched_bp, unmatched_bp/1e6
    }
' "$ASSOC"
echo ""

# --- Target statistics from PAF + index ---

echo "=== Target statistics from PAF (unique target names with any alignment) ==="

# Matched targets: unique target names in PAF column 6
# Target lengths: PAF column 7 (target sequence length)
$CAT "$PAF" | awk '{print $6, $7}' | sort -k1,1 -u | awk '{
    count++; sum += $2
} END {
    printf "Matched target sequences: %d\n", count
    printf "Matched target bp: %d (%.1f Mb)\n", sum, sum/1e6
}' 

# Total targets from index
TOTAL_TARGET_COUNT=$(awk 'NR>1' "$TARGET_IDX" | wc -l)
TOTAL_TARGET_BP=$(awk 'NR>1 {sum+=$2} END {print sum}' "$TARGET_IDX")

# Matched target count and bp from PAF
MATCHED_TARGET=$($CAT "$PAF" | awk '{print $6, $7}' | sort -k1,1 -u | awk '{count++; sum+=$2} END {printf "%d\t%d", count, sum}')
MATCHED_TARGET_COUNT=$(echo "$MATCHED_TARGET" | cut -f1)
MATCHED_TARGET_BP=$(echo "$MATCHED_TARGET" | cut -f2)

UNMATCHED_TARGET_COUNT=$((TOTAL_TARGET_COUNT - MATCHED_TARGET_COUNT))
UNMATCHED_TARGET_BP=$((TOTAL_TARGET_BP - MATCHED_TARGET_BP))

if [ "$UNMATCHED_TARGET_COUNT" -gt 0 ]; then
    UNMATCHED_TARGET_MEAN=$((UNMATCHED_TARGET_BP / UNMATCHED_TARGET_COUNT))
else
    UNMATCHED_TARGET_MEAN=0
fi

PCT_TARGET=$(awk "BEGIN {printf \"%.1f\", ${MATCHED_TARGET_BP}/${TOTAL_TARGET_BP}*100}")

printf "Total target sequences: %d (%d bp, %.1f Mb)\n" "$TOTAL_TARGET_COUNT" "$TOTAL_TARGET_BP" "$(awk "BEGIN {printf \"%.1f\", ${TOTAL_TARGET_BP}/1e6}")"
printf "Matched:   %d of %d\n" "$MATCHED_TARGET_COUNT" "$TOTAL_TARGET_COUNT"
printf "Matched bp: %d of %d (%s%%)\n" "$MATCHED_TARGET_BP" "$TOTAL_TARGET_BP" "$PCT_TARGET"
printf "           %.1f of %.1f Mb\n" "$(awk "BEGIN {printf \"%.1f\", ${MATCHED_TARGET_BP}/1e6}")" "$(awk "BEGIN {printf \"%.1f\", ${TOTAL_TARGET_BP}/1e6}")"
printf "Unmatched: %d (%d bp, %.1f Mb)\n" "$UNMATCHED_TARGET_COUNT" "$UNMATCHED_TARGET_BP" "$(awk "BEGIN {printf \"%.1f\", ${UNMATCHED_TARGET_BP}/1e6}")"
printf "Unmatched mean: %d bp (%.1f kb)\n" "$UNMATCHED_TARGET_MEAN" "$(awk "BEGIN {printf \"%.1f\", ${UNMATCHED_TARGET_MEAN}/1e3}")"
echo ""

# --- Cross-check with no-match files ---

echo "=== Cross-check with D-Genies no-match files ==="
if [ -n "$NO_QUERY" ]; then
    NQ=$(grep -c . "$NO_QUERY" 2>/dev/null || echo 0)
    printf "no_query_matches file:  %d entries\n" "$NQ"
fi
if [ -n "$NO_TARGET" ]; then
    NT=$(grep -c . "$NO_TARGET" 2>/dev/null || echo 0)
    printf "no_target_matches file: %d entries (PAF unmatched: %d)\n" "$NT" "$UNMATCHED_TARGET_COUNT"
fi
echo ""

# --- Alignment block statistics from PAF ---

echo "=== Alignment block statistics from PAF ==="
$CAT "$PAF" | awk -F'\t' '{
    blocks++
    match_bp += $10
    block_len += $11
    for (i = 13; i <= NF; i++) {
        if ($i ~ /^dv:/) {
            split($i, a, ":")
            dv_sum += a[3] * $11
            dv_total += $11
        }
    }
}
END {
    printf "Alignment blocks: %d\n", blocks
    printf "Total matching bp: %d\n", match_bp
    printf "Total block length: %d\n", block_len
    if (blocks > 0) printf "Mean block length (matching bp): %.0f bp\n", match_bp / blocks
    if (block_len > 0) printf "Overall identity (matching/block): %.2f%%\n", match_bp / block_len * 100
    if (dv_total > 0) printf "Weighted mean per-base divergence: %.2f%%\n", (dv_sum / dv_total) * 100
}'
