#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Compute assembly statistics for Table HR1 in Additional file 5.
# Runs seqkit stats on FASTA files and extracts gap/unplaced
# statistics from RagTag AGP files.
#
# Usage:
#   ./compute_assembly_stats.sh <fasta1> [fasta2 ...]
#   ./compute_assembly_stats.sh --agp <agp_file> <total_sequences>
#
# For Table HR1 pre-scaffolding columns:
#   ./compute_assembly_stats.sh hapdup_dual_1.fasta.gz hapdup_dual_2.fasta.gz
#
# For Table HR1 post-scaffolding columns:
#   ./compute_assembly_stats.sh hapdup_dual_1.ragtag.scaffold.fasta.gz hapdup_dual_2.ragtag.scaffold.fasta.gz
#   ./compute_assembly_stats.sh --agp ragtag_scaffold_dual_1/ragtag.scaffold.agp 9852
#   ./compute_assembly_stats.sh --agp ragtag_scaffold_dual_2/ragtag.scaffold.agp 9753

if [ "$1" = "--agp" ]; then
    AGP=$2
    TOTAL_SEQ=$3
    echo "=== AGP statistics: $(basename $AGP) ==="

    UNPLACED=$(grep -c "^contig" "$AGP")
    UNPLACED_BP=$(awk '$5=="W" && $1~/^contig/ {sum+=$3-$2+1} END {print sum}' "$AGP")
    GAPS=$(awk '$5=="N" || $5=="U" {n++} END {print n}' "$AGP")
    GAP_BP=$(awk '($5=="N" || $5=="U") {sum+=$3-$2+1} END {print sum}' "$AGP")

    printf "Unplaced sequences: %d\n" "$UNPLACED"
    printf "Unplaced bases (bp): %d\n" "$UNPLACED_BP"
    printf "Gap sequences: %d\n" "$GAPS"
    printf "Gap bases (bp): %d\n" "$GAP_BP"
else
    echo "=== Assembly statistics (seqkit stats -a) ==="
    for f in "$@"; do
        echo "--- $(basename $f) ---"
        seqkit stats -a "$f"
        echo ""
    done
fi
