#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# run_windowmasker_two_sets.sh
# ============================
# Run WindowMasker on placed and unplaced contig sets using a SHARED
# k-mer frequency table built from the full unscaffolded primary assembly.
#
# WHY shared counts?
# WindowMasker's -mk_counts step calibrates what's "repetitive" based
# on k-mer frequencies in the input. If run independently on each set,
# a repeat-rich set would recalibrate its threshold so repeats become
# the "normal" — masking LESS, not more. By building counts from the
# FULL unscaffolded primary assembly (all P. euphronides contigs before
# any RagTag scaffolding), "repetitive" is defined relative to the
# complete genome.
#
# Workflow:
#   1. Build k-mer counts from full primary assembly (once)
#   2. Apply those counts to placed and unplaced sets with -ustat
#   3. Call parse_windowmasker_results.py to compute per-contig stats
#
# Usage:
#   bash run_windowmasker_two_sets.sh
#
# Date: 2026-02-16

set -euo pipefail

# === Paths — adjust as needed ===
BASEDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"
OUTDIR="${BASEDIR}/windowmasker"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"

PLACED="${BASEDIR}/placed_contigs.fasta"
UNPLACED="${BASEDIR}/unplaced_contigs.fasta"

# Reference for k-mer counts: full unscaffolded primary assembly
FULL_ASM="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"

# === Check dependencies ===
echo "=== WindowMasker two-set comparison (shared counts) ==="
echo "Start: $(date)"
echo ""

if ! command -v windowmasker &> /dev/null; then
    echo "ERROR: windowmasker not found on PATH"
    echo "Install: sudo apt install ncbi-blast+"
    exit 1
fi

for f in "$FULL_ASM" "$PLACED" "$UNPLACED"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

mkdir -p "$OUTDIR"

# ============================================================
# Step 1: Build shared k-mer counts from FULL primary assembly
# ============================================================
SHARED_COUNTS="${OUTDIR}/wm_shared_counts.txt"

if [ -f "$SHARED_COUNTS" ]; then
    echo "Shared counts file exists, reusing: $SHARED_COUNTS"
else
    echo "[Step 1] Building shared k-mer counts from full primary assembly ..."
    echo "  Input: $FULL_ASM"
    T0=$(date +%s)
    windowmasker -in "$FULL_ASM" \
        -mk_counts \
        -out "$SHARED_COUNTS" \
        2> "${OUTDIR}/wm_mk_counts.log"
    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
fi
echo ""

# ============================================================
# Step 2: Apply shared counts to each set
# ============================================================
echo "[Step 2] Applying shared counts to placed and unplaced sets ..."
echo ""

for LABEL in placed unplaced; do
    case $LABEL in
        placed)   FASTA="$PLACED" ;;
        unplaced) FASTA="$UNPLACED" ;;
    esac

    INTERVALS="${OUTDIR}/wm_${LABEL}_intervals.txt"

    if [ -f "$INTERVALS" ]; then
        echo "  [$LABEL] Intervals file exists, reusing: $INTERVALS"
        continue
    fi

    echo "  [$LABEL] Running windowmasker -ustat ..."
    T0=$(date +%s)
    windowmasker -in "$FASTA" \
        -ustat "$SHARED_COUNTS" \
        -outfmt interval \
        -out "$INTERVALS" \
        2> "${OUTDIR}/wm_${LABEL}_ustat.log"
    T1=$(date +%s)
    echo "    Done in $((T1 - T0)) seconds → $INTERVALS"
done
echo ""

# ============================================================
# Step 3: Parse intervals → per-contig masked fraction tables
# ============================================================
echo "[Step 3] Parsing WindowMasker output ..."
echo ""

PARSE_SCRIPT="${SCRIPTDIR}/parse_windowmasker_results.py"
if [ ! -f "$PARSE_SCRIPT" ]; then
    echo "ERROR: $PARSE_SCRIPT not found"
    echo "       Place it in the same directory as this script."
    exit 1
fi

python3 "$PARSE_SCRIPT" \
    --outdir "$OUTDIR" \
    --placed "$PLACED" \
    --unplaced "$UNPLACED"

echo ""
echo "Finished: $(date)"
echo "All output in: $OUTDIR"
echo "=== Done ==="
