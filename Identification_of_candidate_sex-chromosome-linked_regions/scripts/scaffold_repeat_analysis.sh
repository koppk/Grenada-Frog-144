#!/bin/bash
set -euo pipefail

# scaffold_repeat_analysis.sh
#
# Layer 5: Repeat characterization across all 13 scaffolds.
#
# 1. Extract scaffold_1 through scaffold_13 from the scaffolded assembly
# 2. Run WindowMasker using shared k-mer counts (fast, ~minutes)
# 3. Run RepeatMasker with Anura Dfam library (slow, ~hours)
# 4. Parse results into per-Mb repeat density profiles
#
# Uses identical parameters to the contig-level runs in W-chr-workflow
# for direct comparison:
#   WindowMasker: shared counts from full primary assembly, -ustat, -outfmt interval
#   RepeatMasker: -species Anura -pa 6 -gff -xsmall --uncurated
#
# Requires in same directory (SCRIPTDIR):
#   parse_scaffold_repeats.py
#   mann_whitney_repeats.py
#
# Usage:
#   bash scaffold_repeat_analysis.sh <scaffolded_assembly.fasta> \
#       <shared_wm_counts.txt> <z_candidate_regions.tsv> <output_dir>
#
# Example:
#   cd /data/GrenadaFrog144/SexChromosomes/Z-chr-workflow
#   nohup bash scripts/scaffold_repeat_analysis.sh \
#       /data/GrenadaFrog144/Pristimantis_euphronides.genome.fasta \
#       /data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/windowmasker/wm_shared_counts.txt \
#       ./z_candidate_detection/z_candidate_regions.tsv \
#       ./z_candidate_detection/repeat_analysis \
#       > scaffold_repeat_analysis.log 2>&1 &
#
# Author: Kopp K., Pristimantis euphronides genome project

GENOME="${1:?Usage: $0 <scaffolded_assembly.fasta> <shared_wm_counts.txt> <z_candidate_regions.tsv> <output_dir>}"
WM_COUNTS="${2:?Usage: $0 <scaffolded_assembly.fasta> <shared_wm_counts.txt> <z_candidate_regions.tsv> <output_dir>}"
REGIONS="${3:?Usage: $0 <scaffolded_assembly.fasta> <shared_wm_counts.txt> <z_candidate_regions.tsv> <output_dir>}"
OUTDIR="${4:?Usage: $0 <scaffolded_assembly.fasta> <shared_wm_counts.txt> <z_candidate_regions.tsv> <output_dir>}"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"

# RepeatMasker parameters — identical to contig-level run
RM_THREADS=6
RM_SPECIES="Anura"
RM_BIN=$(which RepeatMasker 2>/dev/null || echo "/data/RepeatMasker/RepeatMasker")

mkdir -p "$OUTDIR"

echo "=== Scaffold repeat analysis (Layer 5) ==="
echo "Start: $(date)"
echo "Genome:     $GENOME"
echo "WM counts:  $WM_COUNTS"
echo "Z regions:  $REGIONS"
echo "Output:     $OUTDIR"
echo ""

# ── Step 1: Extract scaffold_1 through scaffold_13 ─────────────
SCAFFOLDS_FA="${OUTDIR}/scaffolds_1_to_13.fasta"

if [ -f "$SCAFFOLDS_FA" ]; then
    echo "[Step 1] Scaffolds FASTA exists, reusing: $SCAFFOLDS_FA"
else
    echo "[Step 1] Extracting scaffold_1 through scaffold_13..."

    # Ensure index exists
    [ ! -f "${GENOME}.fai" ] && samtools faidx "$GENOME"

    # Build region list from index
    SCAFFOLD_LIST=""
    for i in $(seq 1 13); do
        NAME=$(awk -v s="scaffold_${i}" '$1==s {print $1}' "${GENOME}.fai")
        if [ -z "$NAME" ]; then
            echo "  WARNING: scaffold_${i} not found in ${GENOME}.fai"
        else
            SCAFFOLD_LIST="${SCAFFOLD_LIST} ${NAME}"
        fi
    done

    samtools faidx "$GENOME" $SCAFFOLD_LIST > "$SCAFFOLDS_FA"
    samtools faidx "$SCAFFOLDS_FA"

    echo "  Extracted $(grep -c '^>' "$SCAFFOLDS_FA") scaffolds"
    echo "  Total: $(awk '{s+=$2} END {printf "%.1f Mb\n", s/1e6}' "${SCAFFOLDS_FA}.fai")"
fi
echo ""

# ── Step 2: WindowMasker ───────────────────────────────────────
WM_OUTDIR="${OUTDIR}/windowmasker"
WM_INTERVALS="${WM_OUTDIR}/wm_scaffolds_intervals.txt"
mkdir -p "$WM_OUTDIR"

if [ -f "$WM_INTERVALS" ]; then
    echo "[Step 2] WindowMasker intervals exist, reusing: $WM_INTERVALS"
else
    echo "[Step 2] Running WindowMasker with shared k-mer counts..."
    echo "  Shared counts: $WM_COUNTS"
    T0=$(date +%s)

    windowmasker -in "$SCAFFOLDS_FA" \
        -ustat "$WM_COUNTS" \
        -outfmt interval \
        -out "$WM_INTERVALS" \
        2> "${WM_OUTDIR}/wm_scaffolds.log"

    T1=$(date +%s)
    echo "  Done in $(( T1 - T0 )) seconds"
fi
echo ""

# ── Step 3: RepeatMasker ──────────────────────────────────────
RM_OUTDIR="${OUTDIR}/repeatmasker"
RM_DONE="${RM_OUTDIR}/scaffolds_1_to_13.fasta.tbl"
mkdir -p "$RM_OUTDIR"

if [ -f "$RM_DONE" ]; then
    echo "[Step 3] RepeatMasker output exists, reusing: $RM_OUTDIR"
else
    echo "[Step 3] Running RepeatMasker..."
    echo "  Species:  $RM_SPECIES"
    echo "  Threads:  $RM_THREADS"
    echo "  Options:  -gff -xsmall --uncurated"
    T0=$(date +%s)

    ${RM_BIN} \
        -species "${RM_SPECIES}" \
        -pa ${RM_THREADS} \
        -gff \
        -xsmall \
        --uncurated \
        -dir "${RM_OUTDIR}" \
        "${SCAFFOLDS_FA}"

    T1=$(date +%s)
    ELAPSED=$(( T1 - T0 ))
    HOURS=$(( ELAPSED / 3600 ))
    MINUTES=$(( (ELAPSED % 3600) / 60 ))
    echo "  Done in ${HOURS}h ${MINUTES}m"
    echo ""

    # Print summary
    if [ -f "$RM_DONE" ]; then
        echo "--- RepeatMasker summary ---"
        cat "$RM_DONE"
    fi
fi
echo ""

# ── Step 4: Parse into per-Mb profiles and compare ────────────
echo "[Step 4] Parsing repeat data into per-Mb profiles..."
echo ""

python3 "${SCRIPTDIR}/parse_scaffold_repeats.py" \
    "$SCAFFOLDS_FA" \
    "$WM_INTERVALS" \
    "${RM_OUTDIR}/scaffolds_1_to_13.fasta.out" \
    "$REGIONS" \
    "$OUTDIR"

echo ""

# ── Step 5: Statistical test ──────────────────────────────────
echo "[Step 5] Mann-Whitney U tests..."
echo ""
python3 "${SCRIPTDIR}/mann_whitney_repeats.py" "${OUTDIR}/scaffold_repeat_per_mb.tsv" "$OUTDIR"

echo "=== Scaffold repeat analysis complete ==="
echo "End: $(date)"
echo ""
echo "Output files:"
echo "  ${SCAFFOLDS_FA}  (extracted scaffolds)"
echo "  ${WM_INTERVALS}  (WindowMasker intervals)"
echo "  ${RM_OUTDIR}/  (RepeatMasker output)"
echo "  ${OUTDIR}/scaffold_repeat_per_mb.tsv  (per-Mb repeat density)"
echo "  ${OUTDIR}/scaffold_repeat_summary.tsv  (Z vs autosomal comparison)"
echo "  ${OUTDIR}/mann_whitney_repeats.tsv  (statistical test)"
