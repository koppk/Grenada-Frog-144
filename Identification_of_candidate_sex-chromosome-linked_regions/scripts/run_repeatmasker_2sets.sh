#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
# =============================================================================
# run_repeatmasker_2sets.sh
#
# Runs RepeatMasker on placed and unplaced contig sets to classify repeat
# families (SINEs, LINEs, DNA transposons, LTR elements, satellites, etc.)
# for downstream compositional banding and contig characterisation.
#
# Requires:
#   - RepeatMasker installed and on PATH (or adjust RM_BIN below)
#   - Dfam partitions 0 + 12 installed in RepeatMasker/Libraries/famdb/
#   - RMBlast on PATH
#
# Usage:
#   nohup bash run_repeatmasker_2sets.sh > repeatmasker_run.log 2>&1 &
#
# Output per set in OUTDIR:
#   <n>.out          - main annotation table
#   <n>.tbl          - summary statistics (% masked, TE breakdown)
#   <n>.masked       - masked FASTA (Ns replacing repeats)
#   <n>.out.gff      - GFF annotation
# =============================================================================

set -euo pipefail

# ----- CONFIGURATION -----
WORKDIR=/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison
OUTDIR=${WORKDIR}/RepeatMasker_output
THREADS=6           # 6 x 4 cores = 24 cores total (matches nproc 24)
SPECIES="Anura"     # closest available taxon in Dfam for Pristimantis

# RepeatMasker binary (adjust if not on PATH)
RM_BIN=$(which RepeatMasker 2>/dev/null || echo "/data/RepeatMasker/RepeatMasker")

# Input files (placed keeps original, unplaced is >=500bp filtered)
FILES=(
    "placed_contigs.fasta"
    "unplaced_contigs_filtered.fasta"
)

# ----- SETUP -----
mkdir -p "${OUTDIR}"

echo "============================================================"
echo "RepeatMasker batch run — $(date)"
echo "============================================================"
echo "Working directory: ${WORKDIR}"
echo "Output directory:  ${OUTDIR}"
echo "RepeatMasker:      ${RM_BIN}"
echo "Species:           ${SPECIES}"
echo "Threads (-pa):     ${THREADS}"
echo "Input files:       ${FILES[*]}"
echo "============================================================"
echo ""

# ----- CHECK INPUTS -----
for f in "${FILES[@]}"; do
    if [ ! -f "${WORKDIR}/${f}" ]; then
        echo "ERROR: Input file not found: ${WORKDIR}/${f}"
        exit 1
    fi
    echo "Found: ${f} ($(du -h "${WORKDIR}/${f}" | cut -f1))"
done

# Check RepeatMasker is available
if [ ! -x "${RM_BIN}" ]; then
    echo "ERROR: RepeatMasker not found at ${RM_BIN}"
    exit 1
fi

echo ""
echo "RepeatMasker version:"
${RM_BIN} -v
echo ""

# ----- Quick seqkit stats for the log -----
if command -v seqkit &> /dev/null; then
    echo "--- Input file statistics ---"
    for f in "${FILES[@]}"; do
        echo ""
        echo ">>> ${f}"
        seqkit stats -a "${WORKDIR}/${f}"
    done
    echo ""
fi

# ----- RUN REPEATMASKER -----
for f in "${FILES[@]}"; do
    BASENAME="${f%.fasta}"
    SET_OUTDIR="${OUTDIR}/${BASENAME}"
    mkdir -p "${SET_OUTDIR}"

    echo "============================================================"
    echo "[$(date)] Starting RepeatMasker on: ${f}"
    echo "  Output: ${SET_OUTDIR}"
    echo "============================================================"

    START_TIME=$(date +%s)

    ${RM_BIN} \
        -species "${SPECIES}" \
        -pa ${THREADS} \
        -gff \
        -xsmall \
        --uncurated \
        -dir "${SET_OUTDIR}" \
        "${WORKDIR}/${f}"

    END_TIME=$(date +%s)
    ELAPSED=$(( END_TIME - START_TIME ))
    HOURS=$(( ELAPSED / 3600 ))
    MINUTES=$(( (ELAPSED % 3600) / 60 ))

    echo ""
    echo "[$(date)] Finished ${f} in ${HOURS}h ${MINUTES}m"
    echo ""

    # Print the summary table
    TBL_FILE="${SET_OUTDIR}/${f}.tbl"
    if [ -f "${TBL_FILE}" ]; then
        echo "--- Summary for ${BASENAME} ---"
        cat "${TBL_FILE}"
        echo ""
    fi
done

# ----- COMBINED SUMMARY -----
echo "============================================================"
echo "ALL RUNS COMPLETE — $(date)"
echo "============================================================"
echo ""
echo "=== Side-by-side comparison of masking rates ==="
echo ""

for f in "${FILES[@]}"; do
    BASENAME="${f%.fasta}"
    TBL="${OUTDIR}/${BASENAME}/${f}.tbl"
    if [ -f "${TBL}" ]; then
        TOTAL=$(grep "total length:" "${TBL}" | head -1 | awk '{print $3}')
        GC=$(grep "GC level:" "${TBL}" | awk '{print $3}')
        MASKED=$(grep "bases masked:" "${TBL}" | awk '{print $3, $4}')
        echo "${BASENAME}:"
        echo "  Total length: ${TOTAL}"
        echo "  GC level:     ${GC}"
        echo "  Bases masked: ${MASKED}"
        echo ""
    fi
done

echo ""
echo "Detailed .tbl files in: ${OUTDIR}/<set_name>/"
echo "Use these for the fluorochrome banding analysis."
echo ""
echo "Done!"
