#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
set -euo pipefail

###############################################################################
# run_merqury.sh
#
# Reference-free assembly evaluation using Merqury (Rhie et al. 2020).
# Builds a meryl k-mer database from raw reads, then runs Merqury to
# estimate consensus QV, k-mer completeness, and spectral copy-number.
#
# Species:   Pristimantis euphronides (GrenadaFrog144)
# Assembly:  Flye v2.9.5 + Medaka v2.0.1 polished
# Reads:     ONT HAC-basecalled
# Mode:      Non-trio (single individual, no parental data)
#
# Requirements:
#   - meryl, merqury.sh, Rscript (conda env: merqury)
#   - R packages: ggplot2, scales, argparse
#
# Usage:
#   conda activate merqury
#   bash run_merqury.sh 2>&1 | tee merqury.log
#
# Runtime estimate: 2-4 hours (48 threads, 200 GB RAM, ~55 Gb ONT reads)
# Disk estimate:    ~50 GB (meryl database) + ~1 GB (Merqury outputs)
###############################################################################

# --- Configuration -----------------------------------------------------------
THREADS=48
KMER=21
MEMORY=200

BASEDIR="/home/ubuntu/GrenadaFrog144"
READS="${BASEDIR}/GrenadaFrog144_ONT_HAC_all.fastq.gz"
OUTDIR="${BASEDIR}/merqury_output"
MERYL_DB="${OUTDIR}/reads.k${KMER}.meryl"
MERYL_STATS="${OUTDIR}/meryl_statistics.txt"
SUMMARY="${OUTDIR}/merqury_summary_report.txt"
PREFIX="GrenadaFrog144"

# --- Locate assembly (prefer uncompressed) ------------------------------------
if [ -f "${BASEDIR}/final_medaka_polished_assembly_consensus.fasta" ]; then
    ASM="${BASEDIR}/final_medaka_polished_assembly_consensus.fasta"
elif [ -f "${BASEDIR}/final_medaka_polished_assembly_consensus.fasta.gz" ]; then
    ASM="${BASEDIR}/final_medaka_polished_assembly_consensus.fasta.gz"
else
    echo "ERROR: Assembly FASTA not found in ${BASEDIR}" >&2
    exit 1
fi

# --- Pre-flight checks -------------------------------------------------------
echo "=== Pre-flight checks ==="
echo "Date: $(date)"

for tool in meryl merqury.sh Rscript; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: ${tool} not found on PATH. Is the merqury env activated?" >&2
        exit 1
    fi
done

for f in "$READS" "$ASM"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: ${f}" >&2
        exit 1
    fi
done

AVAIL_GB=$(df --output=avail -BG "${BASEDIR}" | tail -1 | tr -dc '0-9')
echo "Available disk: ${AVAIL_GB} GB"
if [ "$AVAIL_GB" -lt 150 ]; then
    echo "WARNING: Only ${AVAIL_GB} GB free. Recommend >=150 GB." >&2
fi

echo "Reads:    ${READS}"
echo "Assembly: ${ASM}"
echo "k-mer:    ${KMER}"
echo "Threads:  ${THREADS}"
echo "Memory:   ${MEMORY} GB"
echo ""

mkdir -p "${OUTDIR}"

# --- Step 1: Build meryl k-mer database --------------------------------------
if [ -d "${MERYL_DB}" ]; then
    echo "=== Step 1: SKIP — meryl database already exists: ${MERYL_DB} ==="
else
    echo "=== Step 1: Building meryl k-mer database (k=${KMER}) ==="
    echo "Start: $(date)"

    meryl count k=${KMER} threads=${THREADS} \
        memory=${MEMORY} \
        "${READS}" \
        output "${MERYL_DB}"

    echo "Finished: $(date)"
    echo "Database size: $(du -sh "${MERYL_DB}" | cut -f1)"
fi
echo ""

# --- Step 2: K-mer statistics (saved to file) ---------------------------------
echo "=== Step 2: Meryl k-mer statistics ==="
meryl statistics "${MERYL_DB}" > "${MERYL_STATS}"
echo "Written to: ${MERYL_STATS}"
echo ""

# --- Step 3: Run Merqury ------------------------------------------------------
echo "=== Step 3: Running Merqury ==="
echo "Start: $(date)"

cd "${OUTDIR}"

# Decompress assembly temporarily if gzipped (Merqury needs uncompressed)
if [[ "${ASM}" == *.gz ]]; then
    ASM_INPUT="${OUTDIR}/assembly_tmp.fasta"
    if [ ! -f "${ASM_INPUT}" ]; then
        echo "Decompressing assembly..."
        gunzip -ck "${ASM}" > "${ASM_INPUT}"
    fi
else
    ASM_INPUT="${ASM}"
fi

merqury.sh "${MERYL_DB}" "${ASM_INPUT}" "${PREFIX}"

echo "Finished: $(date)"

# Clean up temp decompressed assembly
if [ -f "${OUTDIR}/assembly_tmp.fasta" ]; then
    rm "${OUTDIR}/assembly_tmp.fasta"
    echo "Cleaned up temporary decompressed assembly."
fi
echo ""

# --- Step 4: Summary report ---------------------------------------------------
echo "=== Step 4: Writing summary report ==="

{
    echo "============================================================================="
    echo "MERQURY ASSEMBLY EVALUATION — SUMMARY REPORT"
    echo "============================================================================="
    echo ""
    echo "Date:     $(date)"
    echo "Assembly: ${ASM}"
    echo "Reads:    ${READS}"
    echo "k-mer:    ${KMER}"
    echo "Threads:  ${THREADS}"
    echo "Memory:   ${MEMORY} GB"
    echo ""
    echo "============================================================================="
    echo "1. MERYL K-MER DATABASE"
    echo "============================================================================="
    echo ""
    echo "Database: ${MERYL_DB}"
    echo "Size on disk: $(du -sh "${MERYL_DB}" | cut -f1)"
    echo ""
    cat "${MERYL_STATS}"
    echo ""
    echo "============================================================================="
    echo "2. CONSENSUS QUALITY VALUE (QV)"
    echo "============================================================================="
    echo ""
    if [ -f "${OUTDIR}/${PREFIX}.qv" ]; then
        column -t "${OUTDIR}/${PREFIX}.qv"
    else
        echo "WARNING: QV file not found."
    fi
    echo ""
    echo "============================================================================="
    echo "3. K-MER COMPLETENESS"
    echo "============================================================================="
    echo ""
    if [ -f "${OUTDIR}/${PREFIX}.completeness.stats" ]; then
        column -t "${OUTDIR}/${PREFIX}.completeness.stats"
    else
        echo "WARNING: Completeness file not found."
    fi
    echo ""
    echo "============================================================================="
    echo "4. OUTPUT FILES"
    echo "============================================================================="
    echo ""
    ls -lh "${OUTDIR}"/${PREFIX}.* 2>/dev/null || echo "No output files found."
    echo ""
    echo "============================================================================="
    echo "5. DISK USAGE"
    echo "============================================================================="
    echo ""
    echo "Meryl database: $(du -sh "${MERYL_DB}" | cut -f1)"
    echo "Output dir:     $(du -sh "${OUTDIR}" | cut -f1)"
    echo ""
    echo "============================================================================="
    echo "Completed: $(date)"
    echo "============================================================================="
} > "${SUMMARY}" 2>&1

echo "Summary written to: ${SUMMARY}"
echo ""
cat "${SUMMARY}"
