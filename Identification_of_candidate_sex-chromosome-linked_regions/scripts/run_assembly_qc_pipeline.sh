#!/bin/bash
#
# run_assembly_qc_pipeline.sh
# ============================
# Wrapper script for the GrenadaFrog144 primary assembly QC pipeline.
# Runs both analysis scripts in order and produces a combined summary
# suitable for manuscript methods/results and supplementary information.
#
# Pipeline:
#   1. map_reads_primary_assembly_coverage.sh
#      - Maps ONT reads to unscaffolded assembly
#      - Computes per-contig coverage, GC content, read lengths
#      - Classifies placed/unplaced contigs by coverage band
#
#   2. primary_assembly_bam_qc.sh
#      - Assembly contiguity (N50, L50, etc.)
#      - Coverage uniformity (CV per contig)
#      - Mapping quality (MAPQ per contig)
#      - Base quality (Phred per contig)
#      - Heterozygosity (variant density per contig)
#      - Soft-clipping and supplementary rates
#      - Cross-metric quality flags
#
# Dependencies:
#   minimap2, samtools, bcftools, mosdepth, seqkit, parallel, python3
#   Install: mamba install -c bioconda -c conda-forge \
#            minimap2 samtools bcftools mosdepth seqkit parallel
#
# Usage:
#   nohup bash run_assembly_qc_pipeline.sh &> pipeline.log &
#
# Output:
#   Combined summary: /data/GrenadaFrog144/assembly_qc/assembly_qc_combined_summary.txt
#   Individual summaries preserved in their respective output directories.
#

# Author: Kopp K. Pristimantis euphronides genome project.
set -euo pipefail

# Locate scripts relative to this wrapper
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COMBINED_SUMMARY="/data/GrenadaFrog144/assembly_qc/assembly_qc_combined_summary.txt"

echo "========================================================"
echo "  GrenadaFrog144 Assembly QC Pipeline"
echo "  $(date)"
echo "========================================================"
echo ""
echo "Scripts: $SCRIPTDIR"
echo ""

# Check both scripts exist
for script in map_reads_primary_assembly_coverage.sh primary_assembly_bam_qc.sh; do
    if [ ! -f "${SCRIPTDIR}/${script}" ]; then
        echo "ERROR: ${SCRIPTDIR}/${script} not found"
        echo "       All three scripts must be in the same directory."
        exit 1
    fi
done

# ============================================================
# Part 1: Read mapping and coverage analysis
# ============================================================
echo "========================================================"
echo "  PART 1: Read mapping and coverage analysis"
echo "========================================================"
echo ""

bash "${SCRIPTDIR}/map_reads_primary_assembly_coverage.sh"

echo ""
echo "Part 1 complete."
echo ""

# ============================================================
# Part 2: BAM-based assembly QC
# ============================================================
echo "========================================================"
echo "  PART 2: BAM-based assembly QC"
echo "========================================================"
echo ""

bash "${SCRIPTDIR}/primary_assembly_bam_qc.sh"

echo ""
echo "Part 2 complete."
echo ""

# ============================================================
# Part 3: Combined summary for manuscript
# ============================================================
echo "========================================================"
echo "  Generating combined summary"
echo "========================================================"

COV_SUMMARY="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/coverage/coverage_three_set_summary.txt"
QC_SUMMARY="/data/GrenadaFrog144/assembly_qc/assembly_qc_summary.txt"

{
    echo "================================================================"
    echo "  Pristimantis euphronides (GrenadaFrog144)"
    echo "  Primary Assembly Quality Assessment — Combined Report"
    echo "  Generated: $(date)"
    echo "================================================================"
    echo ""
    echo "Sequencing: Oxford Nanopore Technologies (ONT) R10.4.1, HAC basecalling"
    echo "Assembly:   Flye + Medaka polishing (unscaffolded primary assembly)"
    echo "Scaffolding reference: Eleutherodactylus coqui (cross-species, RagTag)"
    echo ""
    echo ""

    echo "################################################################"
    echo "  PART 1: READ MAPPING AND COVERAGE ANALYSIS"
    echo "################################################################"
    echo ""
    if [ -s "$COV_SUMMARY" ]; then
        cat "$COV_SUMMARY"
    else
        echo "  WARNING: $COV_SUMMARY not found or empty"
    fi
    echo ""
    echo ""

    echo "################################################################"
    echo "  PART 2: ASSEMBLY QUALITY METRICS"
    echo "################################################################"
    echo ""
    if [ -s "$QC_SUMMARY" ]; then
        cat "$QC_SUMMARY"
    else
        echo "  WARNING: $QC_SUMMARY not found or empty"
    fi
    echo ""
    echo ""

    echo "################################################################"
    echo "  ALL OUTPUT FILES"
    echo "################################################################"
    echo ""
    echo "  Coverage analysis:"
    COVDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/coverage"
    echo "    ${COVDIR}/primary_assembly.bam"
    echo "    ${COVDIR}/primary_assembly.bam.bai"
    echo "    ${COVDIR}/placed_coverage.tsv"
    echo "    ${COVDIR}/placed_coverage_classified.tsv"
    echo "    ${COVDIR}/unplaced_coverage.tsv"
    echo "    ${COVDIR}/unplaced_coverage_classified.tsv"
    echo "    ${COVDIR}/unmappable_zero_cov_contigs.tsv"
    echo "    ${COVDIR}/missing_zero_cov.txt"
    echo "    ${COVDIR}/samtools_stats.txt"
    echo "    ${COVDIR}/contig_gc_content.tsv"
    echo "    ${COVDIR}/contig_read_lengths.tsv"
    echo "    ${COVDIR}/coverage_three_set_summary.txt"
    echo ""
    echo "  Assembly QC:"
    QCDIR="/data/GrenadaFrog144/assembly_qc"
    echo "    ${QCDIR}/assembly_contiguity.txt"
    echo "    ${QCDIR}/contig_coverage_uniformity.tsv"
    echo "    ${QCDIR}/contig_mapq.tsv"
    echo "    ${QCDIR}/contig_base_quality.tsv"
    echo "    ${QCDIR}/variants.vcf.gz"
    echo "    ${QCDIR}/contig_heterozygosity.tsv"
    echo "    ${QCDIR}/contig_clipping_supplementary.tsv"
    echo "    ${QCDIR}/flagged_contigs.tsv"
    echo "    ${QCDIR}/assembly_qc_summary.txt"
    echo ""
    echo "  Combined summary:"
    echo "    $COMBINED_SUMMARY"
    echo ""
} 2>&1 | tee "$COMBINED_SUMMARY"

echo "  Combined summary: $COMBINED_SUMMARY"
echo ""
echo "========================================================"
echo "  Pipeline complete: $(date)"
echo "========================================================"
