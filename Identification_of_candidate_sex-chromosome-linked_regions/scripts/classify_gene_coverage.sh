#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
set -euo pipefail

# classify_gene_coverage.sh
#
# Three-step pipeline to classify per-gene read coverage across the
# scaffolded P. euphronides genome assembly:
#
#   1. Extract gene coordinates from HANNO annotation (BED format)
#   2. Calculate per-gene coverage with mosdepth
#   3. Classify genes into coverage categories relative to genome-wide
#      median gene-body depth (total_region from mosdepth summary)
#
# Coverage classes:
#   Hemi_0.5x            0.375x – 0.625x of reference (hemizygous)
#   Auto_1.0x            0.75x  – 1.25x  of reference (autosomal)
#   High_Coverage/Repeat  > 1.25x of reference
#   Low_Coverage/Other    < 0.375x of reference
#
# Output is written to the specified output directory.
#
# Note on sequence naming in the scaffolded assembly:
# RagTag scaffolding assigns contigs from the input assembly to reference
# sequences. In our case, contigs placed onto E. coqui chromosomes (NC_*)
# or other E. coqui sequences (NW_*) receive names derived from the
# reference (e.g. NC_*_RagTag). These were renamed to scaffold_1 through
# scaffold_13 for clarity. Contigs that RagTag could not place remain as
# individual sequences with their original names (contig_*) from the
# Flye/Medaka assembly. RagTag does not renumber unplaced contigs — they
# retain their original identifiers. Therefore, by definition of the
# RagTag procedure, any sequence named contig_* in the scaffolded assembly
# was not assigned to a reference chromosome and can be considered unplaced.

THREADS=22

usage() {
    echo "Usage: $(basename "$0") <BESTMODELS_bedDB> <reads_bam> <output_dir>"
    echo ""
    echo "Arguments:"
    echo "  BESTMODELS_bedDB  HANNO annotation file (BESTMODELS-FINAL.renamed.bedDB)"
    echo "  reads_bam         BAM file of reads mapped to scaffolded assembly (e.g. eup_reads.bam)"
    echo "  output_dir        Directory for output files"
    echo ""
    echo "Output files (in output_dir):"
    echo "  all_genes_global_clean.bed              Gene coordinates (BED4)"
    echo "  all_genes_global_clean.regions.bed.gz    Per-gene coverage from mosdepth"
    echo "  all_genes_global_clean.mosdepth.summary.txt"
    echo "  genes_classified_robust.txt              Final classified gene list"
    echo ""
    echo "Requires: awk, mosdepth (or mosdepth_d4), samtools"
    exit 1
}

if [[ $# -ne 3 ]]; then
    usage
fi

BEDDB="$1"
BAM="$2"
OUTDIR="$3"

# Validate inputs
for f in "$BEDDB" "$BAM"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: File not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"

echo "=== Step 1: Extract gene coordinates from HANNO annotation ==="
echo "  Input:  $BEDDB"

awk -F'\t' 'NR>1 && $26 != "-" {print $1"\t"$2"\t"$3"\t"$26}' "$BEDDB" \
    > "${OUTDIR}/all_genes_global_clean.bed"

NGENES=$(wc -l < "${OUTDIR}/all_genes_global_clean.bed")
echo "  Output: ${OUTDIR}/all_genes_global_clean.bed (${NGENES} genes)"

echo ""
echo "=== Step 2: Calculate per-gene coverage with mosdepth ==="
echo "  Input:  $BAM"

# Use mosdepth_d4 if available, otherwise mosdepth
MOSDEPTH="mosdepth"
if command -v mosdepth_d4 &>/dev/null; then
    MOSDEPTH="mosdepth_d4"
fi

"$MOSDEPTH" \
    -t "$THREADS" \
    -n \
    -x \
    --by "${OUTDIR}/all_genes_global_clean.bed" \
    "${OUTDIR}/all_genes_global_clean" \
    "$BAM"

echo "  Output: ${OUTDIR}/all_genes_global_clean.regions.bed.gz"
echo "  Output: ${OUTDIR}/all_genes_global_clean.mosdepth.summary.txt"

echo ""
echo "=== Step 3: Classify genes by coverage ==="

# Reference value: total_region mean from mosdepth summary
REF=$(grep "^total_region" "${OUTDIR}/all_genes_global_clean.mosdepth.summary.txt" \
    | awk '{print $4}')
echo "  Reference coverage (total_region): ${REF}x"

zcat "${OUTDIR}/all_genes_global_clean.regions.bed.gz" \
    | awk -F'\t' -v ref="$REF" '{
    cov = $5
    if (cov >= ref*0.375 && cov <= ref*0.625) class = "Hemi_0.5x"
    else if (cov >= ref*0.75 && cov <= ref*1.25) class = "Auto_1.0x"
    else if (cov > ref*1.25) class = "High_Coverage/Repeat"
    else class = "Low_Coverage/Other"
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"class
}' > "${OUTDIR}/genes_classified_robust.txt"

# Summary
echo ""
echo "=== Classification summary ==="
awk -F'\t' '{c[$6]++} END {for (k in c) printf "  %-25s %d\n", k, c[k]}' \
    "${OUTDIR}/genes_classified_robust.txt" | sort

TOTAL=$(wc -l < "${OUTDIR}/genes_classified_robust.txt")
echo "  -------------------------"
printf "  %-25s %d\n" "Total" "$TOTAL"

# Write summary to file
{
    printf "Coverage_Class\tGene_Count\n"
    awk -F'\t' '{c[$6]++} END {for (k in c) printf "%s\t%d\n", k, c[k]}' \
        "${OUTDIR}/genes_classified_robust.txt" | sort
    printf "Total\t%d\n" "$TOTAL"
} > "${OUTDIR}/genes_classified_robust_summary.tsv"

echo ""
echo "Done."
echo "  ${OUTDIR}/genes_classified_robust.txt"
echo "  ${OUTDIR}/genes_classified_robust_summary.tsv"
