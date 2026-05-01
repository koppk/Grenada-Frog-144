#!/bin/bash
# prepare_contig_sets.sh
# Author: Kopp K. Pristimantis euphronides genome project.
#
# Creates the prerequisite input files for all downstream contigComparison
# analyses from two inputs:
#   1. RagTag AGP (contig-to-scaffold placement)
#   2. Primary (unscaffolded) assembly FASTA
#
# Steps:
#   1. Extract placed/unplaced contig names from AGP
#   2. Extract FASTA subsets
#   3. Length-filter unplaced contigs (≥500 bp) for RepeatMasker
#   4. Per-contig length + GC% (seqkit fx2tab)
#   5. Summary statistics (seqkit stats)
#
# Output (all in OUTDIR):
#   placed_contig_names.txt
#   unplaced_contig_names.txt
#   placed_contigs.fasta
#   unplaced_contigs.fasta
#   unplaced_contigs_filtered.fasta     (≥500 bp for RepeatMasker)
#   placed_contigs.length_GC_cont.tsv
#   unplaced_contigs.length_GC_cont.tsv
#   placed_contigs.stats.tsv
#   unplaced_contigs.stats.tsv
#
# Dependencies: samtools, seqtk, seqkit
#
# Usage:
#   bash prepare_contig_sets.sh

set -euo pipefail

# === Paths ===
AGP="/data/GrenadaFrog144/ragtag.scaffold.renamed.agp"
ASM="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"
OUTDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"
MIN_LEN_CONTIGS=500  # minimum contig length for RepeatMasker

echo "=== Prepare contig sets ==="
echo "Start: $(date)"
echo "  AGP:        $AGP"
echo "  Assembly:   $ASM"
echo "  Output:     $OUTDIR"
echo ""

# === Sanity checks ===
for f in "$AGP" "$ASM"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

for cmd in samtools seqtk seqkit; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd not found on PATH"
        exit 1
    fi
done

mkdir -p "$OUTDIR"

# Index assembly if needed
if [ ! -f "${ASM}.fai" ]; then
    echo "Indexing assembly ..."
    samtools faidx "$ASM"
fi

# ============================================================
# Step 1: Extract placed/unplaced contig names from AGP
# ============================================================
echo "[Step 1] Extracting contig names from AGP ..."

# Placed: AGP type W rows where scaffold name starts with scaffold_
awk '!/^#/ && $5 == "W" && $1 ~ /^scaffold_/ { print $6 }' "$AGP" \
    | sort -u > "${OUTDIR}/placed_contig_names.txt"

# Unplaced: AGP type W rows where scaffold name starts with contig_
# (these are contigs RagTag could not place on any E. coqui chromosome)
awk '!/^#/ && $5 == "W" && $1 ~ /^contig_/ { print $6 }' "$AGP" \
    | sort -u > "${OUTDIR}/unplaced_contig_names.txt"

N_PLACED=$(wc -l < "${OUTDIR}/placed_contig_names.txt")
N_UNPLACED=$(wc -l < "${OUTDIR}/unplaced_contig_names.txt")
N_TOTAL=$(awk 'END{print NR}' "${ASM}.fai")
echo "  Placed contigs: $N_PLACED"
echo "  Unplaced contigs: $N_UNPLACED"
echo "  Total in assembly: $N_TOTAL"
echo "  Check: $((N_PLACED + N_UNPLACED)) (should equal $N_TOTAL)"
echo ""

# ============================================================
# Step 2: Extract FASTA subsets
# ============================================================
echo "[Step 2] Extracting FASTA subsets ..."

seqkit grep -f "${OUTDIR}/placed_contig_names.txt" "$ASM" \
    > "${OUTDIR}/placed_contigs.fasta"
echo "  placed_contigs.fasta: $(grep -c '^>' "${OUTDIR}/placed_contigs.fasta") contigs"

seqkit grep -f "${OUTDIR}/unplaced_contig_names.txt" "$ASM" \
    > "${OUTDIR}/unplaced_contigs.fasta"
echo "  unplaced_contigs.fasta: $(grep -c '^>' "${OUTDIR}/unplaced_contigs.fasta") contigs"
echo ""

# ============================================================
# Step 3: Length-filtered subset (>= MIN_LEN_CONTIGS bp)
# ============================================================
echo "[Step 3] Filtering unplaced contigs >= ${MIN_LEN_CONTIGS} bp ..."

seqkit seq -m "$MIN_LEN_CONTIGS" "${OUTDIR}/unplaced_contigs.fasta" \
    > "${OUTDIR}/unplaced_contigs_filtered.fasta" 2>/dev/null
N_UNPL_FILT=$(grep -c '^>' "${OUTDIR}/unplaced_contigs_filtered.fasta")
echo "  unplaced_contigs_filtered.fasta: $N_UNPL_FILT contigs ($((N_UNPLACED - N_UNPL_FILT)) dropped)"
echo ""

# ============================================================
# Step 4: Per-contig length and GC content
# ============================================================
echo "[Step 4] Computing per-contig length and GC% ..."

for LABEL in placed unplaced; do
    seqkit fx2tab -nlg "${OUTDIR}/${LABEL}_contigs.fasta" \
        > "${OUTDIR}/${LABEL}_contigs.length_GC_cont.tsv"
    echo "  ${LABEL}_contigs.length_GC_cont.tsv: $(wc -l < "${OUTDIR}/${LABEL}_contigs.length_GC_cont.tsv") contigs"
done
echo ""

# ============================================================
# Step 5: Summary statistics
# ============================================================
echo "[Step 5] Computing summary statistics ..."

for LABEL in placed unplaced; do
    seqkit stats -a "${OUTDIR}/${LABEL}_contigs.fasta" \
        > "${OUTDIR}/${LABEL}_contigs.stats.tsv"
    echo "  ${LABEL}_contigs.stats.tsv"
done
echo ""

# ============================================================
# Summary
# ============================================================
echo "=== Output files ==="
for f in placed_contig_names.txt unplaced_contig_names.txt \
         placed_contigs.fasta unplaced_contigs.fasta \
         unplaced_contigs_filtered.fasta \
         placed_contigs.length_GC_cont.tsv unplaced_contigs.length_GC_cont.tsv \
         placed_contigs.stats.tsv unplaced_contigs.stats.tsv; do
    [ -f "${OUTDIR}/$f" ] && echo "  OK   $f" || echo "  SKIP $f"
done

echo ""
echo "Finished: $(date)"
echo "=== Done ==="
