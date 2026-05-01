#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# run_miniprot_both_sets.sh
# =========================
# Run miniprot with Anuran UniProt (taxonomy 8292) on placed and unplaced
# contig sets for comparable gene content assessment.
#
# Complementary to HANNO (E. coqui proteome):
#   - HANNO detects orthologs to the closest reference
#   - miniprot casts a wider net across Anura
#   - Contigs with miniprot hits but no HANNO hits → divergent W-gametolog candidates
#
# Input:
#   placed_contigs.fasta   (8,815 contigs, 1.48 Gb)
#   unplaced_contigs.fasta (25,291 contigs, 271 Mb — full set incl. <500 bp)
#   Anuran UniProt reviewed proteins (taxonomy 8292)
#
# Output:
#   placed_contigs_amphibia_miniprot.gff
#   unplaced_contigs_amphibia_miniprot.gff
#
# Usage:
#   bash run_miniprot_both_sets.sh
#
# Date: 2026-02-18

set -euo pipefail

BASEDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"
PROTDB="/data/GrenadaFrog144/SexChromosomes/Z-chr-workflow/uniprotkb_taxonomy_id_8292_AND_reviewed_2026_02_12.fasta"
THREADS=8

cd "$BASEDIR"

echo "=== Miniprot: Anuran UniProt on placed + unplaced contigs ==="
echo "Start: $(date)"
echo ""

# === Sanity checks ===
for f in placed_contigs.fasta unplaced_contigs.fasta "$PROTDB"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

if ! command -v miniprot &> /dev/null; then
    echo "ERROR: miniprot not found on PATH"
    exit 1
fi

echo "  Protein DB: $PROTDB"
echo "  Proteins:   $(grep -c '^>' "$PROTDB")"
echo "  Threads:    $THREADS"
echo ""

# === Placed contigs ===
echo "[1/2] Running miniprot on placed contigs ..."
echo "  Input: placed_contigs.fasta ($(grep -c '^>' placed_contigs.fasta) contigs)"
T0=$(date +%s)

miniprot -t "$THREADS" --gff -G 500k -j 1 \
    placed_contigs.fasta \
    "$PROTDB" \
    > placed_contigs_amphibia_miniprot.gff

T1=$(date +%s)
echo "  Done in $((T1 - T0)) seconds"
echo "  mRNA hits: $(grep -c 'mRNA' placed_contigs_amphibia_miniprot.gff)"
echo "  Unique proteins: $(grep 'mRNA' placed_contigs_amphibia_miniprot.gff | sed 's/.*Target=//;s/ .*//' | sort -u | wc -l)"
echo ""

# === Unplaced contigs (full set) ===
echo "[2/2] Running miniprot on unplaced contigs (full set) ..."
echo "  Input: unplaced_contigs.fasta ($(grep -c '^>' unplaced_contigs.fasta) contigs)"
T0=$(date +%s)

miniprot -t "$THREADS" --gff -G 500k -j 1 \
    unplaced_contigs.fasta \
    "$PROTDB" \
    > unplaced_contigs_amphibia_miniprot.gff

T1=$(date +%s)
echo "  Done in $((T1 - T0)) seconds"
echo "  mRNA hits: $(grep -c 'mRNA' unplaced_contigs_amphibia_miniprot.gff)"
echo "  Unique proteins: $(grep 'mRNA' unplaced_contigs_amphibia_miniprot.gff | sed 's/.*Target=//;s/ .*//' | sort -u | wc -l)"
echo ""

echo "=== Output files ==="
echo "  ${BASEDIR}/placed_contigs_amphibia_miniprot.gff"
echo "  ${BASEDIR}/unplaced_contigs_amphibia_miniprot.gff"
echo ""
echo "Finished: $(date)"
echo "=== Done ==="
