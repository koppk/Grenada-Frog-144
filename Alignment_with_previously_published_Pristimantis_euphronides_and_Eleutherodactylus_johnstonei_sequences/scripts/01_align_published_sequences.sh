#!/usr/bin/env bash
#
# 01_align_published_sequences.sh
#
# Aligns the three publicly available P. euphronides GenBank sequences
# to the Medaka-polished Flye assembly using Minimap2, extracts the
# top-matching contigs, and prepares subsequences for BLASTn validation.
#
# NCBI web BLASTn steps (not scriptable) are documented at the end
# of this file.
#
# Input:
#   input/EF493427_PriEuph_gene.fasta       RAG-1 partial CDS (578 bp)
#   input/EF493489_PriEuph_gene.fasta       Tyr exon 1 partial CDS (493 bp)
#   input/EF493527_PriEuph_gene.fasta       mt tRNA-Phe/12S/tRNA-Val/16S (2570 bp)
#   Assembly: /data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta
#
# Output (in output/):
#   EF493427_PriEuph_final_medaka_polished_assembly.paf
#   EF493489_PriEuph_final_medaka_polished_assembly.paf
#   EF493527_PriEuph_final_medaka_polished_assembly.paf
#   contig_list.txt
#   mapping_contigs.fasta
#   mapping_contigs.fasta.split/mapping_contigs.part_*.fasta
#   contig_*_part_that_aligned_to_*.fasta
#
# Prerequisites:
#   minimap2 v2.28  [Li, 2018]
#   samtools v1.20  [Danecek, 2021]
#   seqkit          [Shen, 2024]
#
# Author: Kopp K, Pristimantis euphronides genome project

set -euo pipefail

BASEDIR="$(cd "$(dirname "$0")/.." && pwd)"
INPUT="${BASEDIR}/input"
OUTPUT="${BASEDIR}/output"
ASSEMBLY="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"

mkdir -p "${OUTPUT}"

# === Step 1: Minimap2 alignment of each reference to the assembly ===

for acc in EF493427 EF493489 EF493527; do
    query="${INPUT}/${acc}_PriEuph_gene.fasta"
    paf="${OUTPUT}/${acc}_PriEuph_final_medaka_polished_assembly.paf"

    echo "[$(date)] Aligning ${acc} to assembly ..."
    minimap2 -x map-ont "${ASSEMBLY}" "${query}" > "${paf}"
    echo "  Output: ${paf}"
done

# === Step 2: Extract top-matching contig name from each PAF ===
# The top hit per query is the first line (highest alignment score)

> "${OUTPUT}/contig_list.txt"
for paf in "${OUTPUT}"/*_PriEuph_final_medaka_polished_assembly.paf; do
    # Column 6 = target name; take first line (best hit)
    head -1 "${paf}" | cut -f6 >> "${OUTPUT}/contig_list.txt"
done

sort -u "${OUTPUT}/contig_list.txt" -o "${OUTPUT}/contig_list.txt"
echo ""
echo "Top-matching contigs:"
cat "${OUTPUT}/contig_list.txt"

# === Step 3: Extract these contigs from the assembly ===

seqkit grep -f "${OUTPUT}/contig_list.txt" "${ASSEMBLY}" \
    > "${OUTPUT}/mapping_contigs.fasta"

# === Step 4: Split into per-contig FASTAs ===

seqkit split -i "${OUTPUT}/mapping_contigs.fasta" \
    -O "${OUTPUT}/mapping_contigs.fasta.split"

# === Step 5: Extract aligned portions of contigs ===
# Parse PAF coordinates (col 6=target, col 8=tstart, col 9=tend)
# and extract the corresponding subsequence from each contig.

extract_aligned_region() {
    local paf="$1"
    local acc="$2"
    local outfasta="$3"

    # Best hit = first line
    local contig start end
    contig=$(head -1 "${paf}" | cut -f6)
    start=$(head -1 "${paf}" | cut -f8)
    end=$(head -1 "${paf}" | cut -f9)

    echo "  ${acc}: ${contig}:${start}-${end}"
    samtools faidx "${OUTPUT}/mapping_contigs.fasta" \
        "${contig}:$((start+1))-${end}" > "${outfasta}"
}

echo ""
echo "Extracting aligned regions ..."
samtools faidx "${OUTPUT}/mapping_contigs.fasta"

extract_aligned_region \
    "${OUTPUT}/EF493427_PriEuph_final_medaka_polished_assembly.paf" \
    "EF493427" \
    "${OUTPUT}/contig_15875_part_that_aligned_to_EF493427_PriEuph.fasta"

extract_aligned_region \
    "${OUTPUT}/EF493489_PriEuph_final_medaka_polished_assembly.paf" \
    "EF493489" \
    "${OUTPUT}/contig_19362_part_that_aligned_to_EF493489_PriEuph.fasta"

extract_aligned_region \
    "${OUTPUT}/EF493527_PriEuph_final_medaka_polished_assembly.paf" \
    "EF493527" \
    "${OUTPUT}/contig_57691_part_that_aligned_to_EF493527_PriEuph.fasta"

echo ""
echo "[$(date)] Local steps complete."
echo ""
echo "================================================================"
echo "MANUAL STEPS (NCBI web BLASTn interface)"
echo "================================================================"
echo ""
echo "Step 6: BLASTn pairwise alignment vs P. euphronides references"
echo "  Tool: NCBI BLASTn (https://blast.ncbi.nlm.nih.gov/)"
echo "  Mode: 'Align two or more sequences'"
echo "  For each locus:"
echo "    Query:   per-contig FASTA from mapping_contigs.fasta.split/"
echo "    Subject: respective GenBank reference (EF493427/489/527)"
echo "  Save pairwise alignment as:"
echo "    EF493427_PriEuph_gene_contig_15875_*_Alignment.txt"
echo "    EF493489_PriEuph_gene_contig_19362_*_Alignment.txt"
echo "    EF493527_PriEuph_gene_contig_57691_*_Alignment.txt"
echo ""
echo "Step 7: BLASTn vs E. johnstonei (taxid:350008)"
echo "  Tool: NCBI BLASTn"
echo "  Database: Nucleotide collection (nr/nt)"
echo "  Organism: Eleutherodactylus johnstonei (taxid:350008)"
echo "  Query: aligned-portion FASTAs from Step 5"
echo "    (Exception: for RAG-1, use the full contig_15875 as query)"
echo "  Save top-scoring alignment as:"
echo "    OM914617_EleJoh_gene_contig_57691_*_Alignment.txt"
echo "    OM928401_EleJoh_gene_contig_19362_*_Alignment.txt"
echo "    JX298190_EleJoh_gene_contig_15875_*_Alignment.txt"
echo ""
echo "Step 8: Compile summary table from alignment results"
echo "================================================================"
