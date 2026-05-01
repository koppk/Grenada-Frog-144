#!/usr/bin/env bash
#
# 01_map_reads_and_compute_coverage.sh
#
# Maps HAC Oxford Nanopore reads to the primary scaffolded Pristimantis
# euphronides genome assembly using Minimap2 and computes genome-wide
# coverage in non-overlapping 10 kb windows using mosdepth.
#
# Corresponds to Additional file 2, section "Read mapping and genome
# coverage analysis", steps 1 (alignment) and 2a (windowed coverage).
#
# Input:
#   READS    : HAC basecalled Nanopore reads (FASTQ)
#   ASSEMBLY : Primary scaffolded P. euphronides genome (FASTA)
#
# Output:
#   eup_reads.bam                        : sorted, indexed BAM alignment
#   eup_reads.bam.bai                    : BAM index
#   eup_cov.regions.bed                  : 10 kb windowed coverage (BED4)
#   eup_cov.regions.bed.gz               : gzipped copy with tabix index
#   eup_cov.mosdepth.summary.txt         : per-scaffold and total summary
#   eup_cov.mosdepth.global.dist.txt     : cumulative coverage distribution
#   eup_cov.mosdepth.region.dist.txt     : per-region coverage distribution
#   eup_cov.per-base.bed.gz              : per-base coverage
#
# Prerequisites:
#   minimap2 v2.28   [Li, 2018]
#   samtools v1.20   [Danecek, 2021]
#   mosdepth v0.3.11 [Pedersen & Quinlan, 2018]
#
# Author: Kopp K, Pristimantis euphronides genome project

set -euo pipefail

# === Configuration ===
THREADS=24
READS="/data/GrenadaFrog144/GrenadaFrog144_ONT_HAC_all.fastq.gz"
ASSEMBLY="/data/GrenadaFrog144/Pristimantis_euphronides.genome.fasta"
OUTDIR="/data/GrenadaFrog144/coverage"
PREFIX="eup"

# === Verify inputs ===
for f in "$READS" "$ASSEMBLY"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Input file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# === Step 1: Align HAC reads to scaffolded assembly ===
echo "[$(date)] Step 1: Minimap2 alignment (map-ont preset)"
minimap2 -ax map-ont -t "$THREADS" "$ASSEMBLY" "$READS" | \
    samtools sort -@ "$THREADS" -o "${PREFIX}_reads.bam"

echo "[$(date)] Step 1: Indexing BAM"
samtools index -@ "$THREADS" "${PREFIX}_reads.bam"

# === Step 2a: Compute coverage in 10 kb non-overlapping windows ===
echo "[$(date)] Step 2a: mosdepth coverage computation (--by 10000)"
mosdepth --by 10000 -t "$THREADS" "${PREFIX}_cov" "${PREFIX}_reads.bam"

# === Decompress regions BED for downstream scripts ===
if [[ -f "${PREFIX}_cov.regions.bed.gz" && ! -f "${PREFIX}_cov.regions.bed" ]]; then
    echo "[$(date)] Decompressing regions BED"
    gunzip -k "${PREFIX}_cov.regions.bed.gz"
fi

# === Verify outputs ===
echo ""
echo "[$(date)] Output verification:"
for f in "${PREFIX}_reads.bam" \
         "${PREFIX}_reads.bam.bai" \
         "${PREFIX}_cov.regions.bed" \
         "${PREFIX}_cov.mosdepth.summary.txt" \
         "${PREFIX}_cov.mosdepth.global.dist.txt" \
         "${PREFIX}_cov.mosdepth.region.dist.txt"; do
    if [[ -f "$f" ]]; then
        echo "  OK: $f ($(du -h "$f" | cut -f1))"
    else
        echo "  MISSING: $f" >&2
    fi
done

echo ""
echo "[$(date)] Genome-wide summary:"
grep "^total" "${PREFIX}_cov.mosdepth.summary.txt"
echo ""
echo "[$(date)] Done."
