#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# run_inspector.sh
# ================
# Reference-free evaluation of the unscaffolded Medaka-polished Flye assembly
# using Inspector (Chen et al., Genome Biology 2021).
#
# Evaluates:
#   - Consensus quality (QV score)
#   - Structural errors (misjoins, collapses, expansions)
#   - Small-scale errors (substitutions, insertions, deletions)
#   - Read mapping statistics (mapping rate, splitting rate)
#
# Input:
#   - Assembly: Flye + Medaka polished contigs (unscaffolded)
#   - Reads: HAC basecalled ONT reads in FASTQ (with quality scores)
#
# Parameters:
#   --min_contig_length 3000        : evaluate contigs >= 3 kb for consensus
#       errors (lowered from default 10 kb to include small contigs; matches
#       QUAST min contig threshold)
#   --min_contig_length_assemblyerror 50000 : detect structural errors in
#       contigs >= 50 kb (lowered from default 1 Mb because contig N50 is
#       ~302 kb; default would skip the vast majority of the assembly.
#       Structural error detection needs sufficient length to be meaningful)
#   --datatype nanopore : sets p-value threshold to 0.05 for error detection
#
# Disk: Inspector creates read_to_contig.bam (~40-50 GB with quality scores)
#
# BEFORE RUNNING: conda activate /data/envs/inspector-env
#
# Date: 2026-02-20
# ===========================================================================

set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────
ASSEMBLY="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"
READS="/data/GrenadaFrog144/GrenadaFrog144_ONT_HAC_all.fastq.gz"
OUTDIR="/data/GrenadaFrog144/inspector_output"
THREADS=24

# ── Pre-flight checks ─────────────────────────────────────────────
echo "============================================================"
echo "Inspector: reference-free assembly evaluation"
echo "============================================================"
echo ""

errors=0

for f in "$ASSEMBLY" "$READS"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: File not found: $f"
        errors=$((errors + 1))
    else
        size=$(du -h "$f" | cut -f1)
        echo "  OK: $(basename "$f") ($size)"
    fi
done

if [[ $errors -gt 0 ]]; then
    echo ""
    echo "ABORTING: $errors input file(s) not found."
    exit 1
fi

echo ""

# ── Check disk space ──────────────────────────────────────────────
avail=$(df -BG /data | tail -1 | awk '{print $4}' | tr -d 'G')
echo "Available disk space on /data: ${avail} GB"
if [[ "$avail" -lt 60 ]]; then
    echo "WARNING: Less than 60 GB available. Inspector BAM may require ~40-50 GB."
    echo "         Consider freeing space before proceeding."
fi
echo ""

# ── Run Inspector ─────────────────────────────────────────────────
echo "Starting Inspector at $(date)"
echo "Output: $OUTDIR"
echo ""

inspector.py \
    -c "$ASSEMBLY" \
    -r "$READS" \
    --datatype nanopore \
    --min_contig_length 3000 \
    --min_contig_length_assemblyerror 50000 \
    -t "$THREADS" \
    -o "$OUTDIR"

echo ""
echo "Inspector finished at $(date)"
echo ""

# ── Print key results ─────────────────────────────────────────────
if [[ -f "$OUTDIR/summary_statistics" ]]; then
    echo "============================================================"
    echo "SUMMARY"
    echo "============================================================"
    cat "$OUTDIR/summary_statistics"
fi
