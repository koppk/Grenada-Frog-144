#!/bin/bash
set -euo pipefail

# z_candidate_workflow.sh
#
# Automated Z-chromosome candidate detection from mosdepth coverage data.
#
# Layer 1: Statistical scaffold screening (binomial test, Bonferroni-corrected)
# Layer 2: Per-Mb coverage profile extraction for flagged scaffolds
# Layer 3: Boundary detection via optimal binary segmentation + 2 SD extension
# Layer 4: Haplotype divergence comparison (Z-candidate vs autosomal regions)
#
# Requires in same directory (SCRIPTDIR):
#   identify_z_candidate_regions.py   (Layer 1)
#   detect_boundaries.py              (Layer 3)
#   haplotype_identity.sh             (Layer 4)
#
# Usage:
#   bash z_candidate_workflow.sh <eup_cov.regions.bed.gz> <output_dir> \
#       <hap1.fasta> <hap2.fasta> [threads]
#
# Example:
#   cd /data/GrenadaFrog144/SexChromosomes/Z-chr-workflow
#   bash scripts/z_candidate_workflow.sh \
#       /data/GrenadaFrog144/coverage/eup_cov.regions.bed.gz \
#       ./z_candidate_detection/ \
#       /data/GrenadaFrog144/RagTag_postHapDup/ragtag_scaffold_dual_1/PriEup.hap1_genomic.fasta \
#       /data/GrenadaFrog144/RagTag_postHapDup/ragtag_scaffold_dual_2/PriEup.hap2_genomic.fasta
#
# Author: Kopp K., Pristimantis euphronides genome project

BEDGZ="${1:?Usage: $0 <eup_cov.regions.bed.gz> <output_dir> <hap1.fasta> <hap2.fasta> [threads]}"
OUTDIR="${2:?Usage: $0 <eup_cov.regions.bed.gz> <output_dir> <hap1.fasta> <hap2.fasta> [threads]}"
HAP1="${3:?Usage: $0 <eup_cov.regions.bed.gz> <output_dir> <hap1.fasta> <hap2.fasta> [threads]}"
HAP2="${4:?Usage: $0 <eup_cov.regions.bed.gz> <output_dir> <hap1.fasta> <hap2.fasta> [threads]}"
THREADS="${5:-22}"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUTDIR"

echo "=== Z-candidate workflow ==="
echo "Input: $BEDGZ"
echo "Output: $OUTDIR"
echo "HAP1: $HAP1"
echo "HAP2: $HAP2"
echo "Threads: $THREADS"
echo ""

# ── Layer 1: Statistical scaffold screening ────────────────────
echo "[Layer 1] Running identify_z_candidate_regions.py ..."
python3 "${SCRIPTDIR}/identify_z_candidate_regions.py" "$BEDGZ" "$OUTDIR"
echo ""

# Read flagged scaffolds
FLAGGED=$(awk -F'\t' 'NR>1 && $NF=="YES" {print $1}' "${OUTDIR}/scaffold_screening_summary.tsv")

if [ -z "$FLAGGED" ]; then
    echo "No scaffolds flagged. Exiting."
    exit 0
fi

# ── Layer 2: Per-Mb coverage profiles ──────────────────────────
echo "[Layer 2] Extracting per-Mb coverage profiles..."
echo ""

for SCAFF in $FLAGGED; do
    OUTFILE="${OUTDIR}/${SCAFF}_per_mb_coverage.tsv"
    echo "  ${SCAFF} -> ${OUTFILE}"
    printf "Mb\tmean_coverage\n" > "$OUTFILE"
    zcat "$BEDGZ" | awk -v s="$SCAFF" '$1==s' | \
      awk '{mb=int($2/1000000); sum[mb]+=$4; n[mb]++}
           END {for(i=0;i<=500;i++) if(n[i]>0) printf "%d\t%.2f\n", i, sum[i]/n[i]}' | \
      sort -n >> "$OUTFILE"
done

echo ""

# ── Layer 3: Boundary detection ────────────────────────────────
echo "[Layer 3] Detecting half-coverage block boundaries..."
echo "  Method: optimal binary segmentation (2 vs 3 segments, BIC)"
echo "  Outlier capping: 2x median per scaffold"
echo "  Extension: flanking mean - 2 SD threshold"
echo ""

REGIONS="${OUTDIR}/z_candidate_regions.tsv"
printf "scaffold\tstart_Mb\tend_Mb\tlength_Mb\tmean_cov_region\tmean_cov_flanking\tratio\tmodel\n" > "$REGIONS"

for SCAFF in $FLAGGED; do
    INFILE="${OUTDIR}/${SCAFF}_per_mb_coverage.tsv"
    python3 "${SCRIPTDIR}/detect_boundaries.py" "$SCAFF" "$INFILE" "$REGIONS"
done

# ── Layer 3 summary ───────────────────────────────────────────
echo "=== Z-candidate regions ==="
echo ""
tail -n +2 "$REGIONS" | while IFS=$'\t' read -r sc st en ln mc mf rt md; do
    echo "  ${sc}: ${st} - ${en} Mb  (${ln} Mb)"
done
echo ""
TOTAL=$(tail -n +2 "$REGIONS" | awk -F'\t' '{s+=$4} END {print s}')
echo "  Total candidate Z-linked sequence: ${TOTAL} Mb"
echo ""

# ── Layer 4: Haplotype divergence ──────────────────────────────
echo "[Layer 4] Running haplotype identity comparison..."
echo "  (this may take several minutes per region)"
echo ""
bash "${SCRIPTDIR}/haplotype_identity.sh" \
    "$REGIONS" \
    "${OUTDIR}/scaffold_screening_summary.tsv" \
    "$HAP1" "$HAP2" "$OUTDIR" "$THREADS"

# ── Output summary ─────────────────────────────────────────────
echo "Output files:"
echo "  ${OUTDIR}/scaffold_screening_summary.tsv  (Layer 1)"
for SCAFF in $FLAGGED; do
    echo "  ${OUTDIR}/${SCAFF}_per_mb_coverage.tsv  (Layer 2)"
done
echo "  ${REGIONS}  (Layer 3)"
echo "  ${OUTDIR}/haplotype_identity.tsv  (Layer 4)"
echo "  ${OUTDIR}/mann_whitney_identity.tsv  (Layer 4)"
