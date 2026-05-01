#!/bin/bash
set -euo pipefail

# haplotype_identity.sh
#
# Layer 4: Haplotype divergence in Z-candidate vs autosomal regions.
#
# For each Z-candidate region identified in Layer 3, extracts the
# corresponding region from both RagTag-scaffolded haplotype assemblies
# and computes pairwise identity via minimap2 asm5 alignment.
# Compares Z-candidate regions against autosomal flanking regions on the
# same scaffold and against all non-flagged scaffolds (full-length).
# Tests significance via exact one-tailed Mann-Whitney U test.
#
# Reads boundaries from z_candidate_regions.tsv (Layer 3 output).
# Scaffold naming: scaffold_N -> CM0694{04+N}.1_RagTag
#
# Requires in same directory (SCRIPTDIR):
#   mann_whitney_identity.py
#
# Usage:
#   bash haplotype_identity.sh <z_candidate_regions.tsv> \
#       <scaffold_screening_summary.tsv> \
#       <hap1.fasta> <hap2.fasta> <output_dir> [threads]
#
# Author: Kopp K., Pristimantis euphronides genome project

REGIONS="${1:?Usage: $0 <z_candidate_regions.tsv> <scaffold_screening_summary.tsv> <hap1.fasta> <hap2.fasta> <output_dir> [threads]}"
SCREENING="${2:?Usage: $0 <z_candidate_regions.tsv> <scaffold_screening_summary.tsv> <hap1.fasta> <hap2.fasta> <output_dir> [threads]}"
HAP1="${3:?Usage: $0 ... <hap1.fasta> <hap2.fasta> <output_dir> [threads]}"
HAP2="${4:?Usage: $0 ... <hap1.fasta> <hap2.fasta> <output_dir> [threads]}"
OUTDIR="${5:?Usage: $0 ... <output_dir> [threads]}"
THREADS="${6:-22}"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUTDIR"
OUTFILE="${OUTDIR}/haplotype_identity.tsv"
TMPDIR="${OUTDIR}/tmp_hap_identity"
mkdir -p "$TMPDIR"

echo "=== Haplotype identity comparison (Layer 4) ==="
echo "  Regions:    $REGIONS"
echo "  Screening:  $SCREENING"
echo "  HAP1:       $HAP1"
echo "  HAP2:       $HAP2"
echo "  Threads:    $THREADS"
echo ""

# Ensure FAI indices exist
[ ! -f "${HAP1}.fai" ] && samtools faidx "$HAP1"
[ ! -f "${HAP2}.fai" ] && samtools faidx "$HAP2"

# Scaffold_N -> CM accession in RagTag haplotype assemblies
scaffold_to_cm() {
    local scaff_num="${1#scaffold_}"
    printf "CM0%d.1_RagTag" "$((69404 + scaff_num))"
}

# Get scaffold length from FAI
get_length() {
    awk -v s="$1" '$1==s {print $2}' "${HAP1}.fai"
}

# Compute mean pairwise identity between two FASTA regions
# Args: label type region_string
compute_identity() {
    local label="$1"
    local type="$2"
    local region="$3"

    samtools faidx "$HAP1" "$region" > "${TMPDIR}/hap1.fa"
    samtools faidx "$HAP2" "$region" > "${TMPDIR}/hap2.fa"

    minimap2 -t "$THREADS" -cx asm5 "${TMPDIR}/hap1.fa" "${TMPDIR}/hap2.fa" 2>/dev/null | \
        awk -v lab="$label" -v typ="$type" \
        '{sum += $10/$11; n++}
         END {
             if (n > 0) printf "%s\t%s\t%.6f\t%d\n", lab, typ, sum/n, n;
             else        printf "%s\t%s\tNA\t0\n", lab, typ
         }'
}

# Header
printf "Region\tType\tIdentity\tN_alignments\n" > "$OUTFILE"

# ── Flagged scaffolds: Z-candidate + autosomal flanks ──────────
echo "  Processing Z-candidate and flanking autosomal regions..."
echo ""

# Read Z-candidate regions
declare -A flagged
while IFS=$'\t' read -r scaff start_mb end_mb length_mb _rest; do
    [ "$scaff" = "scaffold" ] && continue
    flagged["$scaff"]=1

    cm=$(scaffold_to_cm "$scaff")
    scaff_num="${scaff#scaffold_}"
    scaff_len=$(get_length "$cm")

    # Mb to bp (1-based inclusive, full Mb bins)
    start_bp=$((start_mb * 1000000 + 1))
    end_bp=$(( (end_mb + 1) * 1000000 ))
    [ "$end_bp" -gt "$scaff_len" ] && end_bp="$scaff_len"

    echo "    ${scaff} (${cm}, ${scaff_len} bp)"
    echo "      Z-candidate: ${start_mb}-${end_mb} Mb -> ${start_bp}-${end_bp} bp"

    # Z-candidate region
    compute_identity \
        "Chr${scaff_num}_0.5x_${start_mb}-${end_mb}Mb" \
        "Z_candidate" \
        "${cm}:${start_bp}-${end_bp}" >> "$OUTFILE"

    # Autosomal flanking region(s)
    if [ "$start_mb" -eq 0 ]; then
        # Terminal block at start: autosomal flank after Z region
        auto_start=$((end_bp + 1))
        scaff_len_mb=$((scaff_len / 1000000))
        echo "      Autosomal:    ${end_mb}-${scaff_len_mb} Mb"
        compute_identity \
            "Chr${scaff_num}_auto_${end_mb}-${scaff_len_mb}Mb" \
            "Autosomal" \
            "${cm}:${auto_start}-${scaff_len}" >> "$OUTFILE"
    else
        # Internal block: autosomal flanks on both sides
        auto1_end=$((start_bp - 1))
        echo "      Autosomal:    0-${start_mb} Mb"
        compute_identity \
            "Chr${scaff_num}_auto_0-${start_mb}Mb" \
            "Autosomal" \
            "${cm}:1-${auto1_end}" >> "$OUTFILE"

        auto2_start=$((end_bp + 1))
        if [ "$auto2_start" -lt "$scaff_len" ]; then
            scaff_len_mb=$((scaff_len / 1000000))
            echo "      Autosomal:    ${end_mb}-${scaff_len_mb} Mb"
            compute_identity \
                "Chr${scaff_num}_auto_${end_mb}-${scaff_len_mb}Mb" \
                "Autosomal" \
                "${cm}:${auto2_start}-${scaff_len}" >> "$OUTFILE"
        fi
    fi
    echo ""
done < "$REGIONS"

# ── Non-flagged scaffolds: full-length autosomal reference ─────
echo "  Processing non-flagged scaffolds (full-length autosomal)..."

# Read all scaffolds from screening summary
while IFS=$'\t' read -r scaff _rest; do
    [ "$scaff" = "scaffold" ] && continue
    [ -n "${flagged[$scaff]+x}" ] && continue

    cm=$(scaffold_to_cm "$scaff")
    scaff_num="${scaff#scaffold_}"
    echo "    ${scaff} (${cm})"

    compute_identity "Chr${scaff_num}_full" "Autosomal" "$cm" >> "$OUTFILE"
done < "$SCREENING"

# Clean up
rm -rf "$TMPDIR"

# ── Summary ────────────────────────────────────────────────────
echo ""
echo "=== Haplotype identity results ==="
echo ""
column -t -s$'\t' "$OUTFILE"
echo ""

# ── Mann-Whitney U test ───────────────────────────────────────
python3 "${SCRIPTDIR}/mann_whitney_identity.py" "$OUTFILE" "$OUTDIR"

echo "  Output: $OUTFILE"
echo ""
