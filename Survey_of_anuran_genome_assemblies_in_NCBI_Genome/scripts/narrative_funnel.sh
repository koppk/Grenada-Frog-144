#!/bin/bash
set -euo pipefail

# =============================================================================
# narrative_funnel.sh
#
# Builds the narrative funnel for the Discussion from the merged assembly,
# taxonomy, tech classification, and tissue classification data.
# All output saved to files for reproducibility.
#
# Prerequisites: output from analyze_anuran_genomes.sh and
#                fetch_biosample_tissue.sh in the same output directory.
#
# Usage:
#   bash narrative_funnel.sh -o OUTPUT_DIR
#
# Author: Kopp K, Pristimantis euphronides genome project
# =============================================================================

usage() {
    echo "Usage: $0 -o OUTPUT_DIR"
    echo "  OUTPUT_DIR must contain merged.tsv and assemblies_with_tissue.tsv"
    exit 1
}

OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then usage; fi

MERGED="$OUTPUT_DIR/merged.tsv"
TISSUE="$OUTPUT_DIR/assemblies_with_tissue.tsv"

for f in "$MERGED" "$TISSUE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: $f not found. Run analyze_anuran_genomes.sh and fetch_biosample_tissue.sh first."
        exit 1
    fi
done

REPORT="$OUTPUT_DIR/narrative_funnel.out"

{
echo "============================================================"
echo "NARRATIVE FUNNEL"
echo "============================================================"
echo "Date:    $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Input:   $MERGED"
echo "         $TISSUE"
echo "============================================================"
echo ""

# -------------------------------------------------------------------
# 1. Total assemblies and unique species
# -------------------------------------------------------------------
N_ASM=$(tail -n +2 "$MERGED" | wc -l)
N_SPECIES=$(tail -n +2 "$MERGED" | cut -d$'\t' -f2 | sort -u | wc -l)
echo ">>> 1. Overview"
echo "  Total assemblies: $N_ASM"
echo "  Unique species:   $N_SPECIES"
echo ""

# -------------------------------------------------------------------
# 2. Assembly level breakdown
# -------------------------------------------------------------------
echo ">>> 2. Assembly level"
tail -n +2 "$MERGED" | cut -d$'\t' -f4 | sort | uniq -c | sort -rn
echo ""

N_SCAFF_PLUS=$(tail -n +2 "$MERGED" | awk -F'\t' '$4 == "Scaffold" || $4 == "Chromosome"' | wc -l)
N_CHROM=$(tail -n +2 "$MERGED" | awk -F'\t' '$4 == "Chromosome"' | wc -l)
echo "  Scaffold or Chromosome level: $N_SCAFF_PLUS"
echo "  Chromosome level:             $N_CHROM"
echo ""

# -------------------------------------------------------------------
# 3. Hyloidea
# -------------------------------------------------------------------
echo ">>> 3. Hyloidea assemblies"
N_HYL=$(tail -n +2 "$MERGED" | awk -F'\t' '$9 == "Hyloidea"' | wc -l)
N_HYL_SP=$(tail -n +2 "$MERGED" | awk -F'\t' '$9 == "Hyloidea" {print $2}' | sort -u | wc -l)
echo "  Assemblies: $N_HYL"
echo "  Species:    $N_HYL_SP"
echo "  By level:"
tail -n +2 "$MERGED" | awk -F'\t' '$9 == "Hyloidea" {print $4}' | sort | uniq -c | sort -rn
echo ""

# -------------------------------------------------------------------
# 4. Terrarana — families per Hedges et al. 2008 [19] + Heinicke et al. 2009
# -------------------------------------------------------------------
# Hedges et al. 2008 erected Terrarana with four families:
#   Brachycephalidae, Craugastoridae, Eleutherodactylidae, Strabomantidae
# Heinicke et al. 2009 added Ceuthomantidae.
# Some treatments merge Strabomantidae into Craugastoridae.
# All five included here to cover both frameworks.
TERRARANA_FAMILIES="Brachycephalidae|Craugastoridae|Strabomantidae|Eleutherodactylidae|Ceuthomantidae"

echo ">>> 4. Terrarana assemblies (families: $TERRARANA_FAMILIES)"
echo ""
echo "  All Terrarana assemblies:"
echo -e "  accession\torganism\ttech_class\tlevel\tfamily"
tail -n +2 "$MERGED" | awk -F'\t' -v fams="$TERRARANA_FAMILIES" '$8 ~ fams {
    printf "  %s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $8
}' | sort -t$'\t' -k5
N_TERR=$(tail -n +2 "$MERGED" | awk -F'\t' -v fams="$TERRARANA_FAMILIES" '$8 ~ fams' | wc -l)
N_TERR_SP=$(tail -n +2 "$MERGED" | awk -F'\t' -v fams="$TERRARANA_FAMILIES" '$8 ~ fams {print $2}' | sort -u | wc -l)
echo ""
echo "  Total Terrarana assemblies: $N_TERR"
echo "  Total Terrarana species:    $N_TERR_SP"
echo ""

# -------------------------------------------------------------------
# 5. Pristimantis genus
# -------------------------------------------------------------------
echo ">>> 5. Pristimantis genus"
N_PRIST=$(tail -n +2 "$MERGED" | awk -F'\t' '$2 ~ /^Pristimantis/' | wc -l)
echo "  Assemblies: $N_PRIST"
tail -n +2 "$MERGED" | awk -F'\t' '$2 ~ /^Pristimantis/ {
    printf "  %s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $8
}'
echo ""

# -------------------------------------------------------------------
# 6. Strabomantidae / Craugastoridae (excluding Eleutherodactylidae)
# -------------------------------------------------------------------
echo ">>> 6. Strabomantidae or Craugastoridae only"
N_STRAB=$(tail -n +2 "$MERGED" | awk -F'\t' '$8 ~ /Strabomantidae|Craugastoridae/' | wc -l)
echo "  Assemblies: $N_STRAB"
tail -n +2 "$MERGED" | awk -F'\t' '$8 ~ /Strabomantidae|Craugastoridae/ {
    printf "  %s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $8
}'
echo ""

# -------------------------------------------------------------------
# 7. Tech class breakdown
# -------------------------------------------------------------------
echo ">>> 7. Sequencing technology classification"
echo "  All assemblies:"
tail -n +2 "$MERGED" | cut -d$'\t' -f3 | sort | uniq -c | sort -rn
echo ""
echo "  Scaffold+ level only:"
tail -n +2 "$MERGED" | awk -F'\t' '$4 == "Scaffold" || $4 == "Chromosome" {print $3}' | sort | uniq -c | sort -rn
echo ""

# -------------------------------------------------------------------
# 8. ONT-only assemblies
# -------------------------------------------------------------------
echo ">>> 8. ONT-only assemblies"
echo -e "  accession\torganism\tlevel\tscaffold_N50\ttissue\ttissue_class"
tail -n +2 "$TISSUE" | awk -F'\t' '$3 ~ /^Oxford Nanopore/ || $3 == "ONT" {
    if ($3 !~ /Illumina/ && $3 !~ /PacBio/ && $3 !~ /10x/ && $3 !~ /Chromium/ && $3 !~ /Hi-C/ && $3 !~ /HiC/ && $3 !~ /Arima/)
        printf "  %s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $4, $6, $9, $10
}'
echo ""

# -------------------------------------------------------------------
# 9. Tissue classification (Taberlet 1999)
# -------------------------------------------------------------------
echo ">>> 9. Tissue classification (Taberlet 1999)"
tail -n +2 "$TISSUE" | cut -d$'\t' -f10 | sort | uniq -c | sort -rn
echo ""

# -------------------------------------------------------------------
# 10. Not confirmed destructive: nondestructive + ambiguous at scaffold+ level
# -------------------------------------------------------------------
echo ">>> 10. Nondestructive or ambiguous assemblies at scaffold or chromosome level"
echo "  (animal may have survived — cannot be determined from metadata alone)"
echo ""
echo -e "  accession\torganism\ttech\tlevel\ttissue\tclass"
tail -n +2 "$TISSUE" | awk -F'\t' '($10 == "nondestructive" || $10 == "ambiguous") && ($4 == "Scaffold" || $4 == "Chromosome") {
    printf "  %s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $9, $10
}'
N_ND_SCAFF=$(tail -n +2 "$TISSUE" | awk -F'\t' '$10 == "nondestructive" && ($4 == "Scaffold" || $4 == "Chromosome")' | wc -l)
N_AMB_SCAFF=$(tail -n +2 "$TISSUE" | awk -F'\t' '$10 == "ambiguous" && ($4 == "Scaffold" || $4 == "Chromosome")' | wc -l)
echo ""
echo "  Nondestructive: $N_ND_SCAFF"
echo "  Ambiguous:      $N_AMB_SCAFF"
echo ""

# -------------------------------------------------------------------
# 11. Not confirmed destructive AND ONT-only
# -------------------------------------------------------------------
echo ">>> 11. Nondestructive or ambiguous AND ONT-only"
tail -n +2 "$TISSUE" | awk -F'\t' '
    ($10 == "nondestructive" || $10 == "ambiguous") {
        t = tolower($3)
        gsub(/oxford/, "", t); gsub(/nanopore/, "", t); gsub(/technologies/, "", t)
        gsub(/minion/, "", t); gsub(/promethion/, "", t); gsub(/gridion/, "", t)
        gsub(/flongle/, "", t); gsub(/sequencing/, "", t)
        gsub(/[;,\/\-\.\(\) \t]/, "", t)
        if (t == "" || t ~ /^ont$/)
            printf "  %s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $9, $10
    }'
echo ""

# -------------------------------------------------------------------
# 12. Scaffold+ level by family (breadth of representation)
# -------------------------------------------------------------------
echo ">>> 12. Families with scaffold+ assemblies"
echo -e "  family\tcontig\tscaffold\tchromosome\ttotal"
tail -n +2 "$MERGED" | awk -F'\t' '{
    fam=$8; lvl=$4
    total[fam]++
    if (lvl == "Contig") contig[fam]++
    else if (lvl == "Scaffold") scaffold[fam]++
    else if (lvl == "Chromosome") chromosome[fam]++
} END {
    for (f in total) printf "  %s\t%d\t%d\t%d\t%d\n", f, contig[f]+0, scaffold[f]+0, chromosome[f]+0, total[f]
}' | sort -t$'\t' -k5 -rn
echo ""

# -------------------------------------------------------------------
# 13. Funnel summary
# -------------------------------------------------------------------
echo "============================================================"
echo "FUNNEL SUMMARY"
echo "============================================================"
echo ""
echo "  All anuran genome assemblies in NCBI:           $N_ASM"
echo "  Unique species:                                 $N_SPECIES"
echo "  Scaffold or chromosome level:                   $N_SCAFF_PLUS"
echo "  Chromosome level:                               $N_CHROM"
echo "  Hyloidea assemblies:                            $N_HYL ($N_HYL_SP species)"
echo "  Terrarana assemblies:                            $N_TERR ($N_TERR_SP species)"
echo "  Strabomantidae/Craugastoridae assemblies:        $N_STRAB"
echo "  Pristimantis assemblies:                         $N_PRIST"
echo "  Nondestructive at scaffold+ level:               $N_ND_SCAFF"
echo "  Ambiguous at scaffold+ level:                    $N_AMB_SCAFF"
echo "  Nondestructive or ambiguous + ONT-only:          (see section 11)"
echo ""
echo "============================================================"

} > "$REPORT" 2>&1

echo "Done."
echo "Full report: $REPORT"
