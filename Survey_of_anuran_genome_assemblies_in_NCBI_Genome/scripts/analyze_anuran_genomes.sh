#!/bin/bash
set -euo pipefail

# =============================================================================
# analyze_anuran_genomes.sh
#
# Granular analysis of NCBI Anura genome assemblies for narrative claims.
# Classifies by: sequencing technology, taxonomy (family/superfamily),
# assembly level, and cross-tabulations.
#
# Prerequisites:
#   conda install -c conda-forge ncbi-datasets-cli
#   conda install -c bioconda taxonkit
#   taxonkit needs NCBI taxonomy dump:
#     wget -qO- ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz | tar xz -C ~/.taxonkit/
#
# Usage:
#   bash analyze_anuran_genomes.sh -i anuran_genomes_ncbi.tsv -o output_dir
#
# Author: Kopp K, Pristimantis euphronides genome project
# =============================================================================

usage() {
    echo "Usage: $0 -i INPUT_TSV -o OUTPUT_DIR"
    echo ""
    echo "Required:"
    echo "  -i  Input TSV from scrub_anuran_genomes_ncbi.sh"
    echo "  -o  Output directory"
    exit 1
}

INPUT=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: -i and -o are required."
    usage
fi

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

for cmd in taxonkit; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: '$cmd' not found. Install with: conda install -c bioconda taxonkit"
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

REPORT="$OUTPUT_DIR/anuran_genome_analysis.out"

{
echo "============================================================"
echo "Granular analysis of NCBI Anura genome assemblies"
echo "============================================================"
echo "Date:    $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Input:   $INPUT"
echo "Output:  $OUTPUT_DIR"
echo "============================================================"
echo ""

# -------------------------------------------------------------------
# 1. Deduplicate GCA/GCF pairs
# -------------------------------------------------------------------
# GCA (GenBank) and GCF (RefSeq) accessions with the same assembly name
# represent the same underlying assembly. Keep GCF when both exist.
echo ">>> Step 1: Deduplicating GCA/GCF pairs by assembly name..."
DEDUP="$OUTPUT_DIR/dedup.tsv"
head -1 "$INPUT" > "$DEDUP"
tail -n +2 "$INPUT" | awk -F'\t' '
BEGIN { OFS="\t" }
{
    acc = $1
    asm_name = $11   # assembly name is column 11

    # Use assembly name as dedup key if available, else numeric accession part
    if (asm_name != "" && asm_name != "na") {
        key = asm_name
    } else {
        key = acc
        sub(/^GC[AF]_/, "", key)
    }

    if (!(key in seen)) {
        seen[key] = $0
        prefix[key] = substr(acc, 1, 3)
    } else if (substr(acc, 1, 3) == "GCF" && prefix[key] != "GCF") {
        # GCF (RefSeq) takes priority over GCA (GenBank)
        seen[key] = $0
        prefix[key] = "GCF"
    }
}
END {
    for (key in seen) print seen[key]
}' | sort -t$'\t' -k1,1 >> "$DEDUP"

N_RAW=$(tail -n +2 "$INPUT" | wc -l)
N_DEDUP=$(tail -n +2 "$DEDUP" | wc -l)
N_REMOVED=$((N_RAW - N_DEDUP))
echo "  Raw assemblies:    $N_RAW"
echo "  After dedup:       $N_DEDUP"
echo "  GCA/GCF pairs removed: $N_REMOVED"
echo ""

# -------------------------------------------------------------------
# 2. Extract unique species and look up taxonomy
# -------------------------------------------------------------------
echo ">>> Step 2: Looking up taxonomy via taxonkit..."
SPECIES_LIST="$OUTPUT_DIR/species_list.txt"
TAXONOMY="$OUTPUT_DIR/taxonomy.tsv"

tail -n +2 "$DEDUP" | cut -d$'\t' -f2 | sort -u > "$SPECIES_LIST"
N_SPECIES=$(wc -l < "$SPECIES_LIST")
echo "  Unique species: $N_SPECIES"

# name2taxid then lineage
taxonkit name2taxid "$SPECIES_LIST" 2>/dev/null \
    | taxonkit lineage -i 2 -r -n 2>/dev/null \
    | awk -F'\t' 'BEGIN{OFS="\t"; print "species","taxid","lineage","taxid_lineage","rank"} {print}' \
    > "$TAXONOMY"

# Extract family and superfamily from lineage
TAXMAP="$OUTPUT_DIR/taxmap.tsv"
echo -e "species\tfamily\tsuperfamily\torder" > "$TAXMAP"
tail -n +2 "$TAXONOMY" | awk -F'\t' '
BEGIN { OFS="\t" }
{
    species = $1
    lineage = $3
    family = "unknown"
    superfamily = "unknown"

    n = split(lineage, parts, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", parts[i])
        if (parts[i] ~ /idae$/ && family == "unknown") family = parts[i]
        if (parts[i] ~ /oidea$/ && superfamily == "unknown") superfamily = parts[i]
    }
    print species, family, superfamily, "Anura"
}' >> "$TAXMAP"

echo "  Taxonomy mapped: $TAXMAP"
echo ""

# -------------------------------------------------------------------
# 3. Classify sequencing technology
# -------------------------------------------------------------------
echo ">>> Step 3: Classifying sequencing technology..."

TECHCLASS="$OUTPUT_DIR/tech_classified.tsv"

tail -n +2 "$DEDUP" | awk -F'\t' '
BEGIN {
    OFS="\t"
    print "accession","organism","tech_raw","tech_class","assembly_level","total_length","scaffold_n50","contig_n50"
}
{
    acc=$1; org=$2; level=$4; tech=$5; size=$7; sn50=$8; cn50=$9
    t = tolower(tech)

    has_ont=0; has_illumina=0; has_pacbio=0; has_hic=0
    has_10x=0; has_bionano=0; has_other=0

    if (t ~ /nanopore|minion|gridion|promethion|flongle/ || t ~ /(^|[,; ])ont([,; ]|$)/) has_ont=1
    if (t ~ /illumina|hiseq|novaseq|nextseq|miseq/) has_illumina=1
    if (t ~ /pacbio|sequel|revio|rsii|hifi|clr|smrt/) has_pacbio=1
    if (t ~ /hi-c|hic|arima|dovetail|omni/) has_hic=1
    if (t ~ /10x|chromium|linked.read/) has_10x=1
    if (t ~ /bionano|optical|saphyr/) has_bionano=1
    if (t ~ /stlfr|sanger|454|solid|ion.torrent/) has_other=1

    if (has_ont && !has_illumina && !has_pacbio && !has_hic && !has_10x && !has_bionano && !has_other)
        tc = "ONT_only"
    else if (has_ont)
        tc = "ONT_plus_other"
    else if (has_pacbio && has_hic)
        tc = "PacBio_HiC"
    else if (has_pacbio && has_illumina)
        tc = "PacBio_Illumina"
    else if (has_pacbio)
        tc = "PacBio_only"
    else if (has_illumina)
        tc = "Illumina_only"
    else
        tc = "Other"

    print acc, org, tech, tc, level, size, sn50, cn50
}' > "$TECHCLASS"

echo "  Tech classified: $TECHCLASS"
echo ""

# -------------------------------------------------------------------
# 4. Merge taxonomy with tech classification
# -------------------------------------------------------------------
echo ">>> Step 4: Merging taxonomy and tech classification..."

MERGED="$OUTPUT_DIR/merged.tsv"
# Use awk to join tech classification with taxonomy
awk -F'\t' '
BEGIN { OFS="\t"; header_done=0 }
NR==FNR && FNR>1 {
    # Read taxmap: species -> family, superfamily
    fam[$1] = $2
    sfam[$1] = $3
    next
}
NR==FNR { next }
FNR==1 {
    if (!header_done) {
        print "accession","organism","tech_class","assembly_level","total_length","scaffold_n50","contig_n50","family","superfamily"
        header_done=1
    }
    next
}
{
    org = $2
    f = (org in fam) ? fam[org] : "unknown"
    s = (org in sfam) ? sfam[org] : "unknown"
    print $1, $2, $4, $5, $6, $7, $8, f, s
}' "$TAXMAP" "$TECHCLASS" > "$MERGED"
echo "  Merged: $MERGED"
echo ""

# -------------------------------------------------------------------
# 5. Summary tables
# -------------------------------------------------------------------
echo "============================================================"
echo "SUMMARY TABLES"
echo "============================================================"
echo ""

echo "--- A. All assemblies by technology class ---"
tail -n +2 "$MERGED" | cut -d$'\t' -f3 | sort | uniq -c | sort -rn
echo ""

echo "--- B. All assemblies by assembly level ---"
tail -n +2 "$MERGED" | cut -d$'\t' -f4 | sort | uniq -c | sort -rn
echo ""

echo "--- C. All assemblies by family ---"
tail -n +2 "$MERGED" | cut -d$'\t' -f8 | sort | uniq -c | sort -rn
echo ""

echo "--- D. All assemblies by superfamily ---"
tail -n +2 "$MERGED" | cut -d$'\t' -f9 | sort | uniq -c | sort -rn
echo ""

echo "--- E. ONT-only assemblies (all details) ---"
echo -e "accession\torganism\tlevel\tsize_Mb\tscaffold_N50_Mb\tcontig_N50_Mb\tfamily\tsuperfamily"
tail -n +2 "$MERGED" | awk -F'\t' '$3 == "ONT_only" {printf "%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%s\t%s\n", $1, $2, $4, $5/1e6, $6/1e6, $7/1e6, $8, $9}'
echo ""

echo "--- F. ONT-containing assemblies (ONT + any other tech) ---"
echo -e "accession\torganism\ttech_class\tlevel\tfamily"
tail -n +2 "$MERGED" | awk -F'\t' '$3 ~ /ONT/ {print $1"\t"$2"\t"$3"\t"$4"\t"$8}'
echo ""

echo "--- G. Hyloidea assemblies ---"
echo -e "accession\torganism\ttech_class\tlevel\tfamily"
tail -n +2 "$MERGED" | awk -F'\t' '$9 == "Hyloidea" {print $1"\t"$2"\t"$3"\t"$4"\t"$8}'
echo ""

echo "--- H. Terrarana assemblies (Strabomantidae, Craugastoridae, Eleutherodactylidae, Brachycephalidae) ---"
echo -e "accession\torganism\ttech_class\tlevel\tfamily"
tail -n +2 "$MERGED" | awk -F'\t' '$8 == "Strabomantidae" || $8 == "Craugastoridae" || $8 == "Eleutherodactylidae" || $8 == "Brachycephalidae" {print $1"\t"$2"\t"$3"\t"$4"\t"$8}'
echo ""

echo "--- I. Strabomantidae / Craugastoridae assemblies ---"
echo -e "accession\torganism\ttech_class\tlevel\tfamily"
tail -n +2 "$MERGED" | awk -F'\t' '$8 == "Strabomantidae" || $8 == "Craugastoridae" {print $1"\t"$2"\t"$3"\t"$4"\t"$8}'
echo ""

echo "--- J. Chromosome-level assemblies by technology ---"
echo -e "tech_class\tcount"
tail -n +2 "$MERGED" | awk -F'\t' '$4 == "Chromosome" {print $3}' | sort | uniq -c | sort -rn
echo ""

echo "--- K. Cross-table: assembly level × technology class ---"
echo ""
tail -n +2 "$MERGED" | awk -F'\t' '{key=$4"\t"$3; count[key]++} END {for (k in count) print count[k]"\t"k}' | sort -rn
echo ""

echo "--- L. Scaffold-level or better assemblies by family (breadth of representation) ---"
echo -e "family\tcontig\tscaffold\tchromosome\ttotal"
tail -n +2 "$MERGED" | awk -F'\t' '{
    fam=$8; lvl=$4
    total[fam]++
    if (lvl == "Contig") contig[fam]++
    else if (lvl == "Scaffold") scaffold[fam]++
    else if (lvl == "Chromosome") chromosome[fam]++
} END {
    for (f in total) printf "%s\t%d\t%d\t%d\t%d\n", f, contig[f]+0, scaffold[f]+0, chromosome[f]+0, total[f]
}' | sort -t$'\t' -k5 -rn
echo ""

echo "============================================================"
echo "NARRATIVE CLAIM CHECKLIST"
echo "============================================================"
echo ""

n_ont_only=$(tail -n +2 "$MERGED" | awk -F'\t' '$3 == "ONT_only"' | wc -l)
n_ont_only_cr_note="(IUCN status not in NCBI — check manually)"
n_strab=$(tail -n +2 "$MERGED" | awk -F'\t' '$8 == "Strabomantidae" || $8 == "Craugastoridae"' | wc -l)
n_terrarana=$(tail -n +2 "$MERGED" | awk -F'\t' '$8 == "Strabomantidae" || $8 == "Craugastoridae" || $8 == "Eleutherodactylidae" || $8 == "Brachycephalidae"' | wc -l)
n_hyloidea=$(tail -n +2 "$MERGED" | awk -F'\t' '$9 == "Hyloidea"' | wc -l)
n_ont_scaffolded=$(tail -n +2 "$MERGED" | awk -F'\t' '$3 == "ONT_only" && ($4 == "Scaffold" || $4 == "Chromosome")' | wc -l)
n_ont_scaffold_n50_gt10=$(tail -n +2 "$MERGED" | awk -F'\t' '$3 == "ONT_only" && $6 > 10000000' | wc -l)

echo "1. ONT-only anuran genomes in NCBI:       $n_ont_only"
echo "2. ONT-only with scaffold N50 > 10 Mb:    $n_ont_scaffold_n50_gt10"
echo "3. Strabomantidae/Craugastoridae genomes:  $n_strab"
echo "4. Terrarana genomes (all tech):           $n_terrarana"
echo "5. Hyloidea genomes (all tech):            $n_hyloidea"
echo ""
echo "IUCN CR status: $n_ont_only_cr_note"
echo "Non-lethal sampling: not recorded in NCBI metadata — check publications"
echo ""
echo "============================================================"

} > "$REPORT" 2>&1

# Terminal summary
echo "Done."
echo "Full report: $REPORT"
echo "Merged data: $MERGED"
