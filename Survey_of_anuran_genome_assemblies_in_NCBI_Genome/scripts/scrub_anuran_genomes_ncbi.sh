#!/bin/bash
set -euo pipefail

# =============================================================================
# scrub_anuran_genomes_ncbi.sh
#
# Query NCBI for all Anura genome assemblies and extract sequencing technology,
# assembly level, and genome scope (nuclear vs organellar). Produces a TSV that
# can be filtered for ONT-only nuclear genome assemblies and cross-referenced
# with IUCN Red List status.
#
# Prerequisites:
#   conda install -c conda-forge ncbi-datasets-cli
#
# Usage:
#   bash scrub_anuran_genomes_ncbi.sh -o /path/to/output_dir
#
# Output:
#   anuran_genomes_ncbi.tsv          Full metadata for all Anura assemblies
#   anuran_genomes_ont_only.tsv      Subset: assemblies where sequencing tech
#                                    contains "Oxford Nanopore" and does NOT
#                                    contain "Illumina", "PacBio", or "HiFi"
#   scrub_anuran_genomes_ncbi.out    Run record
#
# Author: K. Kopp, P. euphronides genome project
# =============================================================================

usage() {
    echo "Usage: $0 -o OUTPUT_DIR"
    echo ""
    echo "Required:"
    echo "  -o  Output directory"
    echo ""
    echo "Prerequisites:"
    echo "  conda install -c conda-forge ncbi-datasets-cli"
    exit 1
}

# --- Parse arguments ---
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: -o is required."
    usage
fi

# --- Check prerequisites ---
if ! command -v datasets &>/dev/null; then
    echo "ERROR: 'datasets' command not found."
    echo "Install with: conda install -c conda-forge ncbi-datasets-cli"
    exit 1
fi

if ! command -v dataformat &>/dev/null; then
    echo "ERROR: 'dataformat' command not found."
    echo "Install with: conda install -c conda-forge ncbi-datasets-cli"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

OUT_FULL="$OUTPUT_DIR/anuran_genomes_ncbi.tsv"
OUT_ONT="$OUTPUT_DIR/anuran_genomes_ont_only.tsv"
OUT_RUN="$OUTPUT_DIR/scrub_anuran_genomes_ncbi.out"

{
    echo "============================================================"
    echo "NCBI Anura genome assembly metadata scrub"
    echo "============================================================"
    echo "Date:         $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Output dir:   $OUTPUT_DIR"
    echo "datasets ver: $(datasets version 2>&1 | head -1)"
    echo "============================================================"
    echo ""

    # Anura NCBI Taxonomy ID: 8342
    echo ">>> Querying NCBI for all Anura (taxid 8342) genome assemblies..."

    # Show available fields for diagnostics
    echo "  Available dataformat fields:"
    dataformat tsv genome --help 2>&1 | grep -i "fields\|available" | head -3
    echo ""

    FIELDS="accession,organism-name,organism-infraspecific-strain,assminfo-level,assminfo-sequencing-tech,assminfo-biosample-accession,assmstats-total-sequence-len,assmstats-scaffold-n50,assmstats-contig-n50,assminfo-submitter,assminfo-name"

    if ! datasets summary genome taxon 8342 --as-json-lines 2>/tmp/datasets_err.txt \
        | dataformat tsv genome --fields "$FIELDS" \
        > "$OUT_FULL" 2>/tmp/dataformat_err.txt; then
        echo "  ERROR: dataformat failed. Trying alternative field names..."
        cat /tmp/dataformat_err.txt
        # Fallback: try without organism-infraspecific-strain
        FIELDS="accession,organism-name,assminfo-level,assminfo-sequencing-tech,assmstats-total-sequence-len,assmstats-scaffold-n50,assmstats-contig-n50"
        datasets summary genome taxon 8342 --as-json-lines 2>/dev/null \
            | dataformat tsv genome --fields "$FIELDS" \
            > "$OUT_FULL" 2>&1
    fi

    N_TOTAL=$(tail -n +2 "$OUT_FULL" | wc -l)
    if [[ "$N_TOTAL" -eq 0 ]]; then
        echo "  WARNING: No assemblies retrieved. Check field names."
        echo "  TSV header: $(head -1 "$OUT_FULL")"
        echo "  datasets stderr: $(cat /tmp/datasets_err.txt 2>/dev/null)"
        echo "  dataformat stderr: $(cat /tmp/dataformat_err.txt 2>/dev/null)"
    fi
    echo "  Total Anura assemblies found: $N_TOTAL"
    echo "  Written: $OUT_FULL"
    echo ""

    echo ">>> Filtering for ONT-only nuclear genome assemblies..."

    # Header
    head -1 "$OUT_FULL" > "$OUT_ONT"

    # Filter: whitelist approach — strip all known ONT-related terms from the
    # tech field. If nothing substantive remains, the assembly is ONT-only.
    # This avoids maintaining a blacklist of every possible supplementary tech.
    tail -n +2 "$OUT_FULL" \
        | awk -F'\t' -v tech_col="$(head -1 "$OUT_FULL" | tr '\t' '\n' | grep -n -i "sequencing" | head -1 | cut -d: -f1)" '
            BEGIN { if (tech_col == "") tech_col = 5 }
            {
                tech = tolower($tech_col)
                # Must contain at least one ONT indicator
                if (!(tech ~ /nanopore/ || tech ~ /minion/ || tech ~ /gridion/ || tech ~ /promethion/ || tech ~ /flongle/ || tech ~ /(^|[,; ])ont([,; ]|$)/)) next

                # Strip all ONT-related terms
                stripped = tech
                gsub(/oxford/, "", stripped)
                gsub(/nanopore/, "", stripped)
                gsub(/technologies/, "", stripped)
                gsub(/ont/, "", stripped)
                gsub(/minion/, "", stripped)
                gsub(/gridion/, "", stripped)
                gsub(/promethion/, "", stripped)
                gsub(/flongle/, "", stripped)
                gsub(/mk1[bc]/, "", stripped)
                gsub(/p2[4]?/, "", stripped)
                gsub(/p48/, "", stripped)
                gsub(/sequencing/, "", stripped)

                # Strip separators and whitespace
                gsub(/[;,\/\-\.\(\)]/, "", stripped)
                gsub(/[ \t]+/, "", stripped)

                # If nothing remains, it is ONT-only
                if (stripped == "") print
            }
        ' >> "$OUT_ONT"

    N_ONT=$(tail -n +2 "$OUT_ONT" | wc -l)
    echo "  ONT-only assemblies found: $N_ONT"
    echo "  Written: $OUT_ONT"
    echo ""

    if [[ "$N_ONT" -gt 0 ]]; then
        echo ">>> ONT-only Anura assemblies:"
        echo ""
        tail -n +2 "$OUT_ONT" | while IFS=$'\t' read -r acc org strain level tech biosample size sn50 cn50 submitter asm_name; do
            echo "  $acc  $org"
            echo "    Level: $level"
            echo "    Tech:  $tech"
            echo "    Size:  $size"
            echo ""
        done
    else
        echo "  No ONT-only Anura genome assemblies found in NCBI."
    fi

    echo "============================================================"
    echo "Next steps"
    echo "============================================================"
    echo ""
    echo "1. Check IUCN Red List status for any species in $OUT_ONT"
    echo "   https://www.iucnredlist.org/"
    echo ""
    echo "2. Verify whether assemblies are nuclear (check total size;"
    echo "   mitogenomes are ~15-20 kb) vs whole-genome."
    echo ""
    echo "3. Check whether any assembly listed as ONT-only also used"
    echo "   supplementary short reads not captured in NCBI metadata"
    echo "   (read the original paper)."
    echo "============================================================"

} > "$OUT_RUN" 2>&1

# --- Terminal summary ---
N_TOTAL=$(tail -n +2 "$OUT_FULL" | wc -l)
N_ONT=$(tail -n +2 "$OUT_ONT" | wc -l)

echo "Done."
echo "  Total Anura assemblies: $N_TOTAL"
echo "  ONT-only assemblies:    $N_ONT"
echo ""
echo "Outputs:"
echo "  $OUT_FULL"
echo "  $OUT_ONT"
echo "  $OUT_RUN"
echo ""
if [[ "$N_ONT" -gt 0 ]]; then
    echo "ONT-only assemblies:"
    tail -n +2 "$OUT_ONT" | awk -F'\t' '{ printf "  %-20s %s (%.1f Gb, %s)\n", $1, $2, $7/1e9, $4 }'
fi
