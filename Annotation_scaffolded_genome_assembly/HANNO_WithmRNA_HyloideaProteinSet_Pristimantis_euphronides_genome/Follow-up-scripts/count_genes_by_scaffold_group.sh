#!/bin/bash
set -euo pipefail

# =============================================================================
# count_genes_by_scaffold_group.sh
#
# Count protein-coding gene models in the HANNO BESTMODELS GFF3 per sequence
# category. Categories are deterministic from RagTag naming conventions:
#
#   13 longest scaffolds   sequence ID scaffold_1 ... scaffold_13
#                          (chromosome-scale scaffolds, one per E. coqui
#                          hap1 reference chromosome)
#
#   Scaffolds 14-31        sequence ID scaffold_14 ... scaffold_31
#                          (shorter RagTag scaffolds built against
#                          non-chromosomal reference sequences)
#
#   Unplaced contigs       sequence ID contig_*
#                          (Flye contigs not placed by RagTag)
#
# Input:  HANNO BESTMODELS GFF3, for example the aPriEup1.0_genomic.gff3
#         produced by run_annotation_pipeline.sh.
#
# Output: Two files in OUTPUT_DIR:
#           gene_counts_by_scaffold_group.tsv   Counts and percentages.
#           count_genes_by_scaffold_group.out   Reproducible run record.
#
# Usage:
#   bash count_genes_by_scaffold_group.sh \
#       -g /path/to/aPriEup1.0_genomic.gff3 \
#       -o /path/to/annotation_output
#
# Author: K. Kopp, P. euphronides genome project
# =============================================================================

usage() {
    echo "Usage: $0 -g GFF3 -o OUTPUT_DIR"
    echo ""
    echo "Required:"
    echo "  -g  HANNO BESTMODELS GFF3"
    echo "  -o  Output directory"
    exit 1
}

# --- Parse arguments ---
GFF3=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -g) GFF3="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$GFF3" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: -g and -o are required."
    usage
fi

if [[ ! -f "$GFF3" ]]; then
    echo "ERROR: GFF3 not found: $GFF3"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

OUT_TSV="$OUTPUT_DIR/gene_counts_by_scaffold_group.tsv"
OUT_RUN="$OUTPUT_DIR/count_genes_by_scaffold_group.out"

# =============================================================================
# All detailed output is written to OUT_RUN as a reproducible record.
# Only a brief summary is printed to the terminal at the end.
# =============================================================================
{
    echo "============================================================"
    echo "Gene model counts by scaffold group"
    echo "============================================================"
    echo "Date:         $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "GFF3:         $GFF3"
    echo "Output dir:   $OUTPUT_DIR"
    echo "============================================================"
    echo ""

    echo ">>> Counting 'gene' features per category..."

    # Single awk pass: classify each 'gene' feature by sequence ID pattern.
    # Categories are mutually exclusive and deterministic from RagTag naming.
    read -r N_13 N_14_31 N_CONTIG N_OTHER N_TOTAL < <(
        awk -F'\t' '
            /^#/ { next }
            $3 != "gene" { next }
            {
                total++
                seqid = $1
                if (seqid ~ /^scaffold_([1-9]|1[0-3])$/) {
                    n13++
                } else if (seqid ~ /^scaffold_(1[4-9]|2[0-9]|3[01])$/) {
                    n14_31++
                } else if (seqid ~ /^contig_/) {
                    ncontig++
                } else {
                    nother++
                }
            }
            END { printf "%d %d %d %d %d\n", n13+0, n14_31+0, ncontig+0, nother+0, total+0 }
        ' "$GFF3"
    )

    if [[ "$N_OTHER" -gt 0 ]]; then
        echo "ERROR: $N_OTHER gene features on sequence IDs that match none of"
        echo "  scaffold_1..scaffold_31 or contig_*. Check GFF3 sequence IDs."
        exit 1
    fi

    PCT_13=$(awk "BEGIN { printf \"%.2f\", ($N_13 / $N_TOTAL) * 100 }")
    PCT_14_31=$(awk "BEGIN { printf \"%.2f\", ($N_14_31 / $N_TOTAL) * 100 }")
    PCT_CONTIG=$(awk "BEGIN { printf \"%.2f\", ($N_CONTIG / $N_TOTAL) * 100 }")

    echo "  13 longest scaffolds (scaffold_1..13):   $N_13 gene models"
    echo "  Scaffolds 14-31:                         $N_14_31 gene models"
    echo "  Unplaced contigs (contig_*):             $N_CONTIG gene models"
    echo "  Total:                                   $N_TOTAL gene models"
    echo ""

    echo ">>> Writing output TSV..."
    {
        printf "scaffold_group\tn_gene_models\tpercent_of_total\n"
        printf "13_longest_scaffolds\t%s\t%s\n" "$N_13" "$PCT_13"
        printf "scaffolds_14_31\t%s\t%s\n" "$N_14_31" "$PCT_14_31"
        printf "unplaced_contigs\t%s\t%s\n" "$N_CONTIG" "$PCT_CONTIG"
        printf "total\t%s\t100.00\n" "$N_TOTAL"
    } > "$OUT_TSV"
    echo "  Written: $OUT_TSV"
    echo ""

    echo "============================================================"
    echo "Results"
    echo "============================================================"
    echo "  13 longest scaffolds: $N_13 ($PCT_13%)"
    echo "  Scaffolds 14-31:      $N_14_31 ($PCT_14_31%)"
    echo "  Unplaced contigs:     $N_CONTIG ($PCT_CONTIG%)"
    echo "  Total gene models:    $N_TOTAL"
    echo "============================================================"
} > "$OUT_RUN" 2>&1

# --- Terminal summary (variables persist from the grouped block above) ---
echo "Done."
echo "  13 longest scaffolds: $N_13 ($PCT_13%)"
echo "  Scaffolds 14-31:      $N_14_31 ($PCT_14_31%)"
echo "  Unplaced contigs:     $N_CONTIG ($PCT_CONTIG%)"
echo "  Total gene models:    $N_TOTAL"
echo ""
echo "Outputs:"
echo "  $OUT_TSV"
echo "  $OUT_RUN"
