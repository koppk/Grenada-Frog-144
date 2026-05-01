#!/bin/bash
set -euo pipefail

# =============================================================================
# run_annotation_pipeline.sh
#
# One-command pipeline to convert HANNO BESTMODELS output into a
# RefSeq-style annotation file set suitable for NCBI/Zenodo submission,
# with consistent duplicate removal across all output files.
#
# Input:  HANNO output directory containing:
#           BESTMODELS-FINAL.bedDB
#           BESTMODELS-FINAL.AA.faa
#           BESTMODELS-FINAL.CDS.fa
#           BESTMODELS-FINAL.mRNA.fa
#           BESTMODELS-FINAL.gtf  (optional, copied for reference)
#
# Output: RefSeq-style directory with:
#           aPriEup1.0_genomic.gff3
#           aPriEup1.0_protein.faa
#           aPriEup1.0_cds_from_genomic.fna
#           aPriEup1.0_rna_from_genomic.fna
#           annotation_summary.tsv
#           annotation_summary_duplicates.txt
#           retained_ids.txt
#           duplicate_ids.txt
#
# Usage:
#   bash run_annotation_pipeline.sh \
#       -i /path/to/HANNO_output_dir \
#       -o /path/to/output_dir \
#       [-g /path/to/genome.fasta] \
#       [-s /path/to/scripts_dir] \
#       [-p aPriEup1.0] \
#       [--genome-size 1748533034]
#
# Author: K. Kopp, P. euphronides genome project
# =============================================================================

usage() {
    echo "Usage: $0 -i HANNO_DIR -o OUTPUT_DIR [-g GENOME_FASTA] [-s SCRIPTS_DIR] [-p PREFIX] [--genome-size SIZE]"
    echo ""
    echo "Required:"
    echo "  -i  HANNO output directory (contains BESTMODELS-FINAL.bedDB etc.)"
    echo "  -o  Output directory for RefSeq-style files"
    echo ""
    echo "Optional:"
    echo "  -g  Genome FASTA (for ##sequence-region headers in GFF3)"
    echo "  -s  Directory containing the Python scripts (default: same as this script)"
    echo "  -p  File prefix (default: aPriEup1.0)"
    echo "  --genome-size  Assembly size in bp (default: 1748533034)"
    exit 1
}

# --- Defaults ---
PREFIX="aPriEup1.0"
GENOME_SIZE=1748533034
GENOME_FASTA=""
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"

# --- Parse arguments ---
HANNO_DIR=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) HANNO_DIR="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -g) GENOME_FASTA="$2"; shift 2 ;;
        -s) SCRIPTS_DIR="$2"; shift 2 ;;
        -p) PREFIX="$2"; shift 2 ;;
        --genome-size) GENOME_SIZE="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$HANNO_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: -i and -o are required."
    usage
fi

# --- Validate inputs ---
BEDDB="$HANNO_DIR/BESTMODELS-FINAL.bedDB"
PROTEIN_FA="$HANNO_DIR/BESTMODELS-FINAL.AA.faa"
CDS_FA="$HANNO_DIR/BESTMODELS-FINAL.CDS.fa"
MRNA_FA="$HANNO_DIR/BESTMODELS-FINAL.mRNA.fa"
HANNO_GTF="$HANNO_DIR/BESTMODELS-FINAL.gtf"

for f in "$BEDDB" "$PROTEIN_FA" "$CDS_FA" "$MRNA_FA"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

for script in beddb_to_gff3.py filter_hanno_fasta.py parse_hanno_annotation_summary.py; do
    if [[ ! -f "$SCRIPTS_DIR/$script" ]]; then
        echo "ERROR: Script not found: $SCRIPTS_DIR/$script"
        exit 1
    fi
done

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "Annotation Pipeline: HANNO -> RefSeq-style files"
echo "============================================================"
echo "HANNO dir:    $HANNO_DIR"
echo "Output dir:   $OUTPUT_DIR"
echo "Prefix:       $PREFIX"
echo "Genome size:  $GENOME_SIZE bp"
echo "Scripts dir:  $SCRIPTS_DIR"
[[ -n "$GENOME_FASTA" ]] && echo "Genome FASTA: $GENOME_FASTA"
echo "============================================================"
echo ""

# =============================================
# STEP 1: bedDB -> GFF3 + ID lists
# =============================================
echo ">>> STEP 1: Converting bedDB to GFF3..."
GFF3_OUT="$OUTPUT_DIR/${PREFIX}_genomic.gff3"

GFF3_CMD="python3 $SCRIPTS_DIR/beddb_to_gff3.py $BEDDB"
[[ -n "$GENOME_FASTA" ]] && GFF3_CMD="$GFF3_CMD $GENOME_FASTA"
GFF3_CMD="$GFF3_CMD -o $GFF3_OUT --stats"

eval "$GFF3_CMD"
echo ""

# Verify ID lists were created
RETAINED_IDS="$OUTPUT_DIR/retained_ids.txt"
DUPLICATE_IDS="$OUTPUT_DIR/duplicate_ids.txt"

if [[ ! -f "$RETAINED_IDS" ]]; then
    echo "ERROR: retained_ids.txt was not created by beddb_to_gff3.py"
    exit 1
fi

N_RETAINED=$(wc -l < "$RETAINED_IDS")
N_DUPLICATES=$(grep -cv '^#' "$DUPLICATE_IDS" 2>/dev/null || echo 0)
echo "  Retained: $N_RETAINED genes"
echo "  Duplicates removed: $N_DUPLICATES"
echo ""

# =============================================
# STEP 2: Filter HANNO FASTAs
# =============================================
echo ">>> STEP 2: Filtering HANNO FASTA files..."
python3 "$SCRIPTS_DIR/filter_hanno_fasta.py" \
    --ids "$RETAINED_IDS" \
    --protein "$PROTEIN_FA" \
    --cds "$CDS_FA" \
    --mrna "$MRNA_FA" \
    -o "$OUTPUT_DIR" \
    --prefix "$PREFIX"
echo ""

# =============================================
# STEP 3: Compute annotation summary
# =============================================
echo ">>> STEP 3: Computing annotation summary statistics..."
SUMMARY_OUT="$OUTPUT_DIR/annotation_summary.tsv"
python3 "$SCRIPTS_DIR/parse_hanno_annotation_summary.py" \
    "$BEDDB" \
    --genome-size "$GENOME_SIZE" \
    -o "$SUMMARY_OUT"
echo ""

# =============================================
# STEP 4: Copy HANNO GTF for reference
# =============================================
if [[ -f "$HANNO_GTF" ]]; then
    echo ">>> STEP 4: Copying HANNO original GTF for reference..."
    cp "$HANNO_GTF" "$OUTPUT_DIR/${PREFIX}_HANNO_original.gtf"
    echo "  Copied to ${PREFIX}_HANNO_original.gtf"
    echo ""
fi

# =============================================
# STEP 5: Verify consistency
# =============================================
echo ">>> STEP 5: Verifying file consistency..."
echo ""

# Count genes in GFF3
N_GFF3_GENES=$(grep -c $'\tgene\t' "$GFF3_OUT")

# Count sequences in filtered FASTAs
N_PROT=$(grep -c '^>' "$OUTPUT_DIR/${PREFIX}_protein.faa")
N_CDS=$(grep -c '^>' "$OUTPUT_DIR/${PREFIX}_cds_from_genomic.fna")
N_MRNA=$(grep -c '^>' "$OUTPUT_DIR/${PREFIX}_rna_from_genomic.fna")

# Count from summary
N_SUMMARY=$(head -2 "$SUMMARY_OUT" | tail -1 | cut -f2)

echo "  GFF3 gene features:  $N_GFF3_GENES"
echo "  Protein sequences:   $N_PROT"
echo "  CDS sequences:       $N_CDS"
echo "  mRNA sequences:      $N_MRNA"
echo "  Summary gene count:  $N_SUMMARY"
echo ""

# Check all match
ALL_MATCH=true
for COUNT in "$N_PROT" "$N_CDS" "$N_MRNA" "$N_SUMMARY"; do
    if [[ "$COUNT" -ne "$N_GFF3_GENES" ]]; then
        ALL_MATCH=false
    fi
done

if $ALL_MATCH; then
    echo "   ALL COUNTS MATCH: $N_GFF3_GENES genes across all files"
else
    echo "   WARNING: COUNTS DO NOT MATCH — check outputs!"
    echo "    Expected $N_GFF3_GENES (from GFF3)"
fi
echo ""

# =============================================
# Summary
# =============================================
echo "============================================================"
echo "Pipeline complete. Output files:"
echo "============================================================"
echo ""
echo "RefSeq-style submission files:"
echo "  $GFF3_OUT"
echo "  $OUTPUT_DIR/${PREFIX}_protein.faa"
echo "  $OUTPUT_DIR/${PREFIX}_cds_from_genomic.fna"
echo "  $OUTPUT_DIR/${PREFIX}_rna_from_genomic.fna"
echo ""
echo "Statistics:"
echo "  $SUMMARY_OUT"
echo "  $OUTPUT_DIR/annotation_summary_duplicates.txt"
echo ""
echo "Traceability:"
echo "  $RETAINED_IDS"
echo "  $DUPLICATE_IDS"
echo ""
if [[ -f "$OUTPUT_DIR/${PREFIX}_HANNO_original.gtf" ]]; then
    echo "Reference:"
    echo "  $OUTPUT_DIR/${PREFIX}_HANNO_original.gtf"
    echo ""
fi

echo "============================================================"
