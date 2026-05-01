#!/bin/bash
set -euo pipefail

# =============================================================================
# run_ncbi_reformat.sh
#
# Wrapper: HANNO annotation → NCBI-compliant files (Zenodo + table2asn)
#
# Usage:
#   bash run_ncbi_reformat.sh \
#       -i annotation_output/ \
#       -g genome_scaffolded.fasta \
#       -o ncbi_annotation/ \
#       [--locus-prefix PRIEUP] \
#       [--dbname KoppSGU] \
#       [--assembly-report GCA_965278355.2_aPriEup1.0_assembly_report.txt]
#
# Author: K. Kopp, P. euphronides genome project
# =============================================================================

usage() {
    echo "Usage: $0 -i ANNO_DIR -g GENOME_FASTA -o OUTPUT_DIR [options]"
    echo ""
    echo "Required:"
    echo "  -i  Annotation pipeline output directory"
    echo "  -g  Genome FASTA (for scaffold name extraction)"
    echo "  -o  Output directory"
    echo ""
    echo "Optional:"
    echo "  --locus-prefix PREFIX       Registered locus_tag prefix (default: PRIEUP)"
    echo "  --dbname NAME               Lab identifier for table2asn (default: KoppSGU)"
    echo "  --assembly-report FILE      NCBI assembly_report.txt (preferred for mapping)"
    echo "  --wgs-prefix PREFIX         WGS accession prefix (default: CBDIFN02)"
    echo "  --assembly-prefix NAME      Assembly name (default: aPriEup1.0)"
    echo "  --mode MODE                 zenodo|table2asn|both (default: both)"
    echo "  -s DIR                      Scripts directory"
    exit 1
}

LOCUS_PREFIX="PRIEUP"
DBNAME="KoppSGU"
WGS_PREFIX="CBDIFN02"
ASSEMBLY_PREFIX="aPriEup1.0"
ASSEMBLY_REPORT=""
MODE="both"
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
ANNO_DIR=""
GENOME=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) ANNO_DIR="$2"; shift 2 ;;
        -g) GENOME="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        --locus-prefix) LOCUS_PREFIX="$2"; shift 2 ;;
        --dbname) DBNAME="$2"; shift 2 ;;
        --assembly-report) ASSEMBLY_REPORT="$2"; shift 2 ;;
        --wgs-prefix) WGS_PREFIX="$2"; shift 2 ;;
        --assembly-prefix) ASSEMBLY_PREFIX="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        -s) SCRIPTS_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$ANNO_DIR" || -z "$GENOME" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: -i, -g, and -o are required."
    usage
fi

mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "HANNO → NCBI Reformatter"
echo "============================================================"
echo "  Mode:              $MODE"
echo "  Annotation dir:    $ANNO_DIR"
echo "  Genome FASTA:      $GENOME"
echo "  Locus prefix:      $LOCUS_PREFIX"
if [[ "$MODE" == "table2asn" || "$MODE" == "both" ]]; then
    echo "  DB name:           $DBNAME (table2asn only)"
fi
echo "  Output dir:        $OUTPUT_DIR"
echo "============================================================"

# --- Step 1: Scaffold mapping ---
echo ""
echo ">>> Step 1: Generating scaffold → accession mapping..."
SEQID_MAP="$OUTPUT_DIR/scaffold_to_accession.tsv"

if [[ -n "$ASSEMBLY_REPORT" ]]; then
    python3 "$SCRIPTS_DIR/generate_seqid_map.py" \
        --assembly-report "$ASSEMBLY_REPORT" -o "$SEQID_MAP"
else
    python3 "$SCRIPTS_DIR/generate_seqid_map.py" \
        --genome "$GENOME" --wgs-prefix "$WGS_PREFIX" -o "$SEQID_MAP"
fi

# --- Step 2: Reformat ---
echo ""
echo ">>> Step 2: Reformatting annotation..."

GFF3_IN="$ANNO_DIR/${ASSEMBLY_PREFIX}_genomic.gff3"
CDS_IN="$ANNO_DIR/${ASSEMBLY_PREFIX}_cds_from_genomic.fna"
PROT_IN="$ANNO_DIR/${ASSEMBLY_PREFIX}_protein.faa"
MRNA_IN="$ANNO_DIR/${ASSEMBLY_PREFIX}_rna_from_genomic.fna"

CMD="python3 $SCRIPTS_DIR/reformat_hanno_to_ncbi.py"
CMD="$CMD --gff3 $GFF3_IN --seqid-map $SEQID_MAP"
CMD="$CMD --locus-prefix $LOCUS_PREFIX --dbname $DBNAME"
CMD="$CMD --assembly-prefix $ASSEMBLY_PREFIX --mode $MODE -o $OUTPUT_DIR"

[[ -f "$CDS_IN" ]] && CMD="$CMD --cds $CDS_IN"
[[ -f "$PROT_IN" ]] && CMD="$CMD --protein $PROT_IN"
[[ -f "$MRNA_IN" ]] && CMD="$CMD --mrna $MRNA_IN"

eval "$CMD"

# --- Step 3: Gzip for Zenodo ---
if [[ "$MODE" == "zenodo" || "$MODE" == "both" ]]; then
    ZENODO_DIR="$OUTPUT_DIR"
    [[ "$MODE" == "both" ]] && ZENODO_DIR="$OUTPUT_DIR/zenodo"

    echo ""
    echo ">>> Step 3: Creating gzipped copies..."
    for f in "$ZENODO_DIR"/${ASSEMBLY_PREFIX}_*.{gff,fna,faa} 2>/dev/null; do
        if [[ -f "$f" ]]; then
            gzip -kf "$f"
            echo "  $(basename "$f").gz"
        fi
    done
fi

echo ""
echo "============================================================"
echo "Done. See GUIDE_ncbi_submission.md for next steps."
echo "============================================================"
