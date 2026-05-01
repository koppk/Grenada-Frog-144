#!/bin/bash
# blast_unnamed_w_exclusive_vs_nt.sh
# Author: Kopp K, Pristimantis euphronides genome project
#
# Extracts unnamed W-exclusive gene sequences from the find_W_exclusive_genes.sh
# output (both W-genic and W-genic-weak sets) and queries them against a local
# NCBI nt database to identify homologs undetectable by the HANNO pipeline.
#
# Input:  w_exclusive_final.tsv and w_genes_to_blast.fasta from both gene sets
# Output: BLASTn tabular results with subject titles for manual inspection
#
# Usage: bash /data/GrenadaFrog144/SexChromosomes/W-chr-workflow/W_exclusive_search/scripts/blast_unnamed_w_exclusive_vs_nt.sh

set -euo pipefail

# ---------------------------------------------------------------
# Paths
# ---------------------------------------------------------------
BASEDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/W_exclusive_search"
STRONG_DIR="${BASEDIR}/output_W_exclusive_search/w_genic_strong"
WEAK_DIR="${BASEDIR}/output_W_exclusive_search/w_genic_weak"
OUTDIR="${BASEDIR}/output_unnamed_w_exclusive_vs_nt"
NT_DB="/data/NCBI_nt/nt"
THREADS=24

mkdir -p "${OUTDIR}"

# ---------------------------------------------------------------
# Verify inputs
# ---------------------------------------------------------------
for f in "${STRONG_DIR}/w_exclusive_final.tsv" \
         "${STRONG_DIR}/w_genes_to_blast.fasta" \
         "${WEAK_DIR}/w_exclusive_final.tsv" \
         "${WEAK_DIR}/w_genes_to_blast.fasta"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: File not found: ${f}"
        exit 1
    fi
done

if [ ! -f "${NT_DB}.nhd" ] && [ ! -f "${NT_DB}.000.nhd" ]; then
    echo "ERROR: NCBI nt database not found at ${NT_DB}"
    exit 1
fi

# ---------------------------------------------------------------
# Step 1: Collect unnamed W-exclusive gene IDs from both sets
# ---------------------------------------------------------------
echo "=== Step 1: Collecting unnamed W-exclusive gene IDs ==="

awk -F'\t' 'NR>1 && $2=="-" {print $1}' \
    "${STRONG_DIR}/w_exclusive_final.tsv" \
    "${WEAK_DIR}/w_exclusive_final.tsv" \
    | sort -u > "${OUTDIR}/unnamed_w_exclusive_ids.txt"

TOTAL_IDS=$(wc -l < "${OUTDIR}/unnamed_w_exclusive_ids.txt")
echo "  Unnamed W-exclusive gene IDs: ${TOTAL_IDS}"

# ---------------------------------------------------------------
# Step 2: Extract sequences from existing FASTA files
# ---------------------------------------------------------------
echo ""
echo "=== Step 2: Extracting sequences ==="

# FASTA headers from bedtools getfasta have format:
#   >hanno.g111.t1::contig_3549:87830-88097(+)
# Gene IDs in the filter list are bare: hanno.g111.t1
# Convert IDs to regex patterns (escape dots, anchor to start of header)
sed 's/\./\\./g; s/$/::/' "${OUTDIR}/unnamed_w_exclusive_ids.txt" \
    > "${OUTDIR}/unnamed_w_exclusive_patterns.txt"

cat "${STRONG_DIR}/w_genes_to_blast.fasta" \
    "${WEAK_DIR}/w_genes_to_blast.fasta" \
    > "${OUTDIR}/all_w_genes_combined.fasta"

seqkit grep -r -f "${OUTDIR}/unnamed_w_exclusive_patterns.txt" \
    "${OUTDIR}/all_w_genes_combined.fasta" \
    > "${OUTDIR}/unnamed_w_exclusive.fasta"

FASTA_COUNT=$(grep -c "^>" "${OUTDIR}/unnamed_w_exclusive.fasta")
echo "  Sequences extracted: ${FASTA_COUNT}"

rm -f "${OUTDIR}/all_w_genes_combined.fasta"

# ---------------------------------------------------------------
# Step 3: BLASTn against local NCBI nt
# ---------------------------------------------------------------
echo ""
echo "=== Step 3: BLASTn against NCBI nt (${THREADS} threads) ==="
echo "  This may take several hours."

blastn \
    -query "${OUTDIR}/unnamed_w_exclusive.fasta" \
    -db "${NT_DB}" \
    -evalue 1e-5 \
    -max_target_seqs 3 \
    -num_threads "${THREADS}" \
    -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "${OUTDIR}/unnamed_w_exclusive_vs_nt.txt"

TOTAL_HITS=$(wc -l < "${OUTDIR}/unnamed_w_exclusive_vs_nt.txt")

if [ "${TOTAL_HITS}" -eq 0 ]; then
    echo "  No hits found against NCBI nt."
    QUERIES_WITH_HITS=0
    QUERIES_NO_HIT=${TOTAL_IDS}
else
    QUERIES_WITH_HITS=$(cut -f1 "${OUTDIR}/unnamed_w_exclusive_vs_nt.txt" | sed 's/::.*//' | sort -u | wc -l)
    QUERIES_NO_HIT=$(comm -23 \
        "${OUTDIR}/unnamed_w_exclusive_ids.txt" \
        <(cut -f1 "${OUTDIR}/unnamed_w_exclusive_vs_nt.txt" | sed 's/::.*//' | sort -u) \
        | wc -l)
fi

echo "  Total alignments: ${TOTAL_HITS}"
echo "  Genes with at least one hit: ${QUERIES_WITH_HITS}"
echo "  Genes with no hit: ${QUERIES_NO_HIT}"

# ---------------------------------------------------------------
# Step 4: Flag sex-related hits
# ---------------------------------------------------------------
echo ""
echo "=== Step 4: Screening for sex-related hits ==="

grep -i "sex\|dmrt\|steroid\|oocyte\|gonad\|ovary\|testis\|aromatase\|estrogen\|androgen\|follicle\|germ.cell\|spermat\|meiosis\|fertiliz" \
    "${OUTDIR}/unnamed_w_exclusive_vs_nt.txt" \
    > "${OUTDIR}/sex_related_hits.txt" || true

SEX_HITS=$(wc -l < "${OUTDIR}/sex_related_hits.txt")

if [ "${SEX_HITS}" -gt 0 ]; then
    echo "  Sex-related hits found: ${SEX_HITS}"
    echo ""
    echo "  --- Sex-related hits ---"
    column -t -s$'\t' "${OUTDIR}/sex_related_hits.txt" > "${OUTDIR}/.tmp_sex.txt"
    cat "${OUTDIR}/.tmp_sex.txt"
    rm -f "${OUTDIR}/.tmp_sex.txt"
else
    echo "  No sex-related hits found."
fi

# ---------------------------------------------------------------
# Step 5: Summary of best hits per gene (top 1 per query)
# ---------------------------------------------------------------
echo ""
echo "=== Step 5: Best hit per gene ==="

echo -e "gene_id\tsseqid\tstitle\tpident\tlength\tevalue\tbitscore" \
    > "${OUTDIR}/best_hits_per_gene.tsv"

if [ "${TOTAL_HITS}" -gt 0 ]; then
    cut -f1 "${OUTDIR}/unnamed_w_exclusive_vs_nt.txt" | sed 's/::.*//' | sort -u | \
    while read gid; do
        grep -F "${gid}::" "${OUTDIR}/unnamed_w_exclusive_vs_nt.txt" | sort -k13,13gr | head -1 | \
            awk -F'\t' -v id="${gid}" '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", id, $2, $3, $4, $5, $12, $13}'
    done >> "${OUTDIR}/best_hits_per_gene.tsv"
fi

BEST_COUNT=$(awk 'NR>1' "${OUTDIR}/best_hits_per_gene.tsv" | wc -l)
echo "  Best-hit table: ${BEST_COUNT} genes (see best_hits_per_gene.tsv)"

# ---------------------------------------------------------------
# Report
# ---------------------------------------------------------------
echo ""
echo "================================================================"
echo "  SUMMARY"
echo "================================================================"
echo ""
echo "  Unnamed W-exclusive genes queried: ${TOTAL_IDS}"
echo "  Sequences extracted:               ${FASTA_COUNT}"
echo "  Genes with NCBI nt hit:            ${QUERIES_WITH_HITS}"
echo "  Genes with no hit anywhere:        ${QUERIES_NO_HIT}"
echo "  Sex-related hits:                  ${SEX_HITS}"
echo ""
echo "  Output in: ${OUTDIR}/"
echo "    unnamed_w_exclusive.fasta        - query sequences"
echo "    unnamed_w_exclusive_vs_nt.txt    - full BLASTn results"
echo "    best_hits_per_gene.tsv           - best hit per gene"
echo "    sex_related_hits.txt             - sex-keyword-flagged hits"
echo "    unnamed_w_exclusive_ids.txt      - gene ID list"
echo ""
echo "=== DONE ==="
