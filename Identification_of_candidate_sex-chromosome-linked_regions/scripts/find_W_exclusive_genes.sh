#!/bin/bash
# find_W_exclusive_genes.sh
# Author: Kopp K, Pristimantis euphronides genome project
#
# Pipeline to identify candidate W-exclusive genes:
#   Step 1: Cross-reference W-genic gene names against placed scaffold genes
#   Step 2: Extract sequences of unmatched + unnamed genes from the assembly
#   Step 3: Build BLAST database from placed scaffolds only (shared)
#   Step 4: BLAST unmatched W-genic sequences against placed scaffolds
#   Step 5: Report genes with no BLAST hit as W-exclusive candidates
#
# Processes both W-genic and W-genic-weak sets separately.
# Shared BLAST database and placed gene name list are built once.
#
# Usage: bash /data/GrenadaFrog144/SexChromosomes/W-chr-workflow/W_exclusive_search/scripts/find_W_exclusive_genes.sh

set -euo pipefail

# ---------------------------------------------------------------
# Input files
# ---------------------------------------------------------------
W_GENIC="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/master_table/w_genic_genes/W-genic_genes.tsv"
W_GENIC_WEAK="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/master_table/w_genic_genes/W-genic-weak_genes.tsv"
PLACED="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/gametolog_discovery_hanno7/all_placed_named.tsv"
ASSEMBLY="/data/GrenadaFrog144/Pristimantis_euphronides.genome.fasta"
BASEDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/W_exclusive_search/output_W_exclusive_search"

mkdir -p "${BASEDIR}"

# ---------------------------------------------------------------
# Verify inputs
# ---------------------------------------------------------------
for f in "${W_GENIC}" "${W_GENIC_WEAK}" "${PLACED}" "${ASSEMBLY}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: File not found: ${f}"
        exit 1
    fi
done

if [ ! -f "${ASSEMBLY}.fai" ]; then
    echo "=== Indexing assembly ==="
    samtools faidx "${ASSEMBLY}"
fi

# ---------------------------------------------------------------
# Shared: placed gene names (used by both runs)
# ---------------------------------------------------------------
awk -F'\t' '{print $1}' "${PLACED}" | sort -u > "${BASEDIR}/placed_names.txt"
P_TOTAL=$(wc -l < "${BASEDIR}/placed_names.txt")
echo "=== Placed gene names: ${P_TOTAL} ==="

# ---------------------------------------------------------------
# Shared: Build BLAST database from placed scaffolds only (once)
# ---------------------------------------------------------------
PLACED_DB="${BASEDIR}/placed_scaffolds_db"

if [ ! -f "${PLACED_DB}.ndb" ]; then
    echo "=== Building placed-scaffolds-only BLAST database ==="
    samtools faidx "${ASSEMBLY}" $(grep "^scaffold_" "${ASSEMBLY}.fai" | cut -f1) \
        > "${BASEDIR}/placed_scaffolds.fasta"

    makeblastdb \
        -in "${BASEDIR}/placed_scaffolds.fasta" \
        -dbtype nucl \
        -parse_seqids \
        -out "${PLACED_DB}"
else
    echo "=== Placed scaffolds BLAST database already exists ==="
fi

# ===============================================================
# Function: process one gene set
# ===============================================================
process_gene_set() {
    local INPUT_FILE="$1"
    local LABEL="$2"
    local OUTDIR="${BASEDIR}/${LABEL}"

    mkdir -p "${OUTDIR}"

    echo ""
    echo "################################################################"
    echo "  Processing: ${LABEL}"
    echo "  Input: ${INPUT_FILE}"
    echo "################################################################"

    # -----------------------------------------------------------
    # Step 1: Name-based cross-reference
    # -----------------------------------------------------------
    echo ""
    echo "=== ${LABEL} Step 1: Name-based cross-reference ==="

    awk -F'\t' 'NR>1 && $6!="-" {print $6}' "${INPUT_FILE}" | sort -u > "${OUTDIR}/w_genic_names.txt"
    local W_TOTAL=$(wc -l < "${OUTDIR}/w_genic_names.txt")

    comm -23 "${OUTDIR}/w_genic_names.txt" "${BASEDIR}/placed_names.txt" > "${OUTDIR}/w_no_name_match.txt"
    local W_NO_MATCH=$(wc -l < "${OUTDIR}/w_no_name_match.txt")

    comm -12 "${OUTDIR}/w_genic_names.txt" "${BASEDIR}/placed_names.txt" > "${OUTDIR}/w_genic_with_placed_copy.txt"
    local W_GAMET=$(wc -l < "${OUTDIR}/w_genic_with_placed_copy.txt")

    local W_UNNAMED=$(awk -F'\t' 'NR>1 && $6=="-"' "${INPUT_FILE}" | wc -l)

    echo "  Named genes:                       ${W_TOTAL}"
    echo "  Name match (gametologs):            ${W_GAMET}"
    echo "  No name match (named):              ${W_NO_MATCH}"
    echo "  Unnamed genes:                      ${W_UNNAMED}"
    echo "  Total for BLAST (no match+unnamed): $((W_NO_MATCH + W_UNNAMED))"

    # -----------------------------------------------------------
    # Step 2: Extract sequences
    # -----------------------------------------------------------
    echo ""
    echo "=== ${LABEL} Step 2: Extracting sequences ==="

    > "${OUTDIR}/w_genes_to_blast.bed"
    > "${OUTDIR}/w_genes_lookup.tsv"

    while read gene; do
        awk -F'\t' -v g="${gene}" 'NR>1 && $6==g {
            print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5;
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $4, $6, $1, $2, $3, $5, $8, $9 >> "/dev/stderr"
        }' "${INPUT_FILE}" >> "${OUTDIR}/w_genes_to_blast.bed" \
           2>> "${OUTDIR}/w_genes_lookup.tsv"
    done < "${OUTDIR}/w_no_name_match.txt"

    awk -F'\t' 'NR>1 && $6=="-" {
        print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5
    }' "${INPUT_FILE}" >> "${OUTDIR}/w_genes_to_blast.bed"

    awk -F'\t' 'NR>1 && $6=="-" {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $4, "-", $1, $2, $3, $5, $8, $9
    }' "${INPUT_FILE}" >> "${OUTDIR}/w_genes_lookup.tsv"

    local BLAST_COUNT=$(wc -l < "${OUTDIR}/w_genes_to_blast.bed")
    echo "  BED entries to extract: ${BLAST_COUNT}"

    bedtools getfasta \
        -fi "${ASSEMBLY}" \
        -bed "${OUTDIR}/w_genes_to_blast.bed" \
        -name \
        -s \
        -fo "${OUTDIR}/w_genes_to_blast.fasta"

    local FASTA_COUNT=$(grep -c "^>" "${OUTDIR}/w_genes_to_blast.fasta")
    echo "  Sequences extracted: ${FASTA_COUNT}"

    # -----------------------------------------------------------
    # Step 3: BLAST against placed scaffolds
    # -----------------------------------------------------------
    echo ""
    echo "=== ${LABEL} Step 3: BLASTn against placed scaffolds ==="

    blastn \
        -query "${OUTDIR}/w_genes_to_blast.fasta" \
        -db "${PLACED_DB}" \
        -evalue 1e-5 \
        -max_target_seqs 5 \
        -num_threads 12 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out "${OUTDIR}/blast_results.txt"

    cut -f1 "${OUTDIR}/blast_results.txt" | sed 's/::.*//' | sort -u \
        > "${OUTDIR}/gene_ids_with_hits.txt"
    local HITS=$(wc -l < "${OUTDIR}/gene_ids_with_hits.txt")
    echo "  Genes with at least one placed hit: ${HITS}"

    # -----------------------------------------------------------
    # Step 4: Identify W-exclusive candidates
    # -----------------------------------------------------------
    echo ""
    echo "=== ${LABEL} Step 4: Identifying W-exclusive candidates ==="

    grep "^>" "${OUTDIR}/w_genes_to_blast.fasta" | sed 's/^>//; s/::.*//' | sort -u \
        > "${OUTDIR}/all_blasted_gene_ids.txt"

    comm -23 "${OUTDIR}/all_blasted_gene_ids.txt" "${OUTDIR}/gene_ids_with_hits.txt" \
        > "${OUTDIR}/w_exclusive_gene_ids.txt"

    local W_EXCLUSIVE=$(wc -l < "${OUTDIR}/w_exclusive_gene_ids.txt")

    echo -e "gene_id\tpreferred_name\tcontig\tstart\tend\tstrand\tPFAMs\torflen" \
        > "${OUTDIR}/w_exclusive_final.tsv"

    while read gid; do
        awk -F'\t' -v id="${gid}" 'NR>1 && $4==id {
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $4, $6, $1, $2, $3, $5, $8, $9
        }' "${INPUT_FILE}"
    done < "${OUTDIR}/w_exclusive_gene_ids.txt" >> "${OUTDIR}/w_exclusive_final.tsv"

    # Diverged gametologs table
    echo -e "gene_id\tpreferred_name\tbest_hit_scaffold\tpident\talign_len\tevalue" \
        > "${OUTDIR}/w_diverged_gametologs.tsv"

    while read gid; do
        local best=$(grep "^${gid}::" "${OUTDIR}/blast_results.txt" | sort -k12,12g | head -1)
        if [ -n "${best}" ]; then
            local scaffold=$(echo "${best}" | cut -f2)
            local pident=$(echo "${best}" | cut -f3)
            local alen=$(echo "${best}" | cut -f4)
            local eval=$(echo "${best}" | cut -f11)
            local pname=$(awk -F'\t' -v id="${gid}" '$1==id {print $2; exit}' "${OUTDIR}/w_genes_lookup.tsv")
            echo -e "${gid}\t${pname}\t${scaffold}\t${pident}\t${alen}\t${eval}"
        fi
    done < <(comm -12 "${OUTDIR}/gene_ids_with_hits.txt" "${OUTDIR}/all_blasted_gene_ids.txt") \
        >> "${OUTDIR}/w_diverged_gametologs.tsv"

    # -----------------------------------------------------------
    # Report for this set
    # -----------------------------------------------------------
    echo ""
    echo "================================================================"
    echo "  ${LABEL} SUMMARY"
    echo "================================================================"
    echo ""
    echo "  Named genes:                         ${W_TOTAL}"
    echo "    Name match on placed (gametologs): ${W_GAMET}"
    echo "    No name match (named):             ${W_NO_MATCH}"
    echo "  Unnamed genes:                       ${W_UNNAMED}"
    echo ""
    echo "  Sent to BLAST:                       ${BLAST_COUNT}"
    echo "  With BLAST hit on placed:            ${HITS} (diverged gametologs)"
    echo "  NO BLAST hit (W-exclusive cand.):    ${W_EXCLUSIVE}"
    echo ""
    echo "  W-EXCLUSIVE CANDIDATES (first 50 rows):"
    # Write to temp file to avoid SIGPIPE from head under pipefail
    column -t -s$'\t' "${OUTDIR}/w_exclusive_final.tsv" > "${OUTDIR}/.tmp_display.txt"
    head -51 "${OUTDIR}/.tmp_display.txt"
    rm -f "${OUTDIR}/.tmp_display.txt"
    if [ "${W_EXCLUSIVE}" -gt 50 ]; then
        echo "  ... (${W_EXCLUSIVE} total, see w_exclusive_final.tsv)"
    fi
    echo ""
    echo "  Output in: ${OUTDIR}/"
    echo ""
}

# ===============================================================
# Run both sets
# ===============================================================
process_gene_set "${W_GENIC}" "w_genic_strong"
process_gene_set "${W_GENIC_WEAK}" "w_genic_weak"

echo ""
echo "================================================================"
echo "  ALL DONE"
echo "================================================================"
echo ""
echo "  Results:"
echo "    ${BASEDIR}/w_genic_strong/"
echo "    ${BASEDIR}/w_genic_weak/"
echo ""
