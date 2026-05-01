#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# gametolog_blastn.sh
# ====================
# Pairwise blastn alignment of Z- vs W-gametolog copies for all tiers.
#
# For each gametolog gene (tiers 1a–4, which all have exactly 1 Hemi placed
# + 1 Hemi unplaced copy), extract the genomic region from the unscaffolded
# assembly and align Z (query) against W (subject).
#
# Z as query because: Z-copies are expected to be more complete (better
# assembled placed contigs). Query coverage then reports what fraction of
# the Z-gene is retained on the W-copy. Low qcov = W truncation/degeneration.
#
# Tier 5 (multiple Hemi on one/both sides) is handled separately: all
# pairwise combinations of placed-Hemi × unplaced-Hemi are aligned.
#
# Input:
#   - gene_summary.tsv (from gametolog_discovery_hanno7)
#   - placed_genes.tsv / unplaced_genes.tsv (gene coordinates)
#   - Unscaffolded assembly FASTA (indexed with samtools faidx)
#
# Output:
#   - gametolog_blastn_results.tsv (all pairs, all tiers)
#   - per-tier summary statistics
#   - FASTA files for Z and W regions (in tmp dir, cleaned up unless --keep)
#
# Date: 2026-03-01
# Usage: bash gametolog_blastn.sh

set -euo pipefail

# ============================================================
# PATHS
# ============================================================
GAMETOLOG_DIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/gametolog_discovery_hanno7"
ASSEMBLY="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"

GENE_SUMMARY="${GAMETOLOG_DIR}/gene_summary.tsv"
PLACED_GENES="${GAMETOLOG_DIR}/placed_genes.tsv"
UNPLACED_GENES="${GAMETOLOG_DIR}/unplaced_genes.tsv"

OUTDIR="${GAMETOLOG_DIR}/blastn"
TMPDIR="${OUTDIR}/tmp"

# ============================================================
# SANITY CHECKS
# ============================================================
echo "=== Gametolog Pairwise Blastn ==="
echo "Start: $(date)"
echo ""

for f in "$GENE_SUMMARY" "$PLACED_GENES" "$UNPLACED_GENES" "$ASSEMBLY" "${ASSEMBLY}.fai"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

command -v blastn >/dev/null 2>&1 || { echo "ERROR: blastn not found in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH"; exit 1; }

mkdir -p "$OUTDIR" "$TMPDIR"

# ============================================================
# STEP 1: Build gene-to-coordinates lookup for Hemi copies
# ============================================================
echo "[Step 1] Extracting Hemi gene coordinates ..."

# From placed_genes.tsv (15 cols):
#   col1=contig, col2=gene_start(0-based), col3=gene_end, col4=gene_id,
#   col5=strand, col6=gene_name, col14=cov_class, col15=z_flag
# Filter: Hemi_0.5x only, named genes only

awk -F'\t' '
    NR > 1 && $6 != "UNNAMED" && $14 == "Hemi_0.5x" {
        # Convert 0-based start to 1-based for samtools faidx
        print $6 "\t" $1 "\t" $2+1 "\t" $3 "\t" $5 "\t" $14 "\t" $15
    }
' "$PLACED_GENES" > "${TMPDIR}/placed_hemi_coords.tsv"

# From unplaced_genes.tsv (11 cols):
#   col1=contig, col2=gene_start(0-based), col3=gene_end, col4=gene_id,
#   col5=strand, col6=gene_name, col11=cov_class

awk -F'\t' '
    NR > 1 && $6 != "UNNAMED" && $11 == "Hemi_0.5x" {
        print $6 "\t" $1 "\t" $2+1 "\t" $3 "\t" $5 "\t" $11
    }
' "$UNPLACED_GENES" > "${TMPDIR}/unplaced_hemi_coords.tsv"

N_PLACED_HEMI=$(wc -l < "${TMPDIR}/placed_hemi_coords.tsv")
N_UNPLACED_HEMI=$(wc -l < "${TMPDIR}/unplaced_hemi_coords.tsv")
echo "  Placed Hemi gene entries: ${N_PLACED_HEMI}"
echo "  Unplaced Hemi gene entries: ${N_UNPLACED_HEMI}"

# ============================================================
# STEP 2: Build pairs list from gene_summary.tsv
# ============================================================
echo "[Step 2] Building gametolog pair list ..."

# For tiers 1a–4: exactly 1 placed Hemi + 1 unplaced Hemi per gene.
# For tier 5: multiple Hemi on one/both sides — handled below.
# gene_summary.tsv col12 = tier

# Extract tier 1a–4 gene names and their tiers
awk -F'\t' '
    NR > 1 && $12 ~ /^tier[1-4]/ {
        print $1 "\t" $12
    }
' "$GENE_SUMMARY" > "${TMPDIR}/genes_to_align.tsv"

# Extract tier 5 gene names
awk -F'\t' '
    NR > 1 && $12 == "tier5_multiHemi" {
        print $1 "\t" $12
    }
' "$GENE_SUMMARY" > "${TMPDIR}/genes_tier5.tsv"

N_T1_4=$(wc -l < "${TMPDIR}/genes_to_align.tsv")
N_T5=$(wc -l < "${TMPDIR}/genes_tier5.tsv")
echo "  Genes in tiers 1–4: ${N_T1_4}"
echo "  Genes in tier 5: ${N_T5}"

# ============================================================
# STEP 3: Extract sequences and run blastn for tiers 1a–4
# ============================================================
echo "[Step 3] Running pairwise blastn (tiers 1a–4) ..."

# Output header
echo -e "gene_name\ttier\tZ_contig\tZ_start\tZ_end\tZ_strand\tZ_length\tW_contig\tW_start\tW_end\tW_strand\tW_length\tpct_identity\talignment_length\tmismatches\tgap_opens\tq_start\tq_end\ts_start\ts_end\tevalue\tbitscore\tquery_coverage" \
    > "${OUTDIR}/gametolog_blastn_results.tsv"

PAIR_COUNT=0
NO_HIT_COUNT=0
NO_COORDS_COUNT=0

while IFS=$'\t' read -r GENE TIER; do
    # Look up placed Hemi coordinates for this gene
    Z_LINE=$(awk -F'\t' -v g="$GENE" '$1 == g' "${TMPDIR}/placed_hemi_coords.tsv" | head -1)
    W_LINE=$(awk -F'\t' -v g="$GENE" '$1 == g' "${TMPDIR}/unplaced_hemi_coords.tsv" | head -1)

    if [ -z "$Z_LINE" ] || [ -z "$W_LINE" ]; then
        NO_COORDS_COUNT=$((NO_COORDS_COUNT + 1))
        continue
    fi

    Z_CONTIG=$(echo "$Z_LINE" | cut -f2)
    Z_START=$(echo "$Z_LINE" | cut -f3)
    Z_END=$(echo "$Z_LINE" | cut -f4)
    Z_STRAND=$(echo "$Z_LINE" | cut -f5)
    Z_LEN=$((Z_END - Z_START + 1))

    W_CONTIG=$(echo "$W_LINE" | cut -f2)
    W_START=$(echo "$W_LINE" | cut -f3)
    W_END=$(echo "$W_LINE" | cut -f4)
    W_STRAND=$(echo "$W_LINE" | cut -f5)
    W_LEN=$((W_END - W_START + 1))

    # Extract sequences
    Z_FA="${TMPDIR}/${GENE}_Z.fa"
    W_FA="${TMPDIR}/${GENE}_W.fa"

    samtools faidx "$ASSEMBLY" "${Z_CONTIG}:${Z_START}-${Z_END}" > "$Z_FA" 2>/dev/null
    samtools faidx "$ASSEMBLY" "${W_CONTIG}:${W_START}-${W_END}" > "$W_FA" 2>/dev/null

    # Check extraction worked
    Z_SEQ_LEN=$(awk '/^[^>]/ {s += length($0)} END {print s+0}' "$Z_FA")
    W_SEQ_LEN=$(awk '/^[^>]/ {s += length($0)} END {print s+0}' "$W_FA")

    if [ "$Z_SEQ_LEN" -eq 0 ] || [ "$W_SEQ_LEN" -eq 0 ]; then
        NO_COORDS_COUNT=$((NO_COORDS_COUNT + 1))
        rm -f "$Z_FA" "$W_FA"
        continue
    fi

    # Run blastn: Z as query, W as subject
    BLAST_RESULT=$(blastn -query "$Z_FA" -subject "$W_FA" \
        -outfmt "6 pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 \
        2>/dev/null | head -1)

    if [ -n "$BLAST_RESULT" ]; then
        # Parse blast fields
        PIDENT=$(echo "$BLAST_RESULT" | cut -f1)
        ALEN=$(echo "$BLAST_RESULT" | cut -f2)
        MISMATCH=$(echo "$BLAST_RESULT" | cut -f3)
        GAPOPEN=$(echo "$BLAST_RESULT" | cut -f4)
        QSTART=$(echo "$BLAST_RESULT" | cut -f5)
        QEND=$(echo "$BLAST_RESULT" | cut -f6)
        SSTART=$(echo "$BLAST_RESULT" | cut -f7)
        SEND=$(echo "$BLAST_RESULT" | cut -f8)
        EVALUE=$(echo "$BLAST_RESULT" | cut -f9)
        BITSCORE=$(echo "$BLAST_RESULT" | cut -f10)

        # Query coverage = (qend - qstart + 1) / query length
        QCOV=$(awk -v qs="$QSTART" -v qe="$QEND" -v ql="$Z_LEN" 'BEGIN { printf "%.2f", ((qe-qs+1)/ql)*100 }')

        echo -e "${GENE}\t${TIER}\t${Z_CONTIG}\t${Z_START}\t${Z_END}\t${Z_STRAND}\t${Z_LEN}\t${W_CONTIG}\t${W_START}\t${W_END}\t${W_STRAND}\t${W_LEN}\t${PIDENT}\t${ALEN}\t${MISMATCH}\t${GAPOPEN}\t${QSTART}\t${QEND}\t${SSTART}\t${SEND}\t${EVALUE}\t${BITSCORE}\t${QCOV}" \
            >> "${OUTDIR}/gametolog_blastn_results.tsv"

        PAIR_COUNT=$((PAIR_COUNT + 1))
    else
        # No significant hit
        echo -e "${GENE}\t${TIER}\t${Z_CONTIG}\t${Z_START}\t${Z_END}\t${Z_STRAND}\t${Z_LEN}\t${W_CONTIG}\t${W_START}\t${W_END}\t${W_STRAND}\t${W_LEN}\tno_hit\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0.00" \
            >> "${OUTDIR}/gametolog_blastn_results.tsv"

        NO_HIT_COUNT=$((NO_HIT_COUNT + 1))
        PAIR_COUNT=$((PAIR_COUNT + 1))
    fi

    rm -f "$Z_FA" "$W_FA"

done < "${TMPDIR}/genes_to_align.tsv"

echo "  Tiers 1–4 aligned: ${PAIR_COUNT}"
echo "  No coordinates found: ${NO_COORDS_COUNT}"
echo "  No significant hit: ${NO_HIT_COUNT}"

# ============================================================
# STEP 4: Handle tier 5 (multiple Hemi — all pairwise)
# ============================================================
echo "[Step 4] Running pairwise blastn (tier 5 — all Hemi combinations) ..."

T5_PAIR_COUNT=0
T5_NO_HIT=0

while IFS=$'\t' read -r GENE TIER; do
    # Get ALL placed Hemi entries for this gene
    Z_LINES=$(awk -F'\t' -v g="$GENE" '$1 == g' "${TMPDIR}/placed_hemi_coords.tsv")
    W_LINES=$(awk -F'\t' -v g="$GENE" '$1 == g' "${TMPDIR}/unplaced_hemi_coords.tsv")

    if [ -z "$Z_LINES" ] || [ -z "$W_LINES" ]; then
        continue
    fi

    # Iterate all Z × W combinations
    while IFS=$'\t' read -r _GN Z_CONTIG Z_START Z_END Z_STRAND _COV _ZFLAG; do
        Z_LEN=$((Z_END - Z_START + 1))
        Z_FA="${TMPDIR}/${GENE}_${Z_CONTIG}_Z.fa"
        samtools faidx "$ASSEMBLY" "${Z_CONTIG}:${Z_START}-${Z_END}" > "$Z_FA" 2>/dev/null

        while IFS=$'\t' read -r _GN2 W_CONTIG W_START W_END W_STRAND _COV2; do
            W_LEN=$((W_END - W_START + 1))
            W_FA="${TMPDIR}/${GENE}_${W_CONTIG}_W.fa"
            samtools faidx "$ASSEMBLY" "${W_CONTIG}:${W_START}-${W_END}" > "$W_FA" 2>/dev/null

            BLAST_RESULT=$(blastn -query "$Z_FA" -subject "$W_FA" \
                -outfmt "6 pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
                -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 \
                2>/dev/null | head -1)

            if [ -n "$BLAST_RESULT" ]; then
                PIDENT=$(echo "$BLAST_RESULT" | cut -f1)
                ALEN=$(echo "$BLAST_RESULT" | cut -f2)
                MISMATCH=$(echo "$BLAST_RESULT" | cut -f3)
                GAPOPEN=$(echo "$BLAST_RESULT" | cut -f4)
                QSTART=$(echo "$BLAST_RESULT" | cut -f5)
                QEND=$(echo "$BLAST_RESULT" | cut -f6)
                SSTART=$(echo "$BLAST_RESULT" | cut -f7)
                SEND=$(echo "$BLAST_RESULT" | cut -f8)
                EVALUE=$(echo "$BLAST_RESULT" | cut -f9)
                BITSCORE=$(echo "$BLAST_RESULT" | cut -f10)
                QCOV=$(awk -v qs="$QSTART" -v qe="$QEND" -v ql="$Z_LEN" 'BEGIN { printf "%.2f", ((qe-qs+1)/ql)*100 }')

                echo -e "${GENE}\t${TIER}\t${Z_CONTIG}\t${Z_START}\t${Z_END}\t${Z_STRAND}\t${Z_LEN}\t${W_CONTIG}\t${W_START}\t${W_END}\t${W_STRAND}\t${W_LEN}\t${PIDENT}\t${ALEN}\t${MISMATCH}\t${GAPOPEN}\t${QSTART}\t${QEND}\t${SSTART}\t${SEND}\t${EVALUE}\t${BITSCORE}\t${QCOV}" \
                    >> "${OUTDIR}/gametolog_blastn_results.tsv"
                T5_PAIR_COUNT=$((T5_PAIR_COUNT + 1))
            else
                echo -e "${GENE}\t${TIER}\t${Z_CONTIG}\t${Z_START}\t${Z_END}\t${Z_STRAND}\t${Z_LEN}\t${W_CONTIG}\t${W_START}\t${W_END}\t${W_STRAND}\t${W_LEN}\tno_hit\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0.00" \
                    >> "${OUTDIR}/gametolog_blastn_results.tsv"
                T5_NO_HIT=$((T5_NO_HIT + 1))
                T5_PAIR_COUNT=$((T5_PAIR_COUNT + 1))
            fi

            rm -f "$W_FA"
        done <<< "$W_LINES"

        rm -f "$Z_FA"
    done <<< "$Z_LINES"

done < "${TMPDIR}/genes_tier5.tsv"

echo "  Tier 5 pairs aligned: ${T5_PAIR_COUNT}"
echo "  Tier 5 no hit: ${T5_NO_HIT}"

# ============================================================
# STEP 5: Summary report
# ============================================================
echo "[Step 5] Generating summary report ..."

TOTAL_PAIRS=$((PAIR_COUNT + T5_PAIR_COUNT))

{
    echo "======================================================================"
    echo "GAMETOLOG PAIRWISE BLASTN REPORT"
    echo "======================================================================"
    echo "Date: $(date)"
    echo ""
    echo "Total pairs aligned: ${TOTAL_PAIRS}"
    echo "  Tiers 1–4: ${PAIR_COUNT}"
    echo "  Tier 5:    ${T5_PAIR_COUNT}"
    echo "  No coords: ${NO_COORDS_COUNT}"
    echo ""
    echo "--- Hit / no-hit breakdown by tier ---"
    echo ""
    echo -e "tier\ttotal\thit\tno_hit\tpct_hit"

    for T in $(awk -F'\t' 'NR > 1 { print $2 }' "${OUTDIR}/gametolog_blastn_results.tsv" | sort -u); do
        awk -F'\t' -v tier="$T" '
            NR > 1 && $2 == tier {
                total++
                if ($13 == "no_hit") nohit++
                else hit++
            }
            END {
                if (total == 0) exit
                printf "%s\t%d\t%d\t%d\t%.1f%%\n", tier, total, hit+0, nohit+0, (hit+0)/total*100
            }
        ' "${OUTDIR}/gametolog_blastn_results.tsv"
    done

    echo ""
    echo "--- Per-tier summary (excluding no_hit) ---"
    echo ""
    echo "tier	n_pairs	mean_identity	median_identity	mean_qcov	median_qcov	min_identity	max_identity"

    # Get unique tiers with hits, sorted
    TIERS_WITH_HITS=$(awk -F'\t' 'NR > 1 && $13 != "no_hit" { print $2 }' \
        "${OUTDIR}/gametolog_blastn_results.tsv" | sort -u)

    for T in $TIERS_WITH_HITS; do
        awk -F'\t' -v tier="$T" '
            NR > 1 && $2 == tier && $13 != "no_hit" {
                id[++n] = $13+0
                qc[n] = $23+0
                sum_id += $13+0
                sum_qc += $23+0
                if (n == 1 || $13+0 < min_id) min_id = $13+0
                if (n == 1 || $13+0 > max_id) max_id = $13+0
            }
            END {
                if (n == 0) exit
                mean_id = sum_id / n
                mean_qc = sum_qc / n

                # Sort identity values for median (simple insertion sort)
                for (i = 2; i <= n; i++) {
                    v = id[i]; j = i - 1
                    while (j >= 1 && id[j] > v) { id[j+1] = id[j]; j-- }
                    id[j+1] = v
                }
                for (i = 2; i <= n; i++) {
                    v = qc[i]; j = i - 1
                    while (j >= 1 && qc[j] > v) { qc[j+1] = qc[j]; j-- }
                    qc[j+1] = v
                }

                if (n % 2 == 1) { med_id = id[int(n/2)+1]; med_qc = qc[int(n/2)+1] }
                else { med_id = (id[n/2] + id[n/2+1]) / 2; med_qc = (qc[n/2] + qc[n/2+1]) / 2 }

                printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", \
                    tier, n, mean_id, med_id, mean_qc, med_qc, min_id, max_id
            }
        ' "${OUTDIR}/gametolog_blastn_results.tsv"
    done

    echo ""
    echo "--- No-hit pairs ---"
    awk -F'\t' 'NR > 1 && $13 == "no_hit" { print $1 "\t" $2 }' \
        "${OUTDIR}/gametolog_blastn_results.tsv"
    echo ""

    echo "--- Identity distribution (all tiers with hits) ---"
    echo ""
    echo "Range           Count"
    awk -F'\t' '
        NR > 1 && $13 != "no_hit" {
            id = $13+0
            if (id >= 95)      bin["95-100"]++
            else if (id >= 90) bin["90-95"]++
            else if (id >= 85) bin["85-90"]++
            else if (id >= 80) bin["80-85"]++
            else if (id >= 70) bin["70-80"]++
            else if (id >= 60) bin["60-70"]++
            else               bin["<60"]++
            total++
        }
        END {
            split("95-100,90-95,85-90,80-85,70-80,60-70,<60", ranges, ",")
            for (i = 1; i <= 7; i++) {
                r = ranges[i]
                printf "%-15s %d\n", r, bin[r]+0
            }
        }
    ' "${OUTDIR}/gametolog_blastn_results.tsv"

    echo ""
    echo "======================================================================"
    echo "Output: ${OUTDIR}/gametolog_blastn_results.tsv"
    echo "======================================================================"

} > "${OUTDIR}/gametolog_blastn_report.txt"

cat "${OUTDIR}/gametolog_blastn_report.txt"

# ============================================================
# CLEANUP
# ============================================================
rm -rf "$TMPDIR"

echo ""
echo "=== Gametolog Blastn complete ==="
