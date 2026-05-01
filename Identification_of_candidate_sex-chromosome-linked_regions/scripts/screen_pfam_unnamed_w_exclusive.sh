#!/bin/bash
# screen_pfam_unnamed_w_exclusive.sh
# Author: Kopp K, Pristimantis euphronides genome project
#
# Screens unnamed W-exclusive genes for Pfam protein domains,
# then flags any domains associated with sex determination or
# sex differentiation pathways.
#
# Input:  w_exclusive_final.tsv from both W-genic and W-genic-weak sets
# Output: List of unnamed W-exclusive genes with Pfam domains,
#         flagged for sex-related domains
#
# Usage: bash /data/GrenadaFrog144/SexChromosomes/W-chr-workflow/W_exclusive_search/scripts/screen_pfam_unnamed_w_exclusive.sh

set -euo pipefail

# ---------------------------------------------------------------
# Paths
# ---------------------------------------------------------------
BASEDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/W_exclusive_search"
STRONG_DIR="${BASEDIR}/output_W_exclusive_search/w_genic_strong"
WEAK_DIR="${BASEDIR}/output_W_exclusive_search/w_genic_weak"
OUTDIR="${BASEDIR}/output_pfam_screen"

mkdir -p "${OUTDIR}"

# ---------------------------------------------------------------
# Verify inputs
# ---------------------------------------------------------------
for f in "${STRONG_DIR}/w_exclusive_final.tsv" \
         "${WEAK_DIR}/w_exclusive_final.tsv"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: File not found: ${f}"
        exit 1
    fi
done

# ---------------------------------------------------------------
# Step 1: Extract all unnamed W-exclusive genes from both sets
# ---------------------------------------------------------------
echo "=== Step 1: Collecting unnamed W-exclusive genes ==="

echo -e "gene_id\tpreferred_name\tcontig\tstart\tend\tstrand\tPFAMs\torflen\tsource" \
    > "${OUTDIR}/unnamed_w_exclusive_all.tsv"

awk -F'\t' 'NR>1 && $2=="-" {print $0"\tw_genic"}' \
    "${STRONG_DIR}/w_exclusive_final.tsv" \
    >> "${OUTDIR}/unnamed_w_exclusive_all.tsv"

awk -F'\t' 'NR>1 && $2=="-" {print $0"\tw_genic_weak"}' \
    "${WEAK_DIR}/w_exclusive_final.tsv" \
    >> "${OUTDIR}/unnamed_w_exclusive_all.tsv"

TOTAL=$(awk 'NR>1' "${OUTDIR}/unnamed_w_exclusive_all.tsv" | wc -l)
echo "  Total unnamed W-exclusive genes: ${TOTAL}"

# ---------------------------------------------------------------
# Step 2: Separate genes with and without Pfam domains
# ---------------------------------------------------------------
echo ""
echo "=== Step 2: Pfam domain classification ==="

awk -F'\t' 'NR>1 && $7!="-"' "${OUTDIR}/unnamed_w_exclusive_all.tsv" \
    > "${OUTDIR}/unnamed_with_pfam.tsv"

awk -F'\t' 'NR>1 && $7=="-"' "${OUTDIR}/unnamed_w_exclusive_all.tsv" \
    > "${OUTDIR}/unnamed_no_pfam.tsv"

WITH_PFAM=$(wc -l < "${OUTDIR}/unnamed_with_pfam.tsv")
NO_PFAM=$(wc -l < "${OUTDIR}/unnamed_no_pfam.tsv")

echo "  With Pfam domain(s): ${WITH_PFAM}"
echo "  Without Pfam domain: ${NO_PFAM}"

# ---------------------------------------------------------------
# Step 3: Screen Pfam domains for sex-related protein families
# ---------------------------------------------------------------
echo ""
echo "=== Step 3: Screening for sex-related Pfam domains ==="

# Sex-related Pfam domain keywords:
#   DM        - DM DNA-binding domain (DMRT family, Dm-W)
#   HMG       - High mobility group box (SOX family)
#   Forkhead  - Forkhead domain (FOXL2, FOXC1)
#   zf-C2H2   - Zinc finger (many TFs including sex regulators)
#   Hormone_recep / Androgen_recep - steroid hormone receptors
#   p450      - Cytochrome P450 (steroidogenic enzymes)
#   SDR       - Short-chain dehydrogenase/reductase (HSD enzymes)
#   Aldo_ket_red - Aldo-keto reductase (steroid metabolism)
#   Wnt       - Wnt signaling
#   TGF_beta  - TGF-beta signaling
#   fn3       - Fibronectin type III (cytokine receptors incl. AMHR2)
#   DEAD      - DEAD-box helicase (RNA helicases in germ cells)
#   KH        - KH RNA-binding domain (DAZL family)
#   RRM       - RNA recognition motif (germ cell RNA regulators)
#   Tudor     - Tudor domain (piRNA pathway, TDRD family)
#   Piwi      - Piwi domain (piRNA pathway)
#   Homeobox  - Homeobox (NOBOX and other developmental TFs)

SEX_PFAM_PATTERN="DM\b|HMG|Forkhead|Hormone_recep|Androgen_recep|[Pp]450|SDR_c|Aldo_ket_red|Wnt|TGF_beta|fn3|DEAD|KH_1|RRM|Tudor|Piwi|Homeobox|zf-C4|Ank.*DM|Steroid"

if [ "${WITH_PFAM}" -gt 0 ]; then
    grep -E "${SEX_PFAM_PATTERN}" "${OUTDIR}/unnamed_with_pfam.tsv" \
        > "${OUTDIR}/sex_related_pfam_hits.tsv" || true
else
    > "${OUTDIR}/sex_related_pfam_hits.tsv"
fi

SEX_HITS=$(wc -l < "${OUTDIR}/sex_related_pfam_hits.tsv")

if [ "${SEX_HITS}" -gt 0 ]; then
    echo "  Sex-related Pfam domains found: ${SEX_HITS}"
    echo ""
    echo "  --- Sex-related Pfam hits ---"
    column -t -s$'\t' "${OUTDIR}/sex_related_pfam_hits.tsv" > "${OUTDIR}/.tmp_pfam.txt"
    cat "${OUTDIR}/.tmp_pfam.txt"
    rm -f "${OUTDIR}/.tmp_pfam.txt"
else
    echo "  No sex-related Pfam domains found."
fi

# ---------------------------------------------------------------
# Step 4: List all Pfam domains found (for completeness)
# ---------------------------------------------------------------
echo ""
echo "=== Step 4: All Pfam domains in unnamed W-exclusive genes ==="

if [ "${WITH_PFAM}" -gt 0 ]; then
    awk -F'\t' '{print $7}' "${OUTDIR}/unnamed_with_pfam.tsv" | \
        tr ',' '\n' | sort | uniq -c | sort -rn \
        > "${OUTDIR}/pfam_domain_counts.txt"

    echo "  Unique Pfam domains:"
    cat "${OUTDIR}/pfam_domain_counts.txt"
else
    echo "  (none)"
    > "${OUTDIR}/pfam_domain_counts.txt"
fi

# ---------------------------------------------------------------
# Report
# ---------------------------------------------------------------
echo ""
echo "================================================================"
echo "  SUMMARY"
echo "================================================================"
echo ""
echo "  Unnamed W-exclusive genes:         ${TOTAL}"
echo "  With Pfam domain(s):               ${WITH_PFAM}"
echo "  Without Pfam domain:               ${NO_PFAM}"
echo "  Sex-related Pfam domains:          ${SEX_HITS}"
echo ""
echo "  Output in: ${OUTDIR}/"
echo "    unnamed_w_exclusive_all.tsv       - all unnamed W-exclusive genes"
echo "    unnamed_with_pfam.tsv             - subset with Pfam domains"
echo "    unnamed_no_pfam.tsv               - subset without Pfam domains"
echo "    sex_related_pfam_hits.tsv         - sex-domain flagged genes"
echo "    pfam_domain_counts.txt            - domain frequency table"
echo ""
echo "=== DONE ==="
