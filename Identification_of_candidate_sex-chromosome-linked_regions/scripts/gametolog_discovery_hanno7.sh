#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# gametolog_discovery_hanno7.sh
# =============================
# Identify candidate gametolog pairs from a single HANNO annotation run
# (7 Hyloidea protein + mRNA references) on the full unscaffolded
# P. euphronides primary assembly.
#
# Logic:
#   1. Split HANNO gene models into placed / unplaced by contig name
#   2. For placed genes: map contig → scaffold + position via AGP,
#      flag Z-candidate regions (scaffold_2: 0-42 Mb, scaffold_8: 47-64 Mb)
#   3. Add coverage class (Hemi_0.5x / Auto_1.0x / etc.) to all genes
#   4. For each named gene: count placed Hemi, placed Hemi in Z-region,
#      unplaced Hemi, total placed, total unplaced
#   5. Classify into tiers:
#        Tier 1a: strictly 1 placed (Hemi, Z-region) + 1 unplaced (Hemi),
#                 no paralogs, clean 1:1 contig pair, >=2 genes on pair.
#        Tier 1b: same as 1a but only 1 gene on the contig pair.
#        Tier 1c: gene-level strict 1:1 Hemi Z + Hemi unplaced, no paralogs,
#                 BUT Z-contig is associated with multiple W-contigs (or vice
#                 versa) via different gene names — no 1:1 contig pair.
#        Tier 2:  strictly 1 placed (Hemi, any scaffold) + 1 unplaced (Hemi),
#                 no paralogs.
#        Tier 3:  1 Hemi placed in Z-region + 1 Hemi unplaced, BUT additional
#                 paralog(s) on placed and/or unplaced contigs.
#        Tier 4:  1 Hemi placed (not Z-region) + 1 Hemi unplaced, BUT additional
#                 paralog(s) on placed and/or unplaced contigs.
#        Tier 5:  multiple Hemi copies (>1) on either/both sides (gene family)
#        Tier 6:  Hemi on one side only (speculative)
#        Tier 7:  no Hemi involvement
#   6. Output detailed tables and summary report
#
# A gene qualifies as a gametolog candidate ONLY if found on BOTH a placed
# Hemi contig AND an unplaced Hemi contig. The Z-region criterion strengthens
# candidates but is not required for tier 2.
#
# Input:
#   - HANNO BESTMODELS-FINAL.bedDB (from full unscaffolded assembly)
#   - Placed / unplaced contig name lists
#   - RagTag AGP (contig → scaffold mapping)
#   - Coverage classification TSVs (placed + unplaced)
#
# Output in OUTDIR:
#   placed_genes.tsv          — all placed gene models with scaffold + coverage
#   unplaced_genes.tsv        — all unplaced gene models with coverage
#   gene_summary.tsv          — per-gene counts and tier classification
#   tier1a_1to1pair_multigene.tsv             — strongest candidates
#   tier1b_1to1pair_singlegene.tsv            — clean pair, 1 gene
#   tier1c_ambiguous_contig_pair.tsv          — gene OK but contig pair ambiguous
#   tier2_strict1to1_hemi.tsv                 — strict 1:1, any scaffold
#   tier3_hemiZ_hemiUnpl_with_paralogs.tsv    — Z-region pair but with paralogs
#   tier4_hemi_hemiUnpl_with_paralogs.tsv     — non-Z pair with paralogs
#   tier5_multiHemi.tsv                       — multiple Hemi copies
#   tier6_single_side.tsv                     — Hemi on one side only
#   contig_pair_summary.tsv                   — tier 1 contig pairs ranked
#   gametolog_report.txt                      — summary statistics
#
# Date: 2026-03-01
# Usage: bash gametolog_discovery_hanno7.sh

set -euo pipefail

# ============================================================
# PATHS — adjust as needed
# ============================================================
HANNO_BED="/data/software/HANNO/WithmRNA_HyloideaProteinSet_final_medaka_polished_assembly_consensus.fasta/BESTMODELS-FINAL.bedDB"
PLACED_NAMES="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/placed_contig_names.txt"
UNPLACED_NAMES="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/unplaced_contig_names.txt"
AGP="/data/GrenadaFrog144/ragtag.scaffold.renamed.agp"
COV_PLACED="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/coverage/placed_coverage_classified.tsv"
COV_UNPLACED="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/coverage/unplaced_coverage_classified.tsv"

OUTDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/gametolog_discovery_hanno7"

# Z-candidate region definitions (from per-Mb coverage analysis)
Z_SCF2_START=0
Z_SCF2_END=42000000
Z_SCF8_START=47000000
Z_SCF8_END=64000000

# ============================================================
# SANITY CHECKS
# ============================================================
echo "=== Gametolog Discovery — HANNO 7-Hyloidea ==="
echo "Start: $(date)"
echo ""

for f in "$HANNO_BED" "$PLACED_NAMES" "$UNPLACED_NAMES" "$AGP" \
         "$COV_PLACED" "$COV_UNPLACED"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

mkdir -p "$OUTDIR"

# ============================================================
# STEP 1: Build contig → scaffold lookup from AGP
# ============================================================
echo "[Step 1] Building contig → scaffold lookup from AGP ..."

# AGP format: scaffold \t start \t end \t part \t W/U \t contig \t c_start \t c_end \t orient
# Only W lines are contigs; U lines are gaps
# Output: contig_name \t scaffold \t scaffold_start \t scaffold_end
awk -F'\t' '
    /^#/ { next }
    $5 == "W" {
        print $6 "\t" $1 "\t" $2 "\t" $3
    }
' "$AGP" > "${OUTDIR}/agp_contig_to_scaffold.tsv"

N_AGP=$(wc -l < "${OUTDIR}/agp_contig_to_scaffold.tsv")
echo "  Contigs in AGP: ${N_AGP}"

# ============================================================
# STEP 2: Build coverage class lookups
# ============================================================
echo "[Step 2] Building coverage class lookups ..."

# placed: col1=chrom, col4=mean, col7=coverage_class (skip header)
awk -F'\t' 'NR > 1 { print $1 "\t" $4 "\t" $7 }' "$COV_PLACED" \
    > "${OUTDIR}/placed_cov_lookup.tsv"

# unplaced: col1=contig, col4=mean_cov, col7=coverage_class (skip header)
awk -F'\t' 'NR > 1 { print $1 "\t" $4 "\t" $7 }' "$COV_UNPLACED" \
    > "${OUTDIR}/unplaced_cov_lookup.tsv"

echo "  Placed contigs with coverage: $(wc -l < "${OUTDIR}/placed_cov_lookup.tsv")"
echo "  Unplaced contigs with coverage: $(wc -l < "${OUTDIR}/unplaced_cov_lookup.tsv")"

# ============================================================
# STEP 3: Extract HANNO gene models, split by placed/unplaced
# ============================================================
echo "[Step 3] Splitting HANNO gene models by placed/unplaced ..."

# Extract key fields from HANNO bedDB:
#   col1  = contig (Chrom)
#   col2  = mRNAstart (0-based)
#   col3  = mRNAend
#   col4  = Name (hanno.gXXX.tXXX)
#   col6  = strand
#   col13 = origName (TU = protein-based, MSTRG = RNA-only)
#   col26 = Preferred_name (gene symbol or "-")
#   col41 = orflen
#   col42 = mRNAlen
#
# Skip header (starts with ##)

awk -F'\t' '
    /^##/ { next }
    {
        gene_name = ($26 == "-" || $26 == "") ? "UNNAMED" : toupper($26)
        source = "unknown"
        if ($13 ~ /^TU/)    source = "protein"
        if ($13 ~ /^MSTRG/) source = "RNA-only"
        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" gene_name "\t" source "\t" $41 "\t" $42
    }
' "$HANNO_BED" > "${OUTDIR}/all_genes_raw.tsv"

N_TOTAL=$(wc -l < "${OUTDIR}/all_genes_raw.tsv")
echo "  Total gene models: ${N_TOTAL}"

# Split by placed / unplaced using exact name matching
awk -F'\t' '
    NR == FNR { placed[$1]; next }
    ($1 in placed)
' "$PLACED_NAMES" "${OUTDIR}/all_genes_raw.tsv" > "${OUTDIR}/placed_genes_raw.tsv"

awk -F'\t' '
    NR == FNR { unplaced[$1]; next }
    ($1 in unplaced)
' "$UNPLACED_NAMES" "${OUTDIR}/all_genes_raw.tsv" > "${OUTDIR}/unplaced_genes_raw.tsv"

N_PLACED=$(wc -l < "${OUTDIR}/placed_genes_raw.tsv")
N_UNPLACED=$(wc -l < "${OUTDIR}/unplaced_genes_raw.tsv")
N_NEITHER=$((N_TOTAL - N_PLACED - N_UNPLACED))

echo "  Placed gene models: ${N_PLACED}"
echo "  Unplaced gene models: ${N_UNPLACED}"
if [ "$N_NEITHER" -gt 0 ]; then
    echo "  WARNING: ${N_NEITHER} gene models on contigs not in either list"
fi

# ============================================================
# STEP 4: Add scaffold context to placed genes
# ============================================================
echo "[Step 4] Adding scaffold context and Z-region flag to placed genes ..."

# First: build a single placed contig lookup combining AGP + coverage
# contig → scaffold, scf_start, scf_end, mean_cov, cov_class, z_flag
awk -F'\t' '
    NR == FNR {
        # AGP lookup
        scf[$1] = $2; ss[$1] = $3; se[$1] = $4
        next
    }
    {
        # Coverage lookup — add scaffold info + Z-flag
        contig = $1
        scaffold = (contig in scf) ? scf[contig] : "NA"
        s_start = (contig in ss) ? ss[contig] : 0
        s_end = (contig in se) ? se[contig] : 0

        z_flag = "no"
        if (scaffold == "scaffold_2" && s_start+0 >= '"$Z_SCF2_START"' && s_end+0 <= '"$Z_SCF2_END"') {
            z_flag = "Z-candidate_scf2"
        }
        if (scaffold == "scaffold_8" && s_start+0 >= '"$Z_SCF8_START"' && s_end+0 <= '"$Z_SCF8_END"') {
            z_flag = "Z-candidate_scf8"
        }

        # contig \t mean_cov \t cov_class \t scaffold \t scf_start \t scf_end \t z_flag
        print contig "\t" $2 "\t" $3 "\t" scaffold "\t" s_start "\t" s_end "\t" z_flag
    }
' "${OUTDIR}/agp_contig_to_scaffold.tsv" \
  "${OUTDIR}/placed_cov_lookup.tsv" \
  > "${OUTDIR}/placed_contig_lookup.tsv"

# Now join placed genes with the single lookup
awk -F'\t' '
    NR == FNR {
        # Load lookup: contig → mean_cov, cov_class, scaffold, scf_start, scf_end, z_flag
        mcov[$1] = $2; ccls[$1] = $3; scf[$1] = $4
        ss[$1] = $5; se[$1] = $6; zf[$1] = $7
        next
    }
    {
        contig = $1
        print $0 "\t" \
              ((contig in scf) ? scf[contig] : "NA") "\t" \
              ((contig in ss)  ? ss[contig]  : "NA") "\t" \
              ((contig in se)  ? se[contig]  : "NA") "\t" \
              ((contig in mcov)? mcov[contig]: "NA") "\t" \
              ((contig in ccls)? ccls[contig]: "NA") "\t" \
              ((contig in zf)  ? zf[contig]  : "no")
    }
' "${OUTDIR}/placed_contig_lookup.tsv" \
  "${OUTDIR}/placed_genes_raw.tsv" \
  > "${OUTDIR}/placed_genes.tsv"

# Add header
sed -i '1i\contig\tgene_start\tgene_end\tgene_id\tstrand\tgene_name\tsource\torflen\tmRNAlen\tscaffold\tscf_start\tscf_end\tmean_cov\tcov_class\tz_flag' \
    "${OUTDIR}/placed_genes.tsv"

# Count Z-region genes
N_ZFLAG=$(awk -F'\t' 'NR > 1 && $15 != "no"' "${OUTDIR}/placed_genes.tsv" | wc -l)
echo "  Placed genes in Z-candidate regions: ${N_ZFLAG}"

# ============================================================
# STEP 5: Add coverage class to unplaced genes
# ============================================================
echo "[Step 5] Adding coverage class to unplaced genes ..."

awk -F'\t' '
    NR == FNR { cov[$1] = $2; cls[$1] = $3; next }
    {
        contig = $1
        mean_cov = (contig in cov) ? cov[contig] : "NA"
        cov_class = (contig in cls) ? cls[contig] : "NA"
        print $0 "\t" mean_cov "\t" cov_class
    }
' "${OUTDIR}/unplaced_cov_lookup.tsv" \
  "${OUTDIR}/unplaced_genes_raw.tsv" \
  > "${OUTDIR}/unplaced_genes.tsv"

# Add header
sed -i '1i\contig\tgene_start\tgene_end\tgene_id\tstrand\tgene_name\tsource\torflen\tmRNAlen\tmean_cov\tcov_class' \
    "${OUTDIR}/unplaced_genes.tsv"

# ============================================================
# STEP 6: Per-gene summary — count hits on each side by coverage class
# ============================================================
echo "[Step 6] Computing per-gene summary counts ..."

# For each named gene (excluding UNNAMED), count:
#   placed_total, placed_hemi, placed_hemi_Zregion, placed_auto
#   unplaced_total, unplaced_hemi, unplaced_auto
#
# Also collect: list of contigs and scaffolds for each gene

awk -F'\t' '
    # --- Placed genes (15 columns, skip header) ---
    NR == FNR && FNR > 1 {
        gn = $6
        if (gn == "UNNAMED") next

        p_total[gn]++
        if ($14 == "Hemi_0.5x") {
            p_hemi[gn]++
            if ($15 != "no") p_hemi_z[gn]++
            # Store placed Hemi contig details (contig:scaffold:scf_start-scf_end:z_flag)
            p_hemi_detail[gn] = p_hemi_detail[gn] (p_hemi_detail[gn] ? ";" : "") \
                $1 ":" $10 ":" $11 "-" $12 ":" $15
        }
        if ($14 == "Auto_1.0x") p_auto[gn]++

        # Track all scaffolds this gene appears on
        if (!seen_scf[gn,$10]++) {
            p_scaffolds[gn] = p_scaffolds[gn] (p_scaffolds[gn] ? "," : "") $10
        }
        next
    }

    # --- Unplaced genes (11 columns, skip header) ---
    FNR > 1 {
        gn = $6
        if (gn == "UNNAMED") next

        u_total[gn]++
        if ($11 == "Hemi_0.5x") {
            u_hemi[gn]++
            # Store unplaced Hemi contig details (contig:mean_cov)
            u_hemi_detail[gn] = u_hemi_detail[gn] (u_hemi_detail[gn] ? ";" : "") \
                $1 ":" $10
        }
        if ($11 == "Auto_1.0x") u_auto[gn]++
    }

    END {
        # Header
        print "gene_name\t" \
              "placed_total\tplaced_hemi\tplaced_hemi_Zregion\tplaced_auto\t" \
              "unplaced_total\tunplaced_hemi\tunplaced_auto\t" \
              "placed_scaffolds\tplaced_hemi_details\tunplaced_hemi_details\t" \
              "tier"

        # Iterate over all genes seen on either side
        for (gn in p_total) { genes[gn] }
        for (gn in u_total) { genes[gn] }

        for (gn in genes) {
            pt = p_total[gn]+0
            ph = p_hemi[gn]+0
            phz = p_hemi_z[gn]+0
            pa = p_auto[gn]+0
            ut = u_total[gn]+0
            uh = u_hemi[gn]+0
            ua = u_auto[gn]+0

            # --- Tier classification ---
            tier = "unclassified"

            if (ph > 0 && uh > 0) {
                # Both sides have Hemi hits
                if (ph == 1 && uh == 1 && pt == 1 && ut == 1 && phz == 1) {
                    tier = "tier1_strict1to1_hemiZ"
                } else if (ph == 1 && uh == 1 && pt == 1 && ut == 1 && phz == 0) {
                    tier = "tier2_strict1to1_hemi"
                } else if (ph == 1 && uh == 1 && phz == 1 && (pt > 1 || ut > 1)) {
                    tier = "tier3_hemiZ_hemiUnpl_with_paralogs"
                } else if (ph == 1 && uh == 1 && phz == 0 && (pt > 1 || ut > 1)) {
                    tier = "tier4_hemi_hemiUnpl_with_paralogs"
                } else {
                    tier = "tier5_multiHemi"
                }
            } else if (ph > 0 && uh == 0) {
                # Hemi on placed side only
                if (ut > 0) {
                    tier = "tier6_placedHemi_unplacedNotHemi"
                } else if (phz > 0) {
                    tier = "tier6_placedHemiZ_only"
                } else {
                    tier = "tier6_placedHemi_only"
                }
            } else if (ph == 0 && uh > 0) {
                # Hemi on unplaced side only
                if (pt > 0) {
                    tier = "tier6_unplacedHemi_placedNotHemi"
                } else {
                    tier = "tier6_unplacedHemi_only"
                }
            } else {
                # No Hemi on either side
                if (pt > 0 && ut > 0) {
                    tier = "tier7_both_noHemi"
                } else if (pt > 0) {
                    tier = "tier7_placedOnly_noHemi"
                } else {
                    tier = "tier7_unplacedOnly_noHemi"
                }
            }

            phd = (p_hemi_detail[gn] != "") ? p_hemi_detail[gn] : "-"
            uhd = (u_hemi_detail[gn] != "") ? u_hemi_detail[gn] : "-"
            psc = (p_scaffolds[gn] != "") ? p_scaffolds[gn] : "-"

            print gn "\t" pt "\t" ph "\t" phz "\t" pa "\t" \
                  ut "\t" uh "\t" ua "\t" \
                  psc "\t" phd "\t" uhd "\t" tier
        }
    }
' "${OUTDIR}/placed_genes.tsv" "${OUTDIR}/unplaced_genes.tsv" \
  > "${OUTDIR}/gene_summary_unsorted.tsv"

# Sort: tier1 first, then tier2, etc., then alphabetically within tier
head -1 "${OUTDIR}/gene_summary_unsorted.tsv" > "${OUTDIR}/gene_summary.tsv"
tail -n +2 "${OUTDIR}/gene_summary_unsorted.tsv" | sort -t$'\t' -k12,12 -k1,1 \
    >> "${OUTDIR}/gene_summary.tsv"

rm "${OUTDIR}/gene_summary_unsorted.tsv"

# ============================================================
# STEP 6b: Refine tier 1 by contig-pair analysis
# ============================================================
echo "[Step 6b] Refining tier 1 by contig-pair analysis ..."

# For tier1 genes: extract Z-contig and W-contig names, check 1:1
# correspondence, count genes per pair, reclassify.
#
# Tier 1a: clean 1:1 contig pair (Z-contig→1 W-contig AND W-contig→1 Z-contig),
#          ≥2 tier1 genes on this pair
# Tier 1b: clean 1:1 contig pair, exactly 1 tier1 gene
# Tier 1c: ambiguous — Z-contig pairs with multiple W-contigs or vice versa

awk -F'\t' '
    NR == 1 { header = $0; next }

    # Pass 1: collect tier1 genes and their contig pairs
    $12 == "tier1_strict1to1_hemiZ" {
        # Extract Z-contig name (before first ":" in placed_hemi_details col10)
        split($10, zparts, ":")
        z_contig = zparts[1]

        # Extract W-contig name (before first ":" in unplaced_hemi_details col11)
        split($11, wparts, ":")
        w_contig = wparts[1]

        gene_z[NR] = z_contig
        gene_w[NR] = w_contig
        line[NR] = $0
        tier1_lines[++t1_count] = NR

        # Track which W-contigs each Z-contig is associated with (via gene names)
        if (!seen_zw[z_contig, w_contig]++) {
            z_to_w_count[z_contig]++
            z_to_w_list[z_contig] = z_to_w_list[z_contig] (z_to_w_list[z_contig] ? "," : "") w_contig
        }

        # Track which Z-contigs each W-contig is associated with (via gene names)
        if (!seen_wz[w_contig, z_contig]++) {
            w_to_z_count[w_contig]++
            w_to_z_list[w_contig] = w_to_z_list[w_contig] (w_to_z_list[w_contig] ? "," : "") z_contig
        }

        # Count genes per Z:W pair
        pair = z_contig "\t" w_contig
        pair_gene_count[pair]++
        pair_genes[pair] = pair_genes[pair] (pair_genes[pair] ? "," : "") $1
    }

    # Non-tier1: just store
    $12 != "tier1_strict1to1_hemiZ" {
        other_lines[++other_count] = $0
    }

    END {
        print header

        # Reclassify each tier1 gene
        for (i = 1; i <= t1_count; i++) {
            nr = tier1_lines[i]
            zc = gene_z[nr]
            wc = gene_w[nr]
            pair = zc "\t" wc

            # Check 1:1 in both directions
            z_clean = (z_to_w_count[zc] == 1) ? 1 : 0
            w_clean = (w_to_z_count[wc] == 1) ? 1 : 0

            if (z_clean && w_clean) {
                # Clean 1:1 contig pair
                n_genes = pair_gene_count[pair]
                if (n_genes >= 2) {
                    new_tier = "tier1a_1to1pair_multigene"
                } else {
                    new_tier = "tier1b_1to1pair_singlegene"
                }
            } else {
                new_tier = "tier1c_ambiguous_contig_pair"
            }

            # Replace tier in line
            # Line has 12 tab-separated fields, tier is field 12
            n = split(line[nr], fields, "\t")
            fields[12] = new_tier
            out = fields[1]
            for (f = 2; f <= n; f++) out = out "\t" fields[f]
            print out
        }

        # Print all non-tier1 lines
        for (i = 1; i <= other_count; i++) {
            print other_lines[i]
        }
    }
' "${OUTDIR}/gene_summary.tsv" > "${OUTDIR}/gene_summary_refined.tsv"

# Re-sort
head -1 "${OUTDIR}/gene_summary_refined.tsv" > "${OUTDIR}/gene_summary.tsv"
tail -n +2 "${OUTDIR}/gene_summary_refined.tsv" | sort -t$'\t' -k12,12 -k1,1 \
    >> "${OUTDIR}/gene_summary.tsv"
rm "${OUTDIR}/gene_summary_refined.tsv"

# Count subtiers
N_1A=$(awk -F'\t' '$12 == "tier1a_1to1pair_multigene"' "${OUTDIR}/gene_summary.tsv" | wc -l)
N_1B=$(awk -F'\t' '$12 == "tier1b_1to1pair_singlegene"' "${OUTDIR}/gene_summary.tsv" | wc -l)
N_1C=$(awk -F'\t' '$12 == "tier1c_ambiguous_contig_pair"' "${OUTDIR}/gene_summary.tsv" | wc -l)
echo "  Tier 1a (clean pair, multi-gene):  ${N_1A}"
echo "  Tier 1b (clean pair, single gene): ${N_1B}"
echo "  Tier 1c (ambiguous contig pair):   ${N_1C}"

# Generate contig pair summary table
echo ""
echo "[Step 6c] Generating contig pair summary ..."

awk -F'\t' '
    $12 ~ /^tier1[abc]/ {
        split($10, zp, ":")
        split($11, wp, ":")
        z = zp[1]; w = wp[1]
        pair = z "\t" w
        pair_tier[$1] = $12
        pair_z_detail[$1] = $10
        pair_w_detail[$1] = $11

        genes[pair] = genes[pair] (genes[pair] ? "," : "") $1
        count[pair]++
        if (!z_detail[pair]) z_detail[pair] = $10
        if (!w_detail[pair]) w_detail[pair] = $11

        # Get the subtier — all genes on same pair have same subtier
        subtier[pair] = $12
    }
    END {
        print "n_genes\ttier\tZ_contig_detail\tW_contig_detail\tgene_names"
        for (p in count) {
            print count[p] "\t" subtier[p] "\t" z_detail[p] "\t" w_detail[p] "\t" genes[p]
        }
    }
' "${OUTDIR}/gene_summary.tsv" > "${OUTDIR}/contig_pair_summary_unsorted.tsv"

# Sort data lines only, keep header on top
head -1 "${OUTDIR}/contig_pair_summary_unsorted.tsv" > "${OUTDIR}/contig_pair_summary.tsv"
tail -n +2 "${OUTDIR}/contig_pair_summary_unsorted.tsv" | sort -t$'\t' -k1,1rn -k2,2 \
    >> "${OUTDIR}/contig_pair_summary.tsv"
rm "${OUTDIR}/contig_pair_summary_unsorted.tsv"

N_PAIRS=$(awk -F'\t' 'NR > 1' "${OUTDIR}/contig_pair_summary.tsv" | wc -l)
N_CLEAN=$(awk -F'\t' 'NR > 1 && $2 ~ /tier1[ab]/' "${OUTDIR}/contig_pair_summary.tsv" | wc -l)
echo "  Total tier 1 contig pairs: ${N_PAIRS}"
echo "  Clean 1:1 pairs: ${N_CLEAN}"

# ============================================================
# STEP 7: Extract tier-specific tables
# ============================================================
echo "[Step 7] Extracting tier-specific tables ..."

for TIER in tier1a_1to1pair_multigene tier1b_1to1pair_singlegene tier1c_ambiguous_contig_pair tier2_strict1to1_hemi tier3_hemiZ_hemiUnpl_with_paralogs tier4_hemi_hemiUnpl_with_paralogs tier5_multiHemi; do
    head -1 "${OUTDIR}/gene_summary.tsv" > "${OUTDIR}/${TIER}.tsv"
    awk -F'\t' -v t="$TIER" '$12 == t' "${OUTDIR}/gene_summary.tsv" \
        >> "${OUTDIR}/${TIER}.tsv"
done

# All tier6 variants in one file
head -1 "${OUTDIR}/gene_summary.tsv" > "${OUTDIR}/tier6_single_side.tsv"
awk -F'\t' '$12 ~ /^tier6/' "${OUTDIR}/gene_summary.tsv" \
    >> "${OUTDIR}/tier6_single_side.tsv"

# ============================================================
# STEP 8: Summary report
# ============================================================
echo "[Step 8] Generating summary report ..."

{
    echo "======================================================================"
    echo "GAMETOLOG DISCOVERY REPORT — HANNO 7-Hyloidea on full primary assembly"
    echo "======================================================================"
    echo "Date: $(date)"
    echo ""
    echo "Input:"
    echo "  HANNO: ${HANNO_BED}"
    echo "  AGP:   ${AGP}"
    echo ""
    echo "--- Gene model counts ---"
    echo ""
    echo "Total gene models in HANNO:      ${N_TOTAL}"
    echo "  On placed contigs:             ${N_PLACED}"
    echo "  On unplaced contigs:           ${N_UNPLACED}"
    if [ "$N_NEITHER" -gt 0 ]; then
        echo "  On neither (WARNING):          ${N_NEITHER}"
    fi
    echo ""

    # Named vs unnamed
    N_NAMED_P=$(awk -F'\t' 'NR > 1 && $6 != "UNNAMED"' "${OUTDIR}/placed_genes.tsv" | wc -l)
    N_UNNAMED_P=$(awk -F'\t' 'NR > 1 && $6 == "UNNAMED"' "${OUTDIR}/placed_genes.tsv" | wc -l)
    N_NAMED_U=$(awk -F'\t' 'NR > 1 && $6 != "UNNAMED"' "${OUTDIR}/unplaced_genes.tsv" | wc -l)
    N_UNNAMED_U=$(awk -F'\t' 'NR > 1 && $6 == "UNNAMED"' "${OUTDIR}/unplaced_genes.tsv" | wc -l)

    echo "Named gene models:"
    echo "  Placed:    ${N_NAMED_P} named, ${N_UNNAMED_P} unnamed"
    echo "  Unplaced:  ${N_NAMED_U} named, ${N_UNNAMED_U} unnamed"
    echo ""

    # Coverage class breakdown for placed
    echo "--- Placed genes by coverage class ---"
    awk -F'\t' 'NR > 1 { print $14 }' "${OUTDIR}/placed_genes.tsv" | \
        sort | uniq -c | sort -rn | awk '{ printf "  %-25s %d\n", $2, $1 }'
    echo ""

    # Coverage class breakdown for unplaced
    echo "--- Unplaced genes by coverage class ---"
    awk -F'\t' 'NR > 1 { print $11 }' "${OUTDIR}/unplaced_genes.tsv" | \
        sort | uniq -c | sort -rn | awk '{ printf "  %-25s %d\n", $2, $1 }'
    echo ""

    # Z-region breakdown
    echo "--- Placed genes by Z-region flag ---"
    awk -F'\t' 'NR > 1 { print $15 }' "${OUTDIR}/placed_genes.tsv" | \
        sort | uniq -c | sort -rn | awk '{ printf "  %-25s %d\n", $2, $1 }'
    echo ""

    # Tier counts (unique genes, not gene models)
    echo "--- Gametolog tier classification (unique named genes) ---"
    echo ""
    awk -F'\t' 'NR > 1 { print $12 }' "${OUTDIR}/gene_summary.tsv" | \
        sort | uniq -c | sort -rn | \
        awk '{ printf "  %-40s %d\n", $2, $1 }'
    echo ""

    # Tier 1a gene list
    N_T1A=$(awk -F'\t' 'NR > 1' "${OUTDIR}/tier1a_1to1pair_multigene.tsv" | wc -l)
    echo "--- Tier 1a gametolog candidates (n=${N_T1A}) ---"
    echo "  (strict 1:1, Z-region, clean 1:1 contig pair, >=2 genes on pair)"
    echo ""
    if [ "$N_T1A" -gt 0 ]; then
        awk -F'\t' 'NR > 1 { printf "  %-20s Z: %-55s  W: %s\n", \
            $1, $10, $11 }' "${OUTDIR}/tier1a_1to1pair_multigene.tsv"
    else
        echo "  (none)"
    fi
    echo ""

    # Tier 1b gene list
    N_T1B=$(awk -F'\t' 'NR > 1' "${OUTDIR}/tier1b_1to1pair_singlegene.tsv" | wc -l)
    echo "--- Tier 1b gametolog candidates (n=${N_T1B}) ---"
    echo "  (strict 1:1, Z-region, clean 1:1 contig pair, 1 gene on pair)"
    echo ""
    if [ "$N_T1B" -gt 0 ]; then
        awk -F'\t' 'NR > 1 { printf "  %-20s Z: %-55s  W: %s\n", \
            $1, $10, $11 }' "${OUTDIR}/tier1b_1to1pair_singlegene.tsv"
    else
        echo "  (none)"
    fi
    echo ""

    # Tier 1c gene list
    N_T1C=$(awk -F'\t' 'NR > 1' "${OUTDIR}/tier1c_ambiguous_contig_pair.tsv" | wc -l)
    echo "--- Tier 1c gametolog candidates (n=${N_T1C}) ---"
    echo "  (strict 1:1 gene-level, Z-region, BUT Z- or W-contig associated with multiple contigs on other side via different genes)"
    echo ""
    if [ "$N_T1C" -gt 0 ]; then
        awk -F'\t' 'NR > 1 { printf "  %-20s Z: %-55s  W: %s\n", \
            $1, $10, $11 }' "${OUTDIR}/tier1c_ambiguous_contig_pair.tsv"
    else
        echo "  (none)"
    fi
    echo ""

    # Contig pair summary
    echo "--- Contig pair summary (tier 1a/1b/1c) ---"
    echo ""
    awk -F'\t' 'NR > 1 { printf "  %d genes  %-35s  Z: %-55s  W: %s\n  %*sgenes: %s\n\n", \
        $1, $2, $3, $4, 10, "", $5 }' "${OUTDIR}/contig_pair_summary.tsv"
    echo ""

    # Tier 2 gene list
    N_T2=$(awk -F'\t' 'NR > 1' "${OUTDIR}/tier2_strict1to1_hemi.tsv" | wc -l)
    echo "--- Tier 2 gametolog candidates (n=${N_T2}) ---"
    echo "  (strictly 1 total placed [Hemi, any scaffold] + 1 total unplaced [Hemi], no paralogs)"
    echo ""
    if [ "$N_T2" -gt 0 ]; then
        awk -F'\t' 'NR > 1 { printf "  %-20s placed_scf: %-12s  W-contig: %s\n", \
            $1, $9, $11 }' "${OUTDIR}/tier2_strict1to1_hemi.tsv"
    else
        echo "  (none)"
    fi
    echo ""

    echo "======================================================================"
    echo "Output files in: ${OUTDIR}/"
    echo "======================================================================"
    echo "  placed_genes.tsv                          — all placed gene models (15 cols)"
    echo "  unplaced_genes.tsv                        — all unplaced gene models (11 cols)"
    echo "  gene_summary.tsv                          — per-gene counts + tier (12 cols)"
    echo "  contig_pair_summary.tsv                   — tier 1 contig pairs ranked by gene count"
    echo "  tier1a_1to1pair_multigene.tsv             — strongest (clean pair, >=2 genes)"
    echo "  tier1b_1to1pair_singlegene.tsv            — clean pair, 1 gene"
    echo "  tier1c_ambiguous_contig_pair.tsv          — gene OK but contig pair ambiguous"
    echo "  tier2_strict1to1_hemi.tsv                 — strict 1:1, any scaffold"
    echo "  tier3_hemiZ_hemiUnpl_with_paralogs.tsv    — Z-region pair but with paralogs"
    echo "  tier4_hemi_hemiUnpl_with_paralogs.tsv     — non-Z pair with paralogs"
    echo "  tier5_multiHemi.tsv                       — multiple Hemi copies"
    echo "  tier6_single_side.tsv                     — Hemi on one side only"
    echo "  gametolog_report.txt                      — this report"
    echo ""
    echo "Done: $(date)"

} > "${OUTDIR}/gametolog_report.txt"

# Print report to stdout as well
cat "${OUTDIR}/gametolog_report.txt"

# ============================================================
# CLEANUP intermediate files
# ============================================================
rm -f "${OUTDIR}/all_genes_raw.tsv" \
      "${OUTDIR}/placed_genes_raw.tsv" \
      "${OUTDIR}/unplaced_genes_raw.tsv" \
      "${OUTDIR}/agp_contig_to_scaffold.tsv" \
      "${OUTDIR}/placed_cov_lookup.tsv" \
      "${OUTDIR}/unplaced_cov_lookup.tsv" \
      "${OUTDIR}/placed_contig_lookup.tsv"

echo ""
echo "=== Gametolog Discovery complete ==="
