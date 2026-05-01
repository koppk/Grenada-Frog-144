#!/bin/bash
#
# parse_repeatmasker_two_sets.sh
# ================================
# Parse RepeatMasker .out files for placed and unplaced contig sets.
# Produces per-contig TE class breakdowns, per-set summaries, and
# cross-tabulation with coverage classification.
#
# Input:
#   RepeatMasker .out files (standard format, 3 header lines)
#   Coverage classification TSV (from map_reads_primary_assembly_coverage.sh)
#
# Output in ${OUTDIR}:
#   placed_TE_per_contig.tsv          : per-contig TE bp by class (placed)
#   unplaced_TE_per_contig.tsv        : per-contig TE bp by class (unplaced)
#   TE_summary_by_set.tsv             : aggregate TE proportions per set
#   unplaced_TE_by_coverage_class.tsv : TE profiles cross-tabulated with coverage class
#   placed_TE_by_scaffold_zone.tsv    : Z-candidate vs autosomal TE profiles (via AGP)
#   z_candidate_contig_names.txt      : contig names in Z-candidate regions
#   TE_two_set_report.txt             : human-readable summary report
#
# Usage:
#   bash parse_repeatmasker_two_sets.sh
#

# Author: Kopp K. Pristimantis euphronides genome project.
set -euo pipefail

# === Paths ===
BASEDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"
RMDIR="${BASEDIR}/RepeatMasker_output"
OUTDIR="${BASEDIR}/RepeatMasker_parsed"
COVDIR="${BASEDIR}/coverage"

PLACED_OUT="${RMDIR}/placed_contigs/placed_contigs.fasta.out"
UNPLACED_OUT="${RMDIR}/unplaced_contigs_filtered/unplaced_contigs_filtered.fasta.out"

# Coverage classification from map_reads_primary_assembly_coverage.sh
# Contains ALL unplaced contigs (25224 with coverage + 67 zero-coverage = 25291)
COV_CLASS="${COVDIR}/unplaced_coverage_classified.tsv"

# Contig lengths (from mosdepth summary)
# placed_coverage.tsv: 8815 contigs (all placed, no length filter)
# Header: chrom  length  bases  mean  min  max
PLACED_COV="${COVDIR}/placed_coverage.tsv"

# RagTag AGP for contig-to-scaffold mapping (Z-candidate region analysis)
AGP="/data/GrenadaFrog144/coverage/ragtag.scaffold.renamed.agp"

echo "=== Parse RepeatMasker results ==="
echo "Start: $(date)"
echo ""

mkdir -p "$OUTDIR"

# ============================================================
# Reusable function: parse a RepeatMasker .out file into
# per-contig TE class breakdown
#
# Arguments:
#   $1 RMFILE   - RepeatMasker .out file
#   $2 LABEL    - set name (placed/unplaced)
#   $3 LENFILE  - TSV with contig name in col 1, length in col 2 (has header)
#   $4 OUTFILE  - output TSV path
#   $5 MIN_LEN  - minimum contig length filter (0 = no filter)
#
# Output columns (13):
#   contig  total_masked_bp  SINE  LINE  LTR  DNA  RC_Helitron
#   Satellite  Simple_Low  Unknown  Other  contig_length  masked_pct
#
# LEFT JOIN from LENFILE → TE data: every contig in LENFILE
# (passing MIN_LEN filter) gets a row. Contigs without TE hits
# get all-zero columns. This ensures correct denominators.
# ============================================================
parse_rm_out() {
    local RMFILE="$1"
    local LABEL="$2"
    local LENFILE="$3"
    local OUTFILE="$4"
    local MIN_LEN="${5:-0}"

    if [ ! -f "$RMFILE" ]; then
        echo "  WARNING: $RMFILE not found, skipping $LABEL"
        return 1
    fi

    echo "  Parsing $LABEL: $RMFILE"

    # RepeatMasker .out has 3 header lines, then space-delimited data
    # Key columns (1-indexed after cleanup):
    #   5  = query sequence (contig name)
    #   6  = begin in query
    #   7  = end in query
    #   11 = repeat class/family (e.g. "LINE/L1", "SINE/tRNA", "DNA/hAT")
    #
    # We extract: contig, masked_bp (end - begin + 1), major_class
    # Major class = part before "/" in class/family

    awk '
    NR <= 3 { next }                       # skip 3 header lines
    /^\s*$/ { next }                       # skip blank lines
    {
        # Collapse whitespace, fields are space-separated
        gsub(/^ +/, "")
        n = split($0, f, /[ \t]+/)
        if (n < 11) next

        contig = f[5]
        qbegin = f[6] + 0
        qend   = f[7] + 0
        masked_bp = qend - qbegin + 1
        if (masked_bp < 1) masked_bp = 1

        class_family = f[11]

        # Extract major class (before /)
        split(class_family, cf, "/")
        major = cf[1]

        # Normalize major classes
        if (major == "SINE")            cls = "SINE"
        else if (major == "LINE")       cls = "LINE"
        else if (major == "LTR")        cls = "LTR"
        else if (major == "DNA")        cls = "DNA"
        else if (major == "RC")         cls = "RC_Helitron"
        else if (major == "Satellite")  cls = "Satellite"
        else if (major == "Simple_repeat" || major == "Low_complexity") cls = "Simple_Low"
        else if (major == "Unknown")    cls = "Unknown"
        else                            cls = "Other"

        bp[contig][cls] += masked_bp
        total[contig] += masked_bp
    }
    END {
        # Header
        printf "contig\ttotal_masked_bp\tSINE\tLINE\tLTR\tDNA\tRC_Helitron\tSatellite\tSimple_Low\tUnknown\tOther\n"
        for (c in total) {
            printf "%s\t%d", c, total[c]
            split("SINE,LINE,LTR,DNA,RC_Helitron,Satellite,Simple_Low,Unknown,Other", classes, ",")
            for (i = 1; i <= 9; i++) {
                v = (bp[c][classes[i]] > 0) ? bp[c][classes[i]] : 0
                printf "\t%d", v
            }
            printf "\n"
        }
    }' "$RMFILE" > "${OUTDIR}/${LABEL}_TE_raw.tmp"

    # LEFT JOIN: read TE data first (file 1), then iterate LENFILE (file 2).
    # Every contig in LENFILE that passes MIN_LEN gets a row.
    # Contigs without TE hits get all-zero columns.
    #
    # File 1: TE raw tmp (has header)
    #   col 1 = contig, col 2 = total_masked_bp, cols 3-11 = TE classes
    # File 2: LENFILE (has header)
    #   col 1 = contig, col 2 = length
    awk -F'\t' -v min_len="$MIN_LEN" '
    NR == FNR {
        # File 1: TE raw tmp — store entire line keyed by contig name
        if (FNR > 1) te_line[$1] = $0
        next
    }
    # File 2: LENFILE
    FNR == 1 { next }
    {
        contig = $1
        clen = $2 + 0
        if (clen < min_len) next

        if (te_line[contig] != "") {
            # Has TE data
            split(te_line[contig], f, "\t")
            total = f[2] + 0
            printf "%s", te_line[contig]
        } else {
            # No TE hits — all zeros
            printf "%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0", contig
            total = 0
        }
        pct = (clen > 0) ? total / clen * 100 : 0
        printf "\t%d\t%.4f\n", clen, pct
    }
    BEGIN {
        printf "contig\ttotal_masked_bp\tSINE\tLINE\tLTR\tDNA\tRC_Helitron\tSatellite\tSimple_Low\tUnknown\tOther\tcontig_length\tmasked_pct\n"
    }' "${OUTDIR}/${LABEL}_TE_raw.tmp" "$LENFILE" > "$OUTFILE"

    rm -f "${OUTDIR}/${LABEL}_TE_raw.tmp"

    local N_TOTAL=$(tail -n+2 "$OUTFILE" | wc -l)
    local N_WITH_TE=$(awk -F'\t' 'NR>1 && $2>0' "$OUTFILE" | wc -l)
    local N_NO_TE=$(awk -F'\t' 'NR>1 && $2==0' "$OUTFILE" | wc -l)
    echo "    Total contigs:    $N_TOTAL"
    echo "    With TE hits:     $N_WITH_TE"
    echo "    No TE hits:       $N_NO_TE"
    echo "    Output: $OUTFILE"
    return 0
}

# ============================================================
# Step 1: Parse each set
# ============================================================
echo "[Step 1] Parsing RepeatMasker output per contig ..."
echo ""

PLACED_OK=false
UNPLACED_OK=false

# Placed: all placed contigs, no length filter (placed_contigs.fasta
# included all placed contigs regardless of size)
if parse_rm_out "$PLACED_OUT" "placed" "$PLACED_COV" "${OUTDIR}/placed_TE_per_contig.tsv" 0; then
    PLACED_OK=true
fi

# Unplaced: COV_CLASS has all 25291 unplaced contigs (25224 with coverage
# + 67 zero-coverage contigs). RepeatMasker input was >=500bp filtered
# (25079 contigs), so apply min_len=500 to match.
if parse_rm_out "$UNPLACED_OUT" "unplaced" "$COV_CLASS" "${OUTDIR}/unplaced_TE_per_contig.tsv" 500; then
    UNPLACED_OK=true
fi

echo ""

# ============================================================
# Step 2: Per-set summary (aggregate TE proportions)
# ============================================================
echo "[Step 2] Computing per-set TE summaries ..."

summarize_set() {
    local TEFILE="$1"
    local LABEL="$2"

    if [ ! -f "$TEFILE" ]; then
        echo "  Skipping $LABEL (file not found)"
        return
    fi

    # Output: 13 columns
    # col 2 = total_masked_bp, cols 3-11 = TE classes, col 12 = contig_length
    awk -F'\t' -v label="$LABEL" '
    NR > 1 {
        n++
        total_masked += $2
        total_len += $12
        for (i = 3; i <= 11; i++) sums[i] += $i
    }
    END {
        pct = (total_len > 0) ? total_masked / total_len * 100 : 0
        printf "%s\t%d\t%d\t%.2f", label, n, total_len, pct
        for (i = 3; i <= 11; i++) {
            cpct = (total_len > 0) ? sums[i] / total_len * 100 : 0
            printf "\t%.2f", cpct
        }
        printf "\n"
    }' "$TEFILE"
}

{
    echo -e "set\tn_contigs\ttotal_bp\tmasked_pct\tSINE_pct\tLINE_pct\tLTR_pct\tDNA_pct\tRC_Helitron_pct\tSatellite_pct\tSimple_Low_pct\tUnknown_pct\tOther_pct"
    if $PLACED_OK; then
        summarize_set "${OUTDIR}/placed_TE_per_contig.tsv" "placed"
    fi
    if $UNPLACED_OK; then
        summarize_set "${OUTDIR}/unplaced_TE_per_contig.tsv" "unplaced"
    fi
} > "${OUTDIR}/TE_summary_by_set.tsv"

echo "  Output: ${OUTDIR}/TE_summary_by_set.tsv"
echo ""

# ============================================================
# Step 3: Cross-tabulate unplaced TE profiles with coverage class
# ============================================================
echo "[Step 3] Cross-tabulating unplaced TE with coverage classes ..."

if $UNPLACED_OK && [ -f "$COV_CLASS" ]; then
    # Join: coverage classification (col 1=contig, col 7=class) with TE per-contig
    # Output: 13 columns (same as per-contig TSV)
    # col 2 = total_masked_bp, cols 3-11 = TE classes, col 12 = contig_length
    awk -F'\t' '
    NR == FNR && FNR > 1 {
        covclass[$1] = $7
        next
    }
    FNR == 1 { next }
    {
        contig = $1
        cls = (covclass[contig] != "") ? covclass[contig] : "unclassified"
        n[cls]++
        total_len[cls] += $12
        total_masked[cls] += $2
        for (i = 3; i <= 11; i++) sums[cls, i] += $i
    }
    END {
        printf "coverage_class\tn_contigs\ttotal_bp\tmasked_pct\tSINE_pct\tLINE_pct\tLTR_pct\tDNA_pct\tRC_Helitron_pct\tSatellite_pct\tSimple_Low_pct\tUnknown_pct\tOther_pct\n"
        split("Hemi_0.5x,Auto_1.0x,High_Coverage/Repeat,Low_Coverage/Other,Zero,unclassified", order, ",")
        for (o = 1; o <= 6; o++) {
            c = order[o]
            if (n[c] > 0) {
                pct = (total_len[c] > 0) ? total_masked[c] / total_len[c] * 100 : 0
                printf "%s\t%d\t%d\t%.2f", c, n[c], total_len[c], pct
                for (i = 3; i <= 11; i++) {
                    cpct = (total_len[c] > 0) ? sums[c, i] / total_len[c] * 100 : 0
                    printf "\t%.2f", cpct
                }
                printf "\n"
            }
        }
    }' "$COV_CLASS" "${OUTDIR}/unplaced_TE_per_contig.tsv" \
        > "${OUTDIR}/unplaced_TE_by_coverage_class.tsv"

    echo "  Output: ${OUTDIR}/unplaced_TE_by_coverage_class.tsv"
else
    echo "  Skipping (unplaced TE or coverage classification not available)"
fi
echo ""

# ============================================================
# Step 3b: Split placed contigs into Z-candidate vs autosomal
#          using RagTag AGP scaffold coordinates
# ============================================================
echo "[Step 3b] Splitting placed contigs by Z-candidate vs autosomal ..."

if $PLACED_OK && [ -f "$AGP" ]; then
    # Extract contig names in Z-candidate regions from AGP
    # Z-candidate: scaffold_2 0-42 Mb, scaffold_8 47-64 Mb
    awk -F'\t' '$5 == "W" {
        if ($1 == "scaffold_2" && $2 < 42000000)
            print $6
        else if ($1 == "scaffold_8" && $2 >= 47000000 && $3 <= 64000000)
            print $6
    }' "$AGP" | sort -u > "${OUTDIR}/z_candidate_contig_names.txt"

    N_Z=$(wc -l < "${OUTDIR}/z_candidate_contig_names.txt")
    echo "  Z-candidate contigs (scf2 0-42Mb + scf8 47-64Mb): $N_Z"

    # Classify placed contigs and compute TE profiles per zone
    # col 2 = total_masked_bp, cols 3-11 = TE classes, col 12 = contig_length
    awk -F'\t' '
    NR == FNR { zcontigs[$1]; next }
    FNR == 1 { next }
    {
        contig = $1
        cls = (contig in zcontigs) ? "Z-candidate" : "Autosomal"
        n[cls]++
        total_len[cls] += $12
        total_masked[cls] += $2
        for (i = 3; i <= 11; i++) sums[cls, i] += $i
    }
    END {
        printf "scaffold_zone\tn_contigs\ttotal_bp\tmasked_pct\tSINE_pct\tLINE_pct\tLTR_pct\tDNA_pct\tRC_Helitron_pct\tSatellite_pct\tSimple_Low_pct\tUnknown_pct\tOther_pct\n"
        split("Z-candidate,Autosomal", order, ",")
        for (o = 1; o <= 2; o++) {
            c = order[o]
            if (n[c] > 0) {
                pct = (total_len[c] > 0) ? total_masked[c] / total_len[c] * 100 : 0
                printf "%s\t%d\t%d\t%.2f", c, n[c], total_len[c], pct
                for (i = 3; i <= 11; i++) {
                    cpct = (total_len[c] > 0) ? sums[c, i] / total_len[c] * 100 : 0
                    printf "\t%.2f", cpct
                }
                printf "\n"
            }
        }
    }' "${OUTDIR}/z_candidate_contig_names.txt" \
       "${OUTDIR}/placed_TE_per_contig.tsv" \
       > "${OUTDIR}/placed_TE_by_scaffold_zone.tsv"

    echo "  Output: ${OUTDIR}/placed_TE_by_scaffold_zone.tsv"
else
    echo "  Skipping (placed TE or AGP not available)"
fi
echo ""

# ============================================================
# Step 4: Human-readable report
# ============================================================
echo "[Step 4] Generating summary report ..."

REPORT="${OUTDIR}/TE_two_set_report.txt"

{
    echo "RepeatMasker Analysis Report"
    echo "============================"
    echo "Date: $(date)"
    echo "RepeatMasker output: ${RMDIR}"
    echo ""

    echo "=== Per-set TE summary ==="
    echo ""
    if [ -f "${OUTDIR}/TE_summary_by_set.tsv" ]; then
        column -t -s$'\t' "${OUTDIR}/TE_summary_by_set.tsv"
    fi
    echo ""

    if [ -f "${OUTDIR}/unplaced_TE_by_coverage_class.tsv" ]; then
        echo "=== Unplaced contigs: TE profile by coverage class ==="
        echo ""
        column -t -s$'\t' "${OUTDIR}/unplaced_TE_by_coverage_class.tsv"
        echo ""
    fi

    if [ -f "${OUTDIR}/placed_TE_by_scaffold_zone.tsv" ]; then
        echo "=== Placed contigs: TE profile by scaffold zone ==="
        echo "  Z-candidate: scaffold_2 (0-42 Mb) + scaffold_8 (47-64 Mb)"
        echo "  Autosomal: all other placed contigs"
        echo ""
        column -t -s$'\t' "${OUTDIR}/placed_TE_by_scaffold_zone.tsv"
        echo ""
    fi

    echo "=== Output files ==="
    echo "  ${OUTDIR}/placed_TE_per_contig.tsv"
    echo "  ${OUTDIR}/unplaced_TE_per_contig.tsv"
    echo "  ${OUTDIR}/TE_summary_by_set.tsv"
    echo "  ${OUTDIR}/unplaced_TE_by_coverage_class.tsv"
    echo "  ${OUTDIR}/placed_TE_by_scaffold_zone.tsv"
    echo "  ${OUTDIR}/z_candidate_contig_names.txt"
    echo "  ${OUTDIR}/TE_two_set_report.txt"
} 2>&1 | tee "$REPORT"

echo ""
echo "Finished: $(date)"
echo "=== Done ==="
