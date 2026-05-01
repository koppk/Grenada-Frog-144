#!/bin/bash
#
# map_reads_primary_assembly_coverage.sh
# =======================================
# Map all ONT reads to the unscaffolded primary assembly and compute
# per-contig coverage using mosdepth for placed vs unplaced comparison.
#
# Purpose: Coverage-based characterisation of placed vs unplaced contigs
# - Compare coverage distributions between placed and unplaced contigs
# - Placed contigs (scaffolded to E. coqui) serve as coverage reference
# - Unplaced contigs (not scaffolded) include candidate sex-chromosome-linked
#   and other unscaffolded sequences
#
# The UNSCAFFOLDED primary assembly is used (not the RagTag-scaffolded
# genome) so that coverage is computed per original contig without
# scaffolding artifacts.
#
# Input:
#   GrenadaFrog144_ONT_HAC_all.fastq.gz   (FASTQ — quality scores preserved)
#   final_medaka_polished_assembly_consensus.fasta.gz  (gzipped reference)
#   placed_contig_names.txt / unplaced_contig_names.txt (from RagTag)
#
# Dependencies:
#   minimap2, samtools, mosdepth, seqkit
#   Install: mamba install -c bioconda -c conda-forge minimap2 samtools mosdepth seqkit
#
# Output in coverage/:
#   primary_assembly.bam / .bam.bai       : sorted BAM + index
#   primary_cov.mosdepth.summary.txt      : per-contig coverage
#   placed_coverage.tsv                   : coverage for placed contigs
#   placed_coverage_classified.tsv        : placed contigs with coverage class
#   unplaced_coverage.tsv                 : coverage for unplaced contigs with reads
#   unplaced_coverage_classified.tsv      : unplaced contigs with coverage class
#   unmappable_zero_cov_contigs.tsv      : contigs with zero read coverage
#   missing_zero_cov.txt                  : names of zero-coverage contigs
#   samtools_stats.txt                    : comprehensive mapping statistics
#   contig_gc_content.tsv                 : per-contig GC%
#   contig_read_lengths.tsv               : per-contig read length stats
#   coverage_three_set_summary.txt        : summary statistics
#
# Three output sets:
#   1. Placed contigs (8,815) — coverage stats + classification
#   2. Unplaced contigs with mappable reads — coverage stats + classification
#   3. Zero-coverage contigs — zero coverage, reported separately
#
# Coverage classification:
#   Autosomal baseline (A) = bp-weighted mean of placed contigs
#   Hemi_0.5x:            0.375A–0.625A  (hemizygous)
#   Auto_1.0x:            0.75A–1.25A    (autosomal-like)
#   High_Coverage/Repeat:  >1.25A        (elevated coverage)
#   Low_Coverage/Other:    everything else
#   Zero:                  0x coverage    (zero-coverage contigs)
#
# Usage:
#   bash map_reads_primary_assembly_coverage.sh
#

# Author: Kopp K. Pristimantis euphronides genome project.
set -euo pipefail

# === Paths ===
BASEDIR="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"
OUTDIR="${BASEDIR}/coverage"
REF_GZ="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta.gz"
REF="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"
READS="/data/GrenadaFrog144/GrenadaFrog144_ONT_HAC_all.fastq.gz"
THREADS=24

PLACED_NAMES="/data/GrenadaFrog144/placed_contig_names.txt"
UNPLACED_NAMES="/data/GrenadaFrog144/unplaced_contig_names.txt"

# === Sanity checks ===
echo "=== Map reads to primary assembly for placed/unplaced coverage ==="
echo "Start: $(date)"
echo "Threads: $THREADS"
echo ""

# ============================================================
# Step 0: Decompress reference if needed
# ============================================================
if [ -f "$REF" ]; then
    echo "[Step 0] Reference FASTA already decompressed: $REF"
elif [ -f "$REF_GZ" ]; then
    echo "[Step 0] Decompressing reference ..."
    gunzip -k "$REF_GZ"
    echo "  Done: $REF"
else
    echo "ERROR: Neither $REF nor $REF_GZ found"
    exit 1
fi

# Index reference FASTA (needed by bcftools mpileup in BAM QC script)
if [ -f "${REF}.fai" ]; then
    echo "  FASTA index exists: ${REF}.fai"
else
    echo "  Indexing reference with samtools faidx ..."
    require_cmd samtools
    samtools faidx "$REF"
fi
echo ""

# Check all required inputs
for f in "$REF" "$READS" "$PLACED_NAMES" "$UNPLACED_NAMES"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

# Tool availability is checked per-step, only when the step needs to run.
# This prevents unnecessary failures when rerunning with completed outputs.
require_cmd() {
    for cmd in "$@"; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "ERROR: $cmd not found on PATH"
            echo "       Install via: mamba install -c bioconda -c conda-forge $cmd"
            exit 1
        fi
    done
}

mkdir -p "$OUTDIR"

BAM="${OUTDIR}/primary_assembly.bam"
PREFIX="${OUTDIR}/primary_cov"
SUMMARY="${PREFIX}.mosdepth.summary.txt"

# ============================================================
# Step 1: Map reads to primary assembly
# ============================================================
if [ -s "$BAM" ]; then
    echo "[Step 1] BAM exists, skipping: $BAM"
else
    echo "[Step 1] Mapping reads to primary assembly ..."
    require_cmd minimap2 samtools
    echo "  Ref:   $REF"
    echo "  Reads: $READS"
    echo "  Note:  FASTQ input — base quality scores will be stored in BAM"
    T0=$(date +%s)
    minimap2 -x map-ont -a -t "$THREADS" "$REF" "$READS" | \
        samtools sort -@ 12 -m 4G -o "$BAM"
    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
fi
echo ""

# ============================================================
# Step 2: Index BAM
# ============================================================
if [ -s "${BAM}.bai" ]; then
    echo "[Step 2] BAM index exists, skipping"
else
    echo "[Step 2] Indexing BAM ..."
    require_cmd samtools
    samtools index -@ "$THREADS" "$BAM"
fi
echo ""

# ============================================================
# Step 3: mosdepth per-contig coverage
# ============================================================
if [ -s "$SUMMARY" ]; then
    echo "[Step 3] mosdepth output exists, skipping"
else
    echo "[Step 3] Running mosdepth ..."
    require_cmd mosdepth
    T0=$(date +%s)
    mosdepth --fast-mode --no-per-base -t 4 \
        "$PREFIX" "$BAM"
    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
fi
echo ""

# ============================================================
# Step 3b: samtools stats — comprehensive mapping QC
# ============================================================
STATS_FILE="${OUTDIR}/samtools_stats.txt"
if [ -s "$STATS_FILE" ]; then
    echo "[Step 3b] samtools stats exists, skipping"
else
    echo "[Step 3b] Running samtools stats ..."
    require_cmd samtools
    T0=$(date +%s)
    samtools stats -@ "$THREADS" "$BAM" > "$STATS_FILE"
    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
fi
echo ""

# ============================================================
# Step 3c: Per-contig GC content (compositional bias screen)
# ============================================================
GC_OUT="${OUTDIR}/contig_gc_content.tsv"
if [ -s "$GC_OUT" ]; then
    echo "[Step 3c] GC content exists, skipping"
else
    echo "[Step 3c] Computing per-contig GC content ..."
    require_cmd seqkit
    seqkit fx2tab -g -l -n "$REF" | \
        awk -F'\t' 'BEGIN {print "contig\tlength\tGC_pct"} {
            printf "%s\t%s\t%.2f\n", $1, $2, $3
        }' > "$GC_OUT"
    echo "  Output: $GC_OUT"
fi
echo ""

# ============================================================
# Step 3d: Per-contig read length statistics
# ============================================================
RLEN_OUT="${OUTDIR}/contig_read_lengths.tsv"
if [ -s "$RLEN_OUT" ]; then
    echo "[Step 3d] Read length stats exist, skipping"
else
    echo "[Step 3d] Computing per-contig read length statistics ..."
    require_cmd samtools
    T0=$(date +%s)
    samtools view -@ "$THREADS" "$BAM" | \
        awk -F'\t' '{
            contig = $3
            rlen = length($10)
            n[contig]++
            sum[contig] += rlen
            sumsq[contig] += rlen * rlen
            if (!(contig in minl) || rlen < minl[contig]) minl[contig] = rlen
            if (rlen > maxl[contig]) maxl[contig] = rlen
        } END {
            print "contig\tn_reads\tmean_readlen\tsd_readlen\tmin_readlen\tmax_readlen"
            for (c in n) {
                mean = sum[c] / n[c]
                var = (sumsq[c] / n[c]) - (mean * mean)
                if (var < 0) var = 0
                printf "%s\t%d\t%.0f\t%.0f\t%d\t%d\n", c, n[c], mean, sqrt(var), minl[c], maxl[c]
            }
        }' | (read -r header; echo "$header"; sort -t$'\t' -k1,1) > "$RLEN_OUT"
    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
    echo "  Output: $RLEN_OUT"
fi
echo ""

# ============================================================
# Step 4: Split into three sets using exact name matching
# ============================================================
echo "[Step 4] Splitting coverage into three sets ..."
require_cmd seqkit

# Placed contigs: exact match on first column
awk -F'\t' 'NR==FNR {names[$1]; next} FNR==1 || ($1 in names)' \
    "$PLACED_NAMES" "$SUMMARY" | \
    grep -v "_region" > "${OUTDIR}/placed_coverage.tsv"

# Unplaced contigs present in mosdepth output: exact match
awk -F'\t' 'NR==FNR {names[$1]; next} FNR==1 || ($1 in names)' \
    "$UNPLACED_NAMES" "$SUMMARY" | \
    grep -v "_region" > "${OUTDIR}/unplaced_coverage.tsv"

# Identify unplaced contigs MISSING from mosdepth (zero coverage)
awk -F'\t' 'NR==FNR {names[$1]; next} !($1 in names)' \
    <(cut -f1 "${OUTDIR}/unplaced_coverage.tsv") "$UNPLACED_NAMES" \
    > "${OUTDIR}/missing_zero_cov.txt"

N_MISSING=$(wc -l < "${OUTDIR}/missing_zero_cov.txt")

# Create separate file for zero-coverage contigs with their lengths
echo -e "chrom\tlength\tbases\tmean\tmin\tmax" \
    > "${OUTDIR}/unmappable_zero_cov_contigs.tsv"

if [ "$N_MISSING" -gt 0 ]; then
    seqkit grep -f "${OUTDIR}/missing_zero_cov.txt" "$REF" | \
        seqkit fx2tab -nl | \
        awk -F'\t' '{print $1 "\t" $2 "\t0\t0.00\t0\t0"}' \
        >> "${OUTDIR}/unmappable_zero_cov_contigs.tsv"
fi

# Counts
N_PLACED=$(($(wc -l < "${OUTDIR}/placed_coverage.tsv") - 1))
N_UNPLACED=$(($(wc -l < "${OUTDIR}/unplaced_coverage.tsv") - 1))
N_TOTAL=$((N_PLACED + N_UNPLACED + N_MISSING))
N_EXPECTED_PLACED=$(wc -l < "$PLACED_NAMES")
N_EXPECTED_UNPLACED=$(wc -l < "$UNPLACED_NAMES")

echo "  Splitting done."
echo ""

# ============================================================
# Step 5: Summary statistics — all output to stdout AND file
# ============================================================
SUMFILE="${OUTDIR}/coverage_three_set_summary.txt"

{
    echo "Coverage Three-Set Summary"
    echo "=========================="
    echo "Mapping: all ONT reads → unscaffolded primary assembly"
    echo "  Reads: $READS  (FASTQ — quality scores available)"
    echo "  Ref:   $REF"
    echo "  Tool:  minimap2 -x map-ont, mosdepth --fast-mode --no-per-base"
    echo "Date: $(date)"
    echo ""

    # --- samtools stats highlights ---
    echo "=== Global mapping statistics (samtools stats) ==="
    awk -F'\t' '
    /^SN\traw total sequences:/ { printf "  Total reads:            %s\n", $3 }
    /^SN\treads mapped:/ { printf "  Reads mapped:           %s\n", $3 }
    /^SN\treads mapped and paired:/ { next }
    /^SN\treads unmapped:/ { printf "  Reads unmapped:         %s\n", $3 }
    /^SN\ttotal length:/ { printf "  Total read bases:       %s\n", $3 }
    /^SN\tbases mapped:/ && !/cigar/ { printf "  Bases mapped:           %s\n", $3 }
    /^SN\taverage length:/ { printf "  Average read length:    %s\n", $3 }
    /^SN\tmaximum length:/ { printf "  Maximum read length:    %s\n", $3 }
    /^SN\taverage quality:/ { printf "  Average base quality:   %s\n", $3 }
    /^SN\terror rate:/ { printf "  Error rate:             %s\n", $3 }
    /^SN\tbases mapped \(cigar\):/ { printf "  Bases mapped (CIGAR):   %s\n", $3 }
    /^SN\tmismatches:/ { printf "  Mismatches:             %s\n", $3 }
    /^SN\tsupplementary alignments:/ { printf "  Supplementary alns:     %s\n", $3 }
    ' "$STATS_FILE"
    echo ""

    echo "=== Contig set counts ==="
    echo "  Placed contigs:                    $N_PLACED"
    echo "  Unplaced contigs (with coverage):  $N_UNPLACED"
    echo "  Unplaced contigs (zero coverage):  $N_MISSING"
    echo "  Total accounted for:               $N_TOTAL"
    echo ""
    echo "  Expected placed:   $N_EXPECTED_PLACED"
    echo "  Expected unplaced: $N_EXPECTED_UNPLACED (= $N_UNPLACED + $N_MISSING)"
    echo ""
    echo "=== Coverage statistics ==="
    echo ""
    echo "  Per-contig:  unweighted stats of mosdepth per-contig mean (col 4)"
    echo "  Genome-wide: base-pair-weighted mean (total mapped bases / total length)"
    echo ""

    # --- Reusable awk for uniform per-contig + bp-weighted stats ---
    STATS_AWK='
    NR>1 {
        n++; a[n] = $4; bp += $2; mapped_bases += $3
        sum += $4; sumsq += $4*$4
        if (n == 1 || $4 < vmin) vmin = $4
        if (n == 1 || $4 > vmax) vmax = $4
    }
    END {
        if (n == 0) { printf "  n=0\n"; exit }
        mean = sum/n; var = sumsq/n - mean*mean
        if (var < 0) var = 0
        asort(a)
        if (n%2) med = a[(n+1)/2]; else med = (a[n/2] + a[n/2+1]) / 2
        bpw = (bp > 0) ? mapped_bases/bp : 0
        printf "  n=%d, total=%d bp\n", n, bp
        printf "  genome-wide mean=%.2fx (bp-weighted)\n", bpw
        printf "  per-contig:  mean=%.2fx, median=%.2fx, sd=%.2f\n", mean, med, sqrt(var)
        printf "  per-contig:  min=%.2fx, max=%.2fx\n", vmin, vmax
    }'

    # Overall genome (all contigs: placed + unplaced + zero-cov)
    echo "Overall genome (all contigs):"
    cat <(tail -n+2 "${OUTDIR}/placed_coverage.tsv") \
        <(tail -n+2 "${OUTDIR}/unplaced_coverage.tsv") \
        <(tail -n+2 "${OUTDIR}/unmappable_zero_cov_contigs.tsv") | \
        awk -F'\t' 'BEGIN {print "h\tlen\tb\tmean\tmin\tmax"} {print}' | \
        awk -F'\t' "$STATS_AWK"
    echo ""

    echo "Set 1 — Placed contigs:"
    awk -F'\t' "$STATS_AWK" "${OUTDIR}/placed_coverage.tsv"
    echo ""

    echo "Set 2 — Unplaced contigs (with coverage):"
    awk -F'\t' "$STATS_AWK" "${OUTDIR}/unplaced_coverage.tsv"
    echo ""

    echo "Set 3 — Zero-coverage contigs:"
    awk -F'\t' "$STATS_AWK" "${OUTDIR}/unmappable_zero_cov_contigs.tsv"
    if [ "$N_MISSING" -gt 0 ]; then
        echo "  (no mapped reads)"
    fi
    echo ""

    # --- GC content summary by set ---
    echo "=== GC content by contig set ==="
    for LABEL in placed unplaced; do
        NAMEFILE="/data/GrenadaFrog144/${LABEL}_contig_names.txt"
        echo "${LABEL}:"
        awk -F'\t' 'NR==FNR {names[$1]; next} NR>1 && ($1 in names) {
            n++; sum += $3
            if (n == 1 || $3 < vmin) vmin = $3
            if (n == 1 || $3 > vmax) vmax = $3
        } END {
            if (n > 0) printf "  n=%d, mean GC=%.2f%%, min=%.2f%%, max=%.2f%%\n", n, sum/n, vmin, vmax
        }' "$NAMEFILE" "$GC_OUT"
    done
    echo ""

    echo "=== Coverage distribution ==="
    for LABEL in placed unplaced; do
        echo ""
        echo "${LABEL}:"
        awk -F'\t' 'NR>1 {
            c = $4
            if (c < 5) bin="a_<5x"
            else if (c < 15) bin="b_5-15x"
            else if (c < 25) bin="c_15-25x"
            else if (c < 35) bin="d_25-35x"
            else if (c < 50) bin="e_35-50x"
            else if (c < 100) bin="f_50-100x"
            else if (c < 500) bin="g_100-500x"
            else bin="h_>500x"
            count[bin]++
        } END {
            for (b in count) printf "  %10s: %6d\n", substr(b,3), count[b]
        }' "${OUTDIR}/${LABEL}_coverage.tsv" | sort
    done
    echo ""
    echo "=== Hypothesis-driven coverage classification ==="
    echo ""
    echo "Autosomal baseline (A): Set 1 bp-weighted mean coverage"
    echo "Classes defined relative to A (multiplier bands):"
    echo "  Hemi_0.5x:           0.375A <= cov <= 0.625A  (hemizygous)"
    echo "  Auto_1.0x:           0.75A  <= cov <= 1.25A   (autosomal-like)"
    echo "  High_Coverage/Repeat: cov > 1.25A             (elevated coverage)"
    echo "  Low_Coverage/Other:   everything else          (below Hemi or in gaps)"
    echo "  Zero:                 coverage = 0x            (zero-coverage contigs)"
    echo ""

    # Compute autosomal baseline from placed contigs (bp-weighted mean)
    AUTO_BASE=$(awk -F'\t' 'NR>1 { bp += $2; bases += $3 } END { printf "%.2f", bases/bp }' \
        "${OUTDIR}/placed_coverage.tsv")
    echo "  Autosomal baseline (A) = ${AUTO_BASE}x"
    echo "  Note: A is the bp-weighted mean of placed contigs on the unscaffolded primary assembly."
    echo "  Hemi band:  $(awk "BEGIN {printf \"%.2f\", 0.375 * $AUTO_BASE}")x – $(awk "BEGIN {printf \"%.2f\", 0.625 * $AUTO_BASE}")x"
    echo "  Auto band:  $(awk "BEGIN {printf \"%.2f\", 0.75 * $AUTO_BASE}")x – $(awk "BEGIN {printf \"%.2f\", 1.25 * $AUTO_BASE}")x"
    echo "  High:       > $(awk "BEGIN {printf \"%.2f\", 1.25 * $AUTO_BASE}")x"
    echo ""

    # --- Classify UNPLACED contigs ---
    CLASSFILE_UNPLACED="${OUTDIR}/unplaced_coverage_classified.tsv"

    {
        echo -e "contig\tlength\tbases\tmean_cov\tmin\tmax\tcoverage_class"
        # Contigs with coverage
        tail -n+2 "${OUTDIR}/unplaced_coverage.tsv"
        # Zero-coverage contigs
        tail -n+2 "${OUTDIR}/unmappable_zero_cov_contigs.tsv"
    } | awk -F'\t' -v A="$AUTO_BASE" '
    NR==1 { print; next }
    {
        c = $4
        if (c == 0)                                     cls = "Zero"
        else if (c >= 0.375 * A && c <= 0.625 * A)     cls = "Hemi_0.5x"
        else if (c >= 0.75 * A  && c <= 1.25 * A)      cls = "Auto_1.0x"
        else if (c > 1.25 * A)                          cls = "High_Coverage/Repeat"
        else                                            cls = "Low_Coverage/Other"
        print $0 "\t" cls
    }' > "$CLASSFILE_UNPLACED"

    echo "Unplaced contig classification:"
    awk -F'\t' 'NR>1 {
        cls = $7; n[cls]++; bp[cls] += $2
    } END {
        printf "  %-22s %8s %14s\n", "Class", "Contigs", "Total bp"
        printf "  %-22s %8s %14s\n", "-----", "-------", "--------"
        split("Hemi_0.5x,Auto_1.0x,High_Coverage/Repeat,Low_Coverage/Other,Zero", order, ",")
        total_n = 0; total_bp = 0
        for (i = 1; i <= 5; i++) {
            c = order[i]
            if (n[c] > 0) {
                printf "  %-22s %8d %14d\n", c, n[c], bp[c]
                total_n += n[c]; total_bp += bp[c]
            }
        }
        printf "  %-22s %8s %14s\n", "-----", "-------", "--------"
        printf "  %-22s %8d %14d\n", "Total", total_n, total_bp
    }' "$CLASSFILE_UNPLACED"
    echo "  Output: $CLASSFILE_UNPLACED"
    echo ""

    # --- Classify PLACED contigs (same baseline, same bands) ---
    CLASSFILE_PLACED="${OUTDIR}/placed_coverage_classified.tsv"

    awk -F'\t' -v A="$AUTO_BASE" '
    NR==1 { print $0 "\tcoverage_class"; next }
    {
        c = $4
        if (c == 0)                                     cls = "Zero"
        else if (c >= 0.375 * A && c <= 0.625 * A)     cls = "Hemi_0.5x"
        else if (c >= 0.75 * A  && c <= 1.25 * A)      cls = "Auto_1.0x"
        else if (c > 1.25 * A)                          cls = "High_Coverage/Repeat"
        else                                            cls = "Low_Coverage/Other"
        print $0 "\t" cls
    }' "${OUTDIR}/placed_coverage.tsv" > "$CLASSFILE_PLACED"

    echo "Placed contig classification:"
    awk -F'\t' 'NR>1 {
        cls = $NF; n[cls]++; bp[cls] += $2
    } END {
        printf "  %-22s %8s %14s\n", "Class", "Contigs", "Total bp"
        printf "  %-22s %8s %14s\n", "-----", "-------", "--------"
        split("Hemi_0.5x,Auto_1.0x,High_Coverage/Repeat,Low_Coverage/Other,Zero", order, ",")
        total_n = 0; total_bp = 0
        for (i = 1; i <= 5; i++) {
            c = order[i]
            if (n[c] > 0) {
                printf "  %-22s %8d %14d\n", c, n[c], bp[c]
                total_n += n[c]; total_bp += bp[c]
            }
        }
        printf "  %-22s %8s %14s\n", "-----", "-------", "--------"
        printf "  %-22s %8d %14d\n", "Total", total_n, total_bp
    }' "$CLASSFILE_PLACED"
    echo "  Output: $CLASSFILE_PLACED"
    echo ""

    echo "=== Output files ==="
    echo "  ${OUTDIR}/primary_assembly.bam"
    echo "  ${OUTDIR}/primary_assembly.bam.bai"
    echo "  ${OUTDIR}/placed_coverage.tsv"
    echo "  ${OUTDIR}/placed_coverage_classified.tsv"
    echo "  ${OUTDIR}/unplaced_coverage.tsv"
    echo "  ${OUTDIR}/unplaced_coverage_classified.tsv"
    echo "  ${OUTDIR}/unmappable_zero_cov_contigs.tsv"
    echo "  ${OUTDIR}/missing_zero_cov.txt"
    echo "  ${OUTDIR}/coverage_three_set_summary.txt"
    echo "  ${OUTDIR}/samtools_stats.txt"
    echo "  ${OUTDIR}/contig_gc_content.tsv"
    echo "  ${OUTDIR}/contig_read_lengths.tsv"
} 2>&1 | tee "$SUMFILE"

echo ""
echo "Finished: $(date)"
echo "=== Done ==="
