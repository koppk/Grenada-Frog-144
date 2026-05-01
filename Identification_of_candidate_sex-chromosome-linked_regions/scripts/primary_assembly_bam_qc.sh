#!/bin/bash
#
# primary_assembly_bam_qc.sh
# ==========================
# Extract assembly quality metrics from the BAM of all reads mapped
# to the unscaffolded primary assembly.
#
# Requires: primary_assembly.bam created by map_reads_primary_assembly_coverage.sh
#           GNU parallel (for parallelized variant calling)
#
# Analyses:
#   1. Assembly contiguity: N50, L50, total size, contig count, largest
#   2. Coverage uniformity: per-contig coefficient of variation (1 kb windows)
#   3. Mapping quality: per-contig mean MAPQ distribution
#   4. Base quality: per-contig mean read base quality (from FASTQ via BAM)
#   5. Heterozygosity: per-contig variant density (bcftools, quality-aware)
#   6. Soft-clipping & supplementary: per-contig rates
#   7. Alignment identity: mean read-to-assembly identity (NanoStat on BAM)
#
# Output: /data/GrenadaFrog144/assembly_qc/
#

# Author: Kopp K. Pristimantis euphronides genome project.
set -euo pipefail

# === Paths ===
BAM="/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison/coverage/primary_assembly.bam"
REF="/data/GrenadaFrog144/final_medaka_polished_assembly_consensus.fasta"
OUTDIR="/data/GrenadaFrog144/assembly_qc"
THREADS=24

# === Sanity checks ===
echo "=== Primary assembly BAM QC ==="
echo "Start: $(date)"
echo "Threads: $THREADS"
echo ""

if [ ! -f "$BAM" ]; then
    echo "ERROR: $BAM not found"
    echo "       Run map_reads_primary_assembly_coverage.sh first."
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "ERROR: $REF not found"
    echo "       Decompress .fasta.gz first (map_reads script does this)."
    exit 1
fi

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

# ============================================================
# 1. Assembly contiguity statistics
# ============================================================
echo "[1/7] Assembly contiguity statistics ..."

CONTIG_OUT="${OUTDIR}/assembly_contiguity.txt"

if [ -s "$CONTIG_OUT" ]; then
    echo "  Output exists, skipping: $CONTIG_OUT"
else
    echo "  Computing N50, L50, total size ..."
    require_cmd seqkit

    # Get contig lengths from reference FASTA
    seqkit fx2tab -nl "$REF" | awk -F'\t' '{print $2}' | sort -rn \
        > "${OUTDIR}/contig_lengths_sorted.txt"

    awk '
    {
        len[NR] = $1
        total += $1
        n++
    }
    END {
        # N50/L50
        cumul = 0
        for (i = 1; i <= n; i++) {
            cumul += len[i]
            if (cumul >= total / 2) {
                n50 = len[i]
                l50 = i
                break
            }
        }
        # N90/L90
        cumul = 0
        for (i = 1; i <= n; i++) {
            cumul += len[i]
            if (cumul >= total * 0.9) {
                n90 = len[i]
                l90 = i
                break
            }
        }
        printf "Assembly Contiguity Statistics\n"
        printf "==============================\n"
        printf "  Total contigs:      %d\n", n
        printf "  Total size:         %d bp (%.2f Gb)\n", total, total/1e9
        printf "  Largest contig:     %d bp (%.2f Mb)\n", len[1], len[1]/1e6
        printf "  Smallest contig:    %d bp\n", len[n]
        printf "  N50:                %d bp (%.2f Mb)\n", n50, n50/1e6
        printf "  L50:                %d contigs\n", l50
        printf "  N90:                %d bp (%.2f kb)\n", n90, n90/1e3
        printf "  L90:                %d contigs\n", l90
        printf "  Mean contig size:   %d bp (%.2f kb)\n", total/n, total/n/1e3
        printf "\n"

        # Size distribution bins
        for (i = 1; i <= n; i++) {
            if      (len[i] >= 10000000) b["a_>=10Mb"]++
            else if (len[i] >= 1000000)  b["b_1-10Mb"]++
            else if (len[i] >= 100000)   b["c_100kb-1Mb"]++
            else if (len[i] >= 10000)    b["d_10-100kb"]++
            else if (len[i] >= 1000)     b["e_1-10kb"]++
            else                         b["f_<1kb"]++
        }
        printf "  Size distribution:\n"
        split("a_>=10Mb,b_1-10Mb,c_100kb-1Mb,d_10-100kb,e_1-10kb,f_<1kb", order, ",")
        for (i = 1; i <= 6; i++) {
            c = order[i]
            if (b[c] > 0) printf "    %-12s  %6d contigs\n", substr(c,3), b[c]
        }
    }' "${OUTDIR}/contig_lengths_sorted.txt" > "$CONTIG_OUT"

    cat "$CONTIG_OUT"
    echo "  Output: $CONTIG_OUT"
fi
echo ""

# ============================================================
# 2. Coverage uniformity: per-contig CV (1 kb windows)
# ============================================================
echo "[2/7] Per-contig coverage uniformity ..."

QUAL_OUT="${OUTDIR}/contig_coverage_uniformity.tsv"

if [ -s "$QUAL_OUT" ]; then
    echo "  Output exists, skipping: $QUAL_OUT"
else
    echo "  Computing coverage in 1kb windows, then CV per contig ..."
    require_cmd mosdepth

    # Use mosdepth with 1kb windows for within-contig variation
    mosdepth --fast-mode --no-per-base --by 1000 -t 4 \
        "${OUTDIR}/qc_1kb_windows" "$BAM"

    # Compute mean, sd, CV per contig from the windowed output
    zcat "${OUTDIR}/qc_1kb_windows.regions.bed.gz" | \
        awk -F'\t' '{
            contig = $1
            cov = $4
            n[contig]++
            sum[contig] += cov
            sumsq[contig] += cov * cov
        } END {
            print "contig\tn_windows\tmean_cov\tsd_cov\tCV"
            for (c in n) {
                mean = sum[c] / n[c]
                variance = (sumsq[c] / n[c]) - (mean * mean)
                if (variance < 0) variance = 0
                sd = sqrt(variance)
                cv = (mean > 0) ? sd / mean : 0
                printf "%s\t%d\t%.2f\t%.2f\t%.4f\n", c, n[c], mean, sd, cv
            }
        }' | (read -r header; echo "$header"; sort -t$'\t' -k5 -rn) > "$QUAL_OUT"

    echo "  Output: $QUAL_OUT"
fi
echo ""

# ============================================================
# 3. Mapping quality: per-contig mean MAPQ
# ============================================================
echo "[3/7] Per-contig mapping quality ..."

MAPQ_OUT="${OUTDIR}/contig_mapq.tsv"

if [ -s "$MAPQ_OUT" ]; then
    echo "  Output exists, skipping: $MAPQ_OUT"
else
    echo "  Extracting MAPQ per alignment ..."
    require_cmd samtools

    samtools view -@ "$THREADS" "$BAM" | \
        awk -F'\t' '{
            contig = $3
            mapq = $5
            n[contig]++
            sum[contig] += mapq
            if (mapq == 0) q0[contig]++
            if (mapq >= 60) q60[contig]++
        } END {
            print "contig\tn_alignments\tmean_mapq\tfrac_mapq0\tfrac_mapq60"
            for (c in n) {
                mean = sum[c] / n[c]
                f0 = (q0[c] + 0) / n[c]
                f60 = (q60[c] + 0) / n[c]
                printf "%s\t%d\t%.2f\t%.4f\t%.4f\n", c, n[c], mean, f0, f60
            }
        }' | (read -r header; echo "$header"; sort -t$'\t' -k3 -n) > "$MAPQ_OUT"

    echo "  Output: $MAPQ_OUT"
fi
echo ""

# ============================================================
# 4. Base quality: per-contig mean read base quality
# ============================================================
echo "[4/7] Per-contig base quality (from FASTQ quality scores in BAM) ..."

BQ_OUT="${OUTDIR}/contig_base_quality.tsv"

if [ -s "$BQ_OUT" ]; then
    echo "  Output exists, skipping: $BQ_OUT"
else
    echo "  Extracting per-read mean base quality ..."
    require_cmd samtools python3
    T0=$(date +%s)

    # Python is ~50x faster than awk here: sum(bytes) is a single C call
    # vs awk iterating character-by-character over every quality string.
    # For 17M ONT reads averaging ~10kb, that's the difference between
    # hours and minutes.
    samtools view -@ "$THREADS" "$BAM" | \
        python3 -c '
import sys
from collections import defaultdict

n = defaultdict(int)
sum_q = defaultdict(float)
sumsq_q = defaultdict(float)
sum_bases = defaultdict(int)
sum_qbases = defaultdict(int)

for line in sys.stdin:
    f = line.split("\t", 12)
    contig = f[2]
    qual = f[10]
    if qual == "*":
        continue
    qbytes = qual.encode()
    qlen = len(qbytes)
    qsum = sum(qbytes) - 33 * qlen
    mean_q = qsum / qlen

    n[contig] += 1
    sum_q[contig] += mean_q
    sumsq_q[contig] += mean_q * mean_q
    sum_bases[contig] += qlen
    sum_qbases[contig] += qsum

print("contig\tn_reads\tmean_baseQ\tsd_baseQ\tbp_weighted_meanQ")
for c in sorted(n.keys(), key=lambda x: sum_q[x]/n[x]):
    mean = sum_q[c] / n[c]
    var = (sumsq_q[c] / n[c]) - (mean * mean)
    if var < 0: var = 0
    bpw = sum_qbases[c] / sum_bases[c] if sum_bases[c] > 0 else 0
    print(f"{c}\t{n[c]}\t{mean:.2f}\t{var**0.5:.2f}\t{bpw:.2f}")
' > "$BQ_OUT"

    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
    echo "  Output: $BQ_OUT"
fi
echo ""

# ============================================================
# 5. Heterozygosity: per-contig variant density
# ============================================================
echo "[5/7] Per-contig heterozygosity (variant calling) ..."

VCF_OUT="${OUTDIR}/variants.vcf.gz"
HET_OUT="${OUTDIR}/contig_heterozygosity.tsv"

if [ -s "$HET_OUT" ]; then
    echo "  Output exists, skipping: $HET_OUT"
else
    if [ ! -s "$VCF_OUT" ]; then
        echo "  Running bcftools mpileup + call in parallel (${THREADS} chunks) ..."
        echo "  Quality-aware: -Q 20 filters bases below Q20 from pileup"
        require_cmd samtools bcftools parallel

        CHUNKS_DIR="${OUTDIR}/mpileup_chunks"
        mkdir -p "$CHUNKS_DIR"

        # Create region BED files: distribute contigs round-robin across chunks
        samtools view -H "$BAM" | grep "^@SQ" | \
            awk -F'\t' '{
                for (i=1; i<=NF; i++) {
                    if ($i ~ /^SN:/) sn = substr($i, 4)
                    if ($i ~ /^LN:/) ln = substr($i, 4)
                }
                print sn "\t0\t" ln
            }' | awk -v n="$THREADS" -v dir="$CHUNKS_DIR" '{
                chunk = (NR - 1) % n
                print >> dir "/chunk_" chunk ".bed"
            }'

        T0=$(date +%s)

        # Run mpileup+call per chunk in parallel
        # -d 200 caps pileup depth to avoid stalling on collapsed repeats
        seq 0 $((THREADS - 1)) | parallel -j "$THREADS" \
            "bcftools mpileup -f $REF -q 10 -Q 20 -d 200 \
                --regions-file ${CHUNKS_DIR}/chunk_{}.bed $BAM 2>/dev/null | \
             bcftools call -mv -Oz -o ${CHUNKS_DIR}/chunk_{}.vcf.gz 2>/dev/null && \
             bcftools index ${CHUNKS_DIR}/chunk_{}.vcf.gz"

        T1=$(date +%s)
        echo "  Parallel mpileup done in $((T1 - T0)) seconds"

        # Concatenate all chunks
        echo "  Concatenating VCF chunks ..."
        bcftools concat -a ${CHUNKS_DIR}/chunk_*.vcf.gz -Oz -o "$VCF_OUT"
        bcftools index "$VCF_OUT"

        # Clean up chunk files
        rm -rf "$CHUNKS_DIR"
        echo "  Chunks cleaned up"
    fi

    echo "  Computing per-contig variant density ..."
    require_cmd samtools bcftools

    # Get contig lengths
    samtools view -H "$BAM" | grep "^@SQ" | \
        awk -F'\t' '{
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SN:/) sn = substr($i, 4)
                if ($i ~ /^LN:/) ln = substr($i, 4)
            }
            print sn "\t" ln
        }' > "${OUTDIR}/contig_lengths_from_bam.txt"

    # Count variants per contig
    bcftools query -f '%CHROM\n' "$VCF_OUT" | \
        sort | uniq -c | awk '{print $2 "\t" $1}' > "${OUTDIR}/variants_per_contig.txt"

    # Count het and hom-alt variants separately
    bcftools query -f '%CHROM\t[%GT]\n' "$VCF_OUT" | \
        awk -F'\t' '{
            contig = $1
            gt = $2
            total[contig]++
            if (gt == "0/1" || gt == "0|1") het[contig]++
            else hom[contig]++
        } END {
            for (c in total) printf "%s\t%d\t%d\t%d\n", c, total[c], het[c]+0, hom[c]+0
        }' | sort -k1,1 > "${OUTDIR}/variants_het_hom_per_contig.txt"

    # Join: contig, length, n_variants, n_het, n_hom, variants_per_kb, het_per_kb
    join -t$'\t' -a1 \
        <(sort -k1,1 "${OUTDIR}/contig_lengths_from_bam.txt") \
        <(sort -k1,1 "${OUTDIR}/variants_het_hom_per_contig.txt") | \
        awk -F'\t' 'BEGIN {print "contig\tlength\tn_variants\tn_het\tn_hom\tvariants_per_kb\thet_per_kb"} {
            len = $2
            nvar = ($3 == "") ? 0 : $3
            nhet = ($4 == "") ? 0 : $4
            nhom = ($5 == "") ? 0 : $5
            per_kb = (len > 0) ? nvar / (len / 1000) : 0
            het_kb = (len > 0) ? nhet / (len / 1000) : 0
            printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\n", $1, len, nvar, nhet, nhom, per_kb, het_kb
        }' | (read -r header; echo "$header"; sort -t$'\t' -k6 -rn) > "$HET_OUT"

    echo "  Output: $HET_OUT"
fi
echo ""

# ============================================================
# 6. Soft-clipping and supplementary alignment rates
# ============================================================
echo "[6/7] Per-contig soft-clipping and supplementary rates ..."

CLIP_OUT="${OUTDIR}/contig_clipping_supplementary.tsv"

if [ -s "$CLIP_OUT" ]; then
    echo "  Output exists, skipping: $CLIP_OUT"
else
    echo "  Parsing CIGAR strings and flags ..."
    require_cmd samtools python3
    T0=$(date +%s)

    samtools view -@ "$THREADS" "$BAM" | \
        python3 -c '
import sys, re
from collections import defaultdict

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
CONSUMES_QUERY = set("MIS=X")

n = defaultdict(int)
supp = defaultdict(int)
sec = defaultdict(int)
clip_frac = defaultdict(float)
hard_sum = defaultdict(int)
total_sum = defaultdict(int)

for line in sys.stdin:
    f = line.split("\t", 7)
    flag = int(f[1])
    contig = f[2]
    cigar = f[5]

    n[contig] += 1

    if flag & 2048:
        supp[contig] += 1
    if flag & 256:
        sec[contig] += 1

    total_bases = 0
    soft_bases = 0
    hard_bases = 0
    for length_str, op in CIGAR_RE.findall(cigar):
        length = int(length_str)
        if op in CONSUMES_QUERY:
            total_bases += length
        if op == "S":
            soft_bases += length
        elif op == "H":
            hard_bases += length

    if total_bases > 0:
        clip_frac[contig] += soft_bases / total_bases
    hard_sum[contig] += hard_bases
    total_sum[contig] += total_bases + hard_bases

print("contig\tn_alignments\tfrac_supplementary\tfrac_secondary\tmean_softclip_frac\tmean_hardclip_frac")
for c in sorted(n.keys(), key=lambda x: clip_frac[x]/n[x], reverse=True):
    f_supp = supp[c] / n[c]
    f_sec = sec[c] / n[c]
    mean_clip = clip_frac[c] / n[c]
    hard_frac = hard_sum[c] / total_sum[c] if total_sum[c] > 0 else 0
    print(f"{c}\t{n[c]}\t{f_supp:.4f}\t{f_sec:.4f}\t{mean_clip:.4f}\t{hard_frac:.4f}")
' > "$CLIP_OUT"

    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
    echo "  Output: $CLIP_OUT"
fi
echo ""

# ============================================================
# 7. Alignment identity: NanoStat on BAM
# ============================================================
echo "[7/7] Read-to-assembly alignment identity (NanoStat) ..."

NANOSTAT_OUT="${OUTDIR}/primary_assembly_nanostat.txt"

if [ -s "$NANOSTAT_OUT" ]; then
    echo "  Output exists, skipping: $NANOSTAT_OUT"
else
    require_cmd NanoStat
    T0=$(date +%s)

    NanoStat --bam "$BAM" --threads "$THREADS" -n "$NANOSTAT_OUT"

    T1=$(date +%s)
    echo "  Done in $((T1 - T0)) seconds"
    echo "  Output: $NANOSTAT_OUT"
fi
echo ""

# ============================================================
# Summary report — output to stdout AND file
# ============================================================
SUMFILE="${OUTDIR}/assembly_qc_summary.txt"

{
    echo "============================================================"
    echo "  Primary Assembly Quality Assessment"
    echo "  Pristimantis euphronides (GrenadaFrog144)"
    echo "  $(date)"
    echo "============================================================"
    echo ""
    echo "Reference: $REF"
    echo "Reads mapped: all ONT HAC reads → unscaffolded primary assembly"
    echo ""

    # === 1. Assembly contiguity ===
    echo "============================================================"
    echo "1. ASSEMBLY CONTIGUITY"
    echo "============================================================"
    if [ -s "$CONTIG_OUT" ]; then
        cat "$CONTIG_OUT"
    fi
    echo ""

    # === 2. Coverage uniformity ===
    echo "============================================================"
    echo "2. COVERAGE UNIFORMITY (1 kb windows, coefficient of variation)"
    echo "============================================================"
    echo ""
    echo "  CV < 0.5:  uniform coverage (well-assembled unique sequence)"
    echo "  CV 0.5-1:  moderate variation (structural variants, CNVs, or"
    echo "             heterozygous indels within contig)"
    echo "  CV > 1:    highly variable (collapsed repeats, chimeric joins,"
    echo "             or misassembly)"
    echo ""
    if [ -s "$QUAL_OUT" ]; then
        awk -F'\t' 'NR>1 {
            n++; sum+=$5
            if ($5 < 0.5) low++
            else if ($5 < 1.0) mid++
            else high++
        } END {
            printf "  Overall:  mean CV = %.4f across %d contigs\n", sum/n, n
            printf "  CV < 0.5 (uniform):    %6d contigs (%5.1f%%)\n", low+0, (low+0)/n*100
            printf "  CV 0.5-1 (moderate):   %6d contigs (%5.1f%%)\n", mid+0, (mid+0)/n*100
            printf "  CV > 1   (variable):   %6d contigs (%5.1f%%)\n", high+0, (high+0)/n*100
        }' "$QUAL_OUT"

        echo ""
        echo "  Top 10 most variable contigs (potential misassemblies):"
        awk -F'\t' 'NR>1 && NR<=11 {
            printf "    %-20s  %d windows  mean_cov=%.1fx  CV=%.4f\n", $1, $2, $3, $5
        }' "$QUAL_OUT"
    fi
    echo ""

    # === 3. Mapping quality ===
    echo "============================================================"
    echo "3. MAPPING QUALITY (MAPQ)"
    echo "============================================================"
    echo ""
    echo "  MAPQ 60:   unique, confident placement"
    echo "  MAPQ 0:    multi-mapped (repetitive sequence)"
    echo "  Low mean:  contig contains mostly repetitive sequence"
    echo ""
    if [ -s "$MAPQ_OUT" ]; then
        awk -F'\t' 'NR>1 {
            n++; sum_mapq+=$3; sum_f0+=$4; sum_f60+=$5
            if ($3 >= 50) hq++
            else if ($3 >= 20) mq++
            else lq++
        } END {
            printf "  Overall:  mean MAPQ = %.2f across %d contigs\n", sum_mapq/n, n
            printf "  Mean frac MAPQ=0:   %.4f (multi-mapped reads)\n", sum_f0/n
            printf "  Mean frac MAPQ>=60: %.4f (uniquely mapped reads)\n", sum_f60/n
            printf "\n"
            printf "  Mean MAPQ >= 50 (high-confidence): %6d contigs (%5.1f%%)\n", hq+0, (hq+0)/n*100
            printf "  Mean MAPQ 20-50 (moderate):        %6d contigs (%5.1f%%)\n", mq+0, (mq+0)/n*100
            printf "  Mean MAPQ < 20  (repetitive):      %6d contigs (%5.1f%%)\n", lq+0, (lq+0)/n*100
        }' "$MAPQ_OUT"

        echo ""
        echo "  Bottom 10 contigs by MAPQ (most repetitive):"
        awk -F'\t' 'NR>1 && NR<=11 {
            printf "    %-20s  %d alns  mean_MAPQ=%.2f  frac_Q0=%.4f\n", $1, $2, $3, $4
        }' "$MAPQ_OUT"
    fi
    echo ""

    # === 4. Base quality ===
    echo "============================================================"
    echo "4. READ BASE QUALITY (Phred scores from FASTQ)"
    echo "============================================================"
    echo ""
    echo "  ONT HAC v4.3.0 typically produces Q15-Q25 modal quality."
    echo "  Mean Q > 20: good; Q > 30: excellent (likely from R10.4.1)"
    echo "  Contigs attracting only low-Q reads may indicate"
    echo "  problematic regions (homopolymers, low-complexity)."
    echo ""
    if [ -s "$BQ_OUT" ]; then
        awk -F'\t' 'NR>1 {
            n++; sum+=$3; sum_bpw+=$5
            if ($3 >= 30) hq++
            else if ($3 >= 20) mq++
            else lq++
        } END {
            printf "  Overall:  mean per-read Q = %.2f across %d contigs\n", sum/n, n
            printf "  Overall:  mean bp-weighted Q = %.2f\n", sum_bpw/n
            printf "\n"
            printf "  Mean Q >= 30 (excellent): %6d contigs (%5.1f%%)\n", hq+0, (hq+0)/n*100
            printf "  Mean Q 20-30 (good):      %6d contigs (%5.1f%%)\n", mq+0, (mq+0)/n*100
            printf "  Mean Q < 20  (low):       %6d contigs (%5.1f%%)\n", lq+0, (lq+0)/n*100
        }' "$BQ_OUT"
    fi
    echo ""

    # === 5. Heterozygosity ===
    echo "============================================================"
    echo "5. HETEROZYGOSITY (variant density from bcftools mpileup/call)"
    echo "============================================================"
    echo ""
    echo "  Variant calling: bcftools mpileup -q10 -Q20 -d200 | call -mv"
    echo "  Base quality filtering is meaningful (FASTQ input)."
    echo ""
    echo "  Expected for diploid amphibian: 2-10 het/kb typical"
    echo "  Elevated het/kb: possible collapsed haplotype or"
    echo "    paralogous mapping (especially if combined with high"
    echo "    coverage and low MAPQ)"
    echo "  Elevated hom-alt: assembly errors or fixed differences"
    echo "    from reference haplotype"
    echo ""
    if [ -s "$HET_OUT" ]; then
        awk -F'\t' 'NR>1 {
            n++; sum_var+=$6; sum_het+=$7
            nvar += $3; nhet += $4; nhom += $5; totlen += $2
            if ($6 == 0) zero++
            else if ($6 < 5) low++
            else if ($6 < 15) mid++
            else high++
        } END {
            printf "  Total variants called:  %d (%d het, %d hom-alt)\n", nvar, nhet, nhom
            printf "  Genome-wide density:    %.2f variants/kb (%.2f het/kb)\n", nvar/(totlen/1000), nhet/(totlen/1000)
            printf "  Per-contig mean:        %.2f variants/kb (%.2f het/kb)\n", sum_var/n, sum_het/n
            if (nvar > 0) printf "  Het/total ratio:        %.4f\n", nhet/nvar
            printf "\n"
            printf "  Zero variants:          %6d contigs (%5.1f%%)\n", zero+0, (zero+0)/n*100
            printf "  < 5 var/kb (low):       %6d contigs (%5.1f%%)\n", low+0, (low+0)/n*100
            printf "  5-15 var/kb (typical):  %6d contigs (%5.1f%%)\n", mid+0, (mid+0)/n*100
            printf "  > 15 var/kb (elevated): %6d contigs (%5.1f%%)\n", high+0, (high+0)/n*100
        }' "$HET_OUT"

        echo ""
        echo "  Top 10 most heterozygous contigs:"
        awk -F'\t' 'NR>1 && NR<=11 {
            printf "    %-20s  %s bp  %d vars  het/kb=%.2f  hom/kb=%.2f\n", $1, $2, $3, $7, ($2>0) ? $5/($2/1000) : 0
        }' "$HET_OUT"
    fi
    echo ""

    # === 6. Soft-clipping and supplementary alignments ===
    echo "============================================================"
    echo "6. SOFT-CLIPPING AND SUPPLEMENTARY ALIGNMENTS"
    echo "============================================================"
    echo ""
    echo "  Soft-clipping: read extends beyond contig or diverges"
    echo "    from reference (structural variant, misassembly boundary)."
    echo "  Supplementary: read maps to multiple locations (split"
    echo "    alignment across chimeric join or structural variant)."
    echo "  Secondary: read maps equally well elsewhere (repeat)."
    echo ""
    echo "  High soft-clip + high supplementary: strong signal for"
    echo "    misassembly or structural rearrangement."
    echo ""
    if [ -s "$CLIP_OUT" ]; then
        awk -F'\t' 'NR>1 {
            n++; sum_clip+=$5; sum_supp+=$3; sum_sec+=$4; sum_hard+=$6
            if ($5 < 0.1) low_clip++
            else if ($5 < 0.3) mid_clip++
            else high_clip++
        } END {
            printf "  Overall:  mean soft-clip fraction = %.4f\n", sum_clip/n
            printf "  Overall:  mean supplementary rate = %.4f\n", sum_supp/n
            printf "  Overall:  mean secondary rate     = %.4f\n", sum_sec/n
            printf "  Overall:  mean hard-clip fraction  = %.4f\n", sum_hard/n
            printf "\n"
            printf "  Soft-clip < 0.1 (clean):    %6d contigs (%5.1f%%)\n", low_clip+0, (low_clip+0)/n*100
            printf "  Soft-clip 0.1-0.3 (some):   %6d contigs (%5.1f%%)\n", mid_clip+0, (mid_clip+0)/n*100
            printf "  Soft-clip > 0.3 (high):     %6d contigs (%5.1f%%)\n", high_clip+0, (high_clip+0)/n*100
        }' "$CLIP_OUT"

        echo ""
        echo "  Top 10 contigs by soft-clipping (potential misassemblies):"
        awk -F'\t' 'NR>1 && NR<=11 {
            printf "    %-20s  %d alns  clip=%.4f  supp=%.4f  sec=%.4f\n", $1, $2, $5, $3, $4
        }' "$CLIP_OUT"
    fi
    echo ""

    # === Alignment identity ===
    echo "============================================================"
    echo "7. READ-TO-ASSEMBLY ALIGNMENT IDENTITY (NanoStat)"
    echo "============================================================"
    echo ""
    if [ -s "$NANOSTAT_OUT" ]; then
        grep -E "^(Mean read quality|Average percent identity|Fraction of bases aligned|Number of reads|Total bases aligned):" "$NANOSTAT_OUT" | \
            sed 's/^/  /'
    else
        echo "  NanoStat output not found: $NANOSTAT_OUT"
    fi
    echo ""

    # === Cross-metric flagging ===
    echo "============================================================"
    echo "8. CROSS-METRIC QUALITY FLAGS"
    echo "============================================================"
    echo ""
    echo "  Contigs flagged as potentially problematic if they meet"
    echo "  2+ of: CV > 1, mean MAPQ < 10, soft-clip > 0.3, het/kb > 20"
    echo ""

    if [ -s "$QUAL_OUT" ] && [ -s "$MAPQ_OUT" ] && [ -s "$CLIP_OUT" ] && [ -s "$HET_OUT" ]; then
        # Join all metrics by contig name and flag multi-criterion outliers
        join -t$'\t' \
            <(awk -F'\t' 'NR>1 {print $1"\t"$5}' "$QUAL_OUT" | sort -k1,1) \
            <(awk -F'\t' 'NR>1 {print $1"\t"$3}' "$MAPQ_OUT" | sort -k1,1) | \
        join -t$'\t' - \
            <(awk -F'\t' 'NR>1 {print $1"\t"$5}' "$CLIP_OUT" | sort -k1,1) | \
        join -t$'\t' - \
            <(awk -F'\t' 'NR>1 {print $1"\t"$6}' "$HET_OUT" | sort -k1,1) | \
        awk -F'\t' '
        BEGIN { print "contig\tCV\tmean_MAPQ\tsoftclip_frac\tvar_per_kb\tflags" }
        {
            flags = 0; reasons = ""
            if ($2 > 1)   { flags++; reasons = reasons "highCV " }
            if ($3 < 10)  { flags++; reasons = reasons "lowMAPQ " }
            if ($4 > 0.3) { flags++; reasons = reasons "highClip " }
            if ($5 > 20)  { flags++; reasons = reasons "highHet " }
            if (flags >= 2) printf "%s\t%.4f\t%.2f\t%.4f\t%.2f\t%s\n", $1, $2, $3, $4, $5, reasons
        }' > "${OUTDIR}/flagged_contigs.tsv"

        N_FLAGGED=$(($(wc -l < "${OUTDIR}/flagged_contigs.tsv") - 1))
        echo "  Flagged contigs (2+ criteria): $N_FLAGGED"
        echo "  Output: ${OUTDIR}/flagged_contigs.tsv"
        echo ""
        if [ "$N_FLAGGED" -gt 0 ]; then
            echo "  Top 20 flagged contigs:"
            head -21 "${OUTDIR}/flagged_contigs.tsv" | \
                awk -F'\t' 'NR==1 { printf "    %-20s  %8s  %8s  %8s  %10s  %s\n", $1, $2, $3, $4, $5, $6; next }
                { printf "    %-20s  %8s  %8s  %8s  %10s  %s\n", $1, $2, $3, $4, $5, $6 }'
        fi
    fi
    echo ""

    # === Output file inventory ===
    echo "============================================================"
    echo "OUTPUT FILES"
    echo "============================================================"
    echo ""
    echo "  Per-contig metrics (TSV, header row, sortable):"
    echo "    $QUAL_OUT"
    echo "    $MAPQ_OUT"
    echo "    $BQ_OUT"
    echo "    $HET_OUT"
    echo "    $CLIP_OUT"
    echo ""
    echo "  Flagged contigs:"
    echo "    ${OUTDIR}/flagged_contigs.tsv"
    echo ""
    echo "  Assembly stats:"
    echo "    $CONTIG_OUT"
    echo ""
    echo "  Alignment identity:"
    echo "    $NANOSTAT_OUT"
    echo ""
    echo "  This summary:"
    echo "    $SUMFILE"
    echo ""
} 2>&1 | tee "$SUMFILE"
