#!/bin/bash
#
# compute_arithmetic_mean_phred.sh
#
# Computes the per-read arithmetic mean of Phred quality scores from
# gzipped FASTQ files. For each read, all per-base Phred scores are
# summed and divided by the number of bases; these per-read means are
# then averaged across all reads (per-read mean). In addition, all
# per-base Phred scores across all reads are summed and divided by the
# total number of bases (read-length-weighted mean). These are the same
# calculations applied per contig in the assembly quality assessment.
#
# Input:  One or more gzipped FASTQ files
# Output: arithmetic_mean_phred_summary.tsv
#
# Dependencies: pigz, python3
#
# Usage:
#   bash compute_arithmetic_mean_phred.sh \
#       GrenadaFrog144_ONT_ALL.fastq.gz \
#       GrenadaFrog144_ONT_HAC_all.fastq.gz
#
# Author: Kopp K, Pristimantis euphronides genome project
set -euo pipefail

THREADS=16
OUTFILE="arithmetic_mean_phred_summary.tsv"

if [ $# -eq 0 ]; then
    echo "Usage: $0 <fastq.gz> [fastq.gz ...]"
    exit 1
fi

for cmd in pigz python3; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd not found" >&2
        exit 1
    fi
done

echo -e "file\tn_reads\tper_read_mean_Q\tread_length_weighted_mean_Q\ttotal_bases" > "$OUTFILE"

for FASTQ in "$@"; do
    if [ ! -f "$FASTQ" ]; then
        echo "ERROR: $FASTQ not found" >&2
        exit 1
    fi

    BASENAME=$(basename "$FASTQ")
    echo "Processing: $BASENAME"
    T0=$(date +%s)

    pigz -dc -p "$THREADS" "$FASTQ" \
        | awk 'NR%4==0' \
        | python3 -c "
import sys

n = 0
total_q = 0.0       # sum of per-read means (for per-read mean)
total_bases = 0      # sum of read lengths (for read-length-weighted mean)
total_qbases = 0     # sum of all per-base Phred scores (for read-length-weighted mean)

for line in sys.stdin:
    qual = line.strip().encode()
    qlen = len(qual)

    # FASTQ encodes Phred quality scores as ASCII characters with an
    # offset of 33 (Sanger/Phred+33 encoding). Each ASCII character
    # value minus 33 gives the Phred score for that base.
    # Example: ASCII 'D' = 68, Phred score = 68 - 33 = Q35
    qsum = sum(qual) - 33 * qlen

    # Per-read mean: sum of per-base Phred scores / number of bases
    mean_q = qsum / qlen

    n += 1
    total_q += mean_q        # accumulate per-read means
    total_bases += qlen       # accumulate base count
    total_qbases += qsum      # accumulate Phred score sum

# Per-read mean: average of per-read means (each read counts equally)
# Read-length-weighted mean: total Phred sum / total bases (each base counts equally)
print(f'${BASENAME}\t{n}\t{total_q/n:.2f}\t{total_qbases/total_bases:.2f}\t{total_bases}')
" >> "$OUTFILE"

    T1=$(date +%s)
    echo "  $BASENAME: done in $((T1 - T0))s"
done

echo ""
column -t -s $'\t' "$OUTFILE"
echo ""
echo "Output: $OUTFILE"
