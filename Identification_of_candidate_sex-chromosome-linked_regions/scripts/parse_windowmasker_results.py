#!/usr/bin/env python3
"""
parse_windowmasker_results.py
=============================
Parse WindowMasker interval output for two contig sets and produce:
  - Per-contig masked fraction tables (TSV)
  - Combined summary statistics

Usage:
    python3 parse_windowmasker_results.py \
        --outdir /path/to/windowmasker_output \
        --placed /path/to/placed_contigs.fasta \
        --unplaced /path/to/unplaced_contigs.fasta

Date: 2026-02-16
Author: Kopp K, Pristimantis euphronides genome project
"""

import argparse
import os
import sys
import numpy as np
from collections import defaultdict


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def get_contig_lengths(fasta_path):
    """Parse FASTA to get contig lengths."""
    lengths = {}
    name = None
    length = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
    return lengths


# ---------------------------------------------------------------------------
# WindowMasker interval parsing
# ---------------------------------------------------------------------------

def parse_windowmasker_intervals(filepath):
    """
    Parse WindowMasker -outfmt interval output.
    Format: >seqname followed by lines of 'start - end'
    """
    masked = defaultdict(list)
    current_seq = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_seq = line[1:].split()[0]
            elif ' - ' in line and current_seq:
                parts = line.split(' - ')
                try:
                    start = int(parts[0].strip())
                    end = int(parts[1].strip())
                    masked[current_seq].append((start, end))
                except ValueError:
                    continue
    return dict(masked)


def compute_masked_bp(intervals_list):
    """Sum masked bp from a list of (start, end) intervals."""
    return sum(end - start + 1 for start, end in intervals_list)


# ---------------------------------------------------------------------------
# Process one set
# ---------------------------------------------------------------------------

def process_set(label, fasta_path, outdir):
    """
    Process one contig set: parse its FASTA and WindowMasker intervals,
    write per-contig TSV, return summary dict.
    """
    intervals_file = os.path.join(outdir, f"wm_{label}_intervals.txt")

    if not os.path.exists(fasta_path):
        print(f"  {label}: FASTA not found ({fasta_path}), skipping")
        return None
    if not os.path.exists(intervals_file):
        print(f"  {label}: intervals file not found ({intervals_file}), skipping")
        return None

    print(f"  Processing {label} ...")

    contig_lengths = get_contig_lengths(fasta_path)
    masked_intervals = parse_windowmasker_intervals(intervals_file)

    # Write per-contig table
    output_tsv = os.path.join(outdir, f"wm_{label}_per_contig.tsv")
    total_len = 0
    total_masked = 0
    masked_fracs = []

    with open(output_tsv, 'w') as out:
        out.write("contig\tlength\tmasked_bp\tmasked_fraction\n")
        for contig, length in sorted(contig_lengths.items(),
                                      key=lambda x: x[1], reverse=True):
            m_bp = compute_masked_bp(masked_intervals.get(contig, []))
            frac = m_bp / length if length > 0 else 0.0
            out.write(f"{contig}\t{length}\t{m_bp}\t{frac:.4f}\n")
            total_len += length
            total_masked += m_bp
            masked_fracs.append(frac)

    fracs = np.array(masked_fracs)
    n = len(fracs)

    print(f"    {n:,} contigs, {total_len / 1e6:.0f} Mb total, "
          f"{total_masked / 1e6:.0f} Mb masked ({total_masked / total_len:.1%})")
    print(f"    Per-contig: median={np.median(fracs):.1%}, "
          f"mean={np.mean(fracs):.1%}, sd={np.std(fracs):.1%}")
    print(f"    >80% masked: {np.sum(fracs > 0.80):,}, "
          f">90%: {np.sum(fracs > 0.90):,}, "
          f"<20%: {np.sum(fracs < 0.20):,}")
    print(f"    Output: {output_tsv}")

    return {
        'n': n,
        'total_bp': total_len,
        'total_masked': total_masked,
        'overall_frac': total_masked / total_len if total_len > 0 else 0,
        'median': float(np.median(fracs)),
        'mean': float(np.mean(fracs)),
        'std': float(np.std(fracs)),
        'gt80': int(np.sum(fracs > 0.80)),
        'gt90': int(np.sum(fracs > 0.90)),
        'lt20': int(np.sum(fracs < 0.20)),
    }


# ---------------------------------------------------------------------------
# Write combined summary
# ---------------------------------------------------------------------------

def write_summary(summaries, outdir):
    """Write side-by-side summary table."""
    summary_file = os.path.join(outdir, "wm_all_sets_summary.txt")
    labels = [l for l in ['placed', 'unplaced']
              if l in summaries]

    with open(summary_file, 'w') as out:
        out.write("WindowMasker Two-Set Comparison\n")
        out.write("Counts calibrated from: full primary assembly\n")
        out.write("  (final_medaka_polished_assembly_consensus.fasta)\n")
        out.write("=" * 75 + "\n\n")

        header = f"{'Metric':<30s}" + "".join(f"  {l:>16s}" for l in labels)
        out.write(header + "\n")
        out.write("-" * len(header) + "\n")

        rows = [
            ("Contigs",              'n',            lambda x: f"{x:,}"),
            ("Total (Mb)",           'total_bp',     lambda x: f"{x / 1e6:.0f}"),
            ("Masked (Mb)",          'total_masked', lambda x: f"{x / 1e6:.0f}"),
            ("Overall masked %",     'overall_frac', lambda x: f"{x:.1%}"),
            ("Per-contig median",    'median',       lambda x: f"{x:.1%}"),
            ("Per-contig mean",      'mean',         lambda x: f"{x:.1%}"),
            ("Per-contig sd",        'std',          lambda x: f"{x:.1%}"),
            ("Contigs >80% masked",  'gt80',         lambda x: f"{x:,}"),
            ("Contigs >90% masked",  'gt90',         lambda x: f"{x:,}"),
            ("Contigs <20% masked",  'lt20',         lambda x: f"{x:,}"),
        ]

        for metric_name, key, fmt in rows:
            line = f"{metric_name:<30s}"
            for l in labels:
                line += f"  {fmt(summaries[l][key]):>16s}"
            out.write(line + "\n")

        out.write("\n" + "=" * 75 + "\n")

    print(f"\n  Summary table: {summary_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Parse WindowMasker intervals for contig sets')
    parser.add_argument('--outdir', required=True,
                        help='Directory containing wm_*_intervals.txt files')
    parser.add_argument('--placed', required=True,
                        help='Placed contigs FASTA')
    parser.add_argument('--unplaced', required=True,
                        help='Unplaced contigs FASTA')
    args = parser.parse_args()

    sets = [
        ('placed',         args.placed),
        ('unplaced',       args.unplaced)
    ]

    summaries = {}
    for label, fasta_path in sets:
        result = process_set(label, fasta_path, args.outdir)
        if result is not None:
            summaries[label] = result

    if summaries:
        write_summary(summaries, args.outdir)

    print("  Done.")


if __name__ == '__main__':
    main()
