#!/usr/bin/env python3
"""
filter_hanno_fasta.py

Filter HANNO BESTMODELS FASTA files (protein, CDS, mRNA) to retain only
sequences whose IDs appear in a retained_ids.txt list produced by
beddb_to_gff3.py. This ensures all submission files are consistent
after duplicate removal.

HANNO FASTA headers use the format:
    >hanno.gN.tM  (matching the transcript name in the bedDB)

Usage:
    python3 filter_hanno_fasta.py \\
        --ids retained_ids.txt \\
        --protein BESTMODELS-FINAL.AA.faa \\
        --cds BESTMODELS-FINAL.CDS.fa \\
        --mrna BESTMODELS-FINAL.mRNA.fa \\
        -o output_dir/

Author: K. Kopp, P. euphronides genome project
"""

import argparse
import sys
import os


def load_ids(path):
    """Load retained transcript IDs from a text file (one per line)."""
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                ids.add(line)
    return ids


def filter_fasta(input_path, output_path, retain_ids):
    """
    Filter a FASTA file, keeping only sequences whose header ID
    (first whitespace-delimited token after '>') is in retain_ids.

    Returns (total_seqs, retained_seqs).
    """
    total = 0
    retained = 0
    keep = False

    with open(input_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                total += 1
                seq_id = line[1:].split()[0]
                # HANNO appends strand to FASTA headers: hanno.g1.t1(-)
                # Strip trailing (+) or (-) to match bedDB transcript names
                if seq_id.endswith('(+)') or seq_id.endswith('(-)'):
                    seq_id = seq_id[:-3]
                keep = seq_id in retain_ids
                if keep:
                    retained += 1
                    fout.write(line)
            else:
                if keep:
                    fout.write(line)

    return total, retained


def main():
    ap = argparse.ArgumentParser(
        description='Filter HANNO FASTA files to match deduplicated gene set')
    ap.add_argument('--ids', required=True,
                    help='retained_ids.txt from beddb_to_gff3.py')
    ap.add_argument('--protein', default=None,
                    help='BESTMODELS-FINAL.AA.faa')
    ap.add_argument('--cds', default=None,
                    help='BESTMODELS-FINAL.CDS.fa')
    ap.add_argument('--mrna', default=None,
                    help='BESTMODELS-FINAL.mRNA.fa')
    ap.add_argument('-o', '--outdir', default='.',
                    help='Output directory (default: current directory)')
    ap.add_argument('--prefix', default='aPriEup1.0',
                    help='File name prefix (default: aPriEup1.0)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    retain_ids = load_ids(args.ids)
    print(f'Loaded {len(retain_ids)} retained IDs from {args.ids}',
          file=sys.stderr)

    files_to_filter = []
    if args.protein:
        files_to_filter.append(
            (args.protein, f'{args.prefix}_protein.faa', 'Protein'))
    if args.cds:
        files_to_filter.append(
            (args.cds, f'{args.prefix}_cds_from_genomic.fna', 'CDS'))
    if args.mrna:
        files_to_filter.append(
            (args.mrna, f'{args.prefix}_rna_from_genomic.fna', 'mRNA'))

    if not files_to_filter:
        print('ERROR: No input FASTA files specified. Use --protein, '
              '--cds, and/or --mrna.', file=sys.stderr)
        sys.exit(1)

    for input_path, out_name, label in files_to_filter:
        if not os.path.exists(input_path):
            print(f'WARNING: {input_path} not found, skipping {label}',
                  file=sys.stderr)
            continue

        output_path = os.path.join(args.outdir, out_name)
        total, retained = filter_fasta(input_path, output_path, retain_ids)
        removed = total - retained
        print(f'{label}: {total} -> {retained} sequences '
              f'({removed} duplicates removed) -> {output_path}',
              file=sys.stderr)

        if retained != len(retain_ids):
            print(f'  WARNING: Expected {len(retain_ids)} sequences but '
                  f'retained {retained}. Check ID format match.',
                  file=sys.stderr)

    print('Done.', file=sys.stderr)


if __name__ == '__main__':
    main()
