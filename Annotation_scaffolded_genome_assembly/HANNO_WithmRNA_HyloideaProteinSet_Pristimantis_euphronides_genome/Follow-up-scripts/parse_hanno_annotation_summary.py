#!/usr/bin/env python3
"""
parse_hanno_annotation_summary.py

Compute structural and functional annotation statistics from a HANNO
BESTMODELS bedDB file, after removing coordinate-identical duplicate
gene models (same chrom, start, end, strand).

Input:
    BESTMODELS-FINAL.renamed.bedDB  (HANNO v0.4 output, tab-separated)

Output:
    annotation_summary.tsv
    annotation_summary_duplicates.txt  (list of removed duplicate IDs)

The bedDB is BED12 with appended annotation columns. Relevant fields
(0-based):

     0  Chrom           9  blockCount       24  Description (eggNOG)
     1  mRNAstart      10  blockSizes       25  Preferred_name
     2  mRNAend        11  blockStarts      26  GOs
     3  Name            6  CDSstart         28  KEGG_ko
     5  strand          7  CDSend           29  KEGG_Pathway
    12  origName       21  eggNOG_OGs       37  PFAMs
    23  COG_category

Transcript names: hanno.gN.tM (N = gene locus, M = isoform).
In BESTMODELS output there is one transcript per gene locus.

Duplicate detection: entries with identical (chrom, start, end, strand)
are flagged; the first occurrence is kept, subsequent ones are skipped.
This matches the logic in beddb_to_gff3.py to ensure both scripts
report identical gene counts.

Evidence source from origName (col 12):
    TU*    -> protein-based (miniprot)
    MSTRG* -> transcript-based (minimap2/StringTie)

Usage:
    python3 parse_hanno_annotation_summary.py \\
        BESTMODELS-FINAL.renamed.bedDB \\
        --genome-size 1748533034 \\
        -o annotation_summary.tsv
        
Author: K. Kopp, P. euphronides genome project
        
"""

import argparse
import sys


def parse_sizes(s):
    """Parse comma-separated integer list, ignoring trailing comma."""
    return [int(x) for x in s.rstrip(',').split(',') if x]


def main():
    ap = argparse.ArgumentParser(
        description='Compute annotation summary from HANNO BESTMODELS bedDB '
                    '(with duplicate removal)')
    ap.add_argument('beddb')
    ap.add_argument('--genome-size', type=int, default=1748533034)
    ap.add_argument('-o', '--output', default='annotation_summary.tsv')
    args = ap.parse_args()

    gsize = args.genome_size

    # --- Pass 1: read all entries, detect duplicates ---
    raw_entries = []
    with open(args.beddb) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 38:
                continue
            raw_entries.append(f)

    n_raw = len(raw_entries)

    # Sort by chrom, start, end (matching beddb_to_gff3.py)
    raw_entries.sort(key=lambda f: (f[0], int(f[1]), int(f[2])))

    # Mark duplicates: same chrom, start, end, strand -> keep first, skip rest
    keep = [True] * n_raw
    n_dups = 0
    dup_names = []
    for i in range(1, n_raw):
        prev = raw_entries[i - 1]
        curr = raw_entries[i]
        if (prev[0] == curr[0] and        # chrom
            prev[1] == curr[1] and         # start
            prev[2] == curr[2] and         # end
            prev[5] == curr[5]):           # strand
            keep[i] = False
            n_dups += 1
            dup_names.append(curr[3])

    print(f'Raw entries: {n_raw}', file=sys.stderr)
    print(f'Duplicates removed: {n_dups}', file=sys.stderr)
    print(f'Entries retained: {n_raw - n_dups}', file=sys.stderr)

    # --- Pass 2: accumulate stats on non-duplicate entries ---
    n_genes = 0
    total_gene_span = 0
    total_exons = 0
    total_exon_bp = 0
    total_cds_features = 0
    total_cds_bp = 0

    n_protein = 0
    n_transcript = 0

    n_pref_name = 0
    n_func_desc = 0
    n_eggnog_ogs = 0
    n_pfam = 0
    n_go = 0
    n_cog = 0
    n_kegg_ko = 0
    n_kegg_pathway = 0

    def valid(v):
        return v not in ('-', '', '--')

    for i, f in enumerate(raw_entries):
        if not keep[i]:
            continue

        n_genes += 1

        mrna_start = int(f[1])
        mrna_end = int(f[2])
        cds_start = int(f[6])
        cds_end = int(f[7])
        block_count = int(f[9])
        block_sizes = parse_sizes(f[10])
        block_starts = parse_sizes(f[11])
        orig_name = f[12]

        # Structural: gene/mRNA span
        total_gene_span += (mrna_end - mrna_start)

        # Exons
        total_exons += block_count
        total_exon_bp += sum(block_sizes)

        # CDS: exon blocks intersected with [cds_start, cds_end]
        for j in range(block_count):
            exon_s = mrna_start + block_starts[j]
            exon_e = exon_s + block_sizes[j]
            cs = max(exon_s, cds_start)
            ce = min(exon_e, cds_end)
            if cs < ce:
                total_cds_features += 1
                total_cds_bp += (ce - cs)

        # Evidence source
        if orig_name.startswith('MSTRG'):
            n_transcript += 1
        else:
            n_protein += 1

        # Functional annotation
        eggnog_ogs = f[21] if len(f) > 21 else '-'
        cog_cat    = f[23] if len(f) > 23 else '-'
        desc       = f[24] if len(f) > 24 else '-'
        pref_name  = f[25] if len(f) > 25 else '-'
        gos        = f[26] if len(f) > 26 else '-'
        kegg_ko    = f[28] if len(f) > 28 else '-'
        kegg_pw    = f[29] if len(f) > 29 else '-'
        pfams      = f[37] if len(f) > 37 else '-'

        if valid(pref_name):
            n_pref_name += 1

        has_func = False
        if valid(desc) and 'hypothetical' not in desc.lower():
            has_func = True
        if valid(pfams) or valid(gos) or valid(kegg_ko) or valid(cog_cat):
            has_func = True
        if has_func:
            n_func_desc += 1

        if valid(eggnog_ogs):
            n_eggnog_ogs += 1
        if valid(pfams):
            n_pfam += 1
        if valid(gos):
            n_go += 1
        if valid(cog_cat):
            n_cog += 1
        if valid(kegg_ko):
            n_kegg_ko += 1
        if valid(kegg_pw):
            n_kegg_pathway += 1

    # Derived values
    mean_exons = total_exons / n_genes if n_genes else 0
    mean_cds = total_cds_bp / n_genes if n_genes else 0
    mean_mrna = total_exon_bp / n_genes if n_genes else 0

    pg = lambda bp: 100.0 * bp / gsize
    pn = lambda n: 100.0 * n / n_genes

    # Write TSV
    with open(args.output, 'w') as out:
        w = out.write
        w('Category\tCount\tTotal_length_bp\tPct_of_genome\tPct_of_genes\n')
        w(f'Protein-coding gene models\t{n_genes}\t{total_gene_span}\t{pg(total_gene_span):.2f}\t\n')
        w(f'mRNA\t{n_genes}\t{total_gene_span}\t{pg(total_gene_span):.2f}\t\n')
        w(f'Exons\t{total_exons}\t{total_exon_bp}\t{pg(total_exon_bp):.2f}\t\n')
        w(f'CDS features\t{total_cds_features}\t{total_cds_bp}\t{pg(total_cds_bp):.2f}\t\n')
        w(f'Mean exons per gene\t{mean_exons:.2f}\t\t\t\n')
        w(f'Mean CDS length per gene\t{mean_cds:.0f} bp\t\t\t\n')
        w(f'Mean mRNA length per gene\t{mean_mrna:.0f} bp\t\t\t\n')
        w(f'Protein-based models (miniprot)\t{n_protein}\t\t\t{pn(n_protein):.1f}\n')
        w(f'Transcript-based models (minimap2)\t{n_transcript}\t\t\t{pn(n_transcript):.1f}\n')
        w(f'Genes with gene symbol (Preferred_name)\t{n_pref_name}\t\t\t{pn(n_pref_name):.1f}\n')
        w(f'Genes with functional description\t{n_func_desc}\t\t\t{pn(n_func_desc):.1f}\n')
        w(f'Genes with eggNOG orthologous groups\t{n_eggnog_ogs}\t\t\t{pn(n_eggnog_ogs):.1f}\n')
        w(f'Genes with Pfam domains\t{n_pfam}\t\t\t{pn(n_pfam):.1f}\n')
        w(f'Genes with GO terms\t{n_go}\t\t\t{pn(n_go):.1f}\n')
        w(f'Genes with COG category\t{n_cog}\t\t\t{pn(n_cog):.1f}\n')
        w(f'Genes with KEGG orthologs (KO)\t{n_kegg_ko}\t\t\t{pn(n_kegg_ko):.1f}\n')
        w(f'Genes with KEGG pathway mapping\t{n_kegg_pathway}\t\t\t{pn(n_kegg_pathway):.1f}\n')

    # Write duplicate list if any
    if dup_names:
        dup_file = args.output.replace('.tsv', '_duplicates.txt')
        with open(dup_file, 'w') as df:
            df.write('# Coordinate-identical duplicate gene models removed\n')
            df.write('# (same chrom, start, end, strand as a preceding entry)\n')
            for name in dup_names:
                df.write(name + '\n')
        print(f'Duplicate list: {dup_file}', file=sys.stderr)

    # Stdout summary
    print(f'Genes: {n_genes}  (after removing {n_dups} duplicates from {n_raw})')
    print(f'Gene/mRNA span: {total_gene_span:,} bp ({pg(total_gene_span):.2f}%)')
    print(f'Exons: {total_exons:,}  total bp: {total_exon_bp:,} ({pg(total_exon_bp):.2f}%)')
    print(f'CDS features: {total_cds_features:,}  total bp: {total_cds_bp:,} ({pg(total_cds_bp):.2f}%)')
    print(f'Mean exons/gene: {mean_exons:.2f}')
    print(f'Mean CDS length: {mean_cds:.0f} bp')
    print(f'Mean mRNA length: {mean_mrna:.0f} bp')
    print(f'Protein: {n_protein} ({pn(n_protein):.1f}%)  Transcript: {n_transcript} ({pn(n_transcript):.1f}%)')
    print(f'Gene symbol: {n_pref_name}  Func desc: {n_func_desc}  eggNOG: {n_eggnog_ogs}')
    print(f'Pfam: {n_pfam}  GO: {n_go}  COG: {n_cog}  KEGG KO: {n_kegg_ko}  KEGG pw: {n_kegg_pathway}')
    print(f'Output: {args.output}')


if __name__ == '__main__':
    main()
