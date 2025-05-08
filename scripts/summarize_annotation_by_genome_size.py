#!/usr/bin/env python3
"""
summarize_annotation_by_genome_size.py

Summarizes genome annotation features from a HANNO `.bedDB` file and computes:
  - Gene counts per functional category
  - Total coding length (bp) per category
  - Percentage of the scaffolded genome occupied by those coding bases

Inputs:
  1. HANNO .bedDB file (tab-delimited)
  2. Total scaffolded genome size (in bp)

Usage:
  python3 summarize_annotation_by_genome_size.py BESTMODELS-FINAL.bedDB 1748533034

Output:
  Printed summary table with category-wise gene counts, coding bp, and % genome

Author: K. Kopp / Assistant
Date: 2025
"""

import pandas as pd
import sys

# Check arguments
if len(sys.argv) != 2 and len(sys.argv) != 3:
    print("Usage: python3 summarize_annotation_by_genome_size.py <bedDB_file> <genome_size_bp>")
    sys.exit(1)

beddb_path = sys.argv[1]

try:
    genome_size = int(sys.argv[2])
except:
    print("Error: genome_size must be an integer (e.g., 1748533034).")
    sys.exit(1)

# Define expected columns
colnames = [
    "chrom", "start", "end", "transcript_id", "score", "strand",
    "cds_start", "cds_end", "color", "num_exons", "exon_sizes", "exon_starts",
    "gene_id", "gene_name", "transcript_name", "source", "feature", "frame",
    "tag", "product", "description", "go_terms", "pfam", "kog", "note", "eggNOG",
    "preferred_name", "gene_symbol", "eggnog_function", "eggnog_name", "eggnog_desc",
    "COG_cat", "eggNOG_OGs", "EC_number", "KEGG_KOs", "KEGG_Pathway", "BiGG_Reaction",
    "PFAMs", "busco1", "busco2", "busco3", "busco4", "orflen", "mRNAlen",
    "cds_exons_count", "mRNA_exon_count"
]

# Load and preprocess
beddb = pd.read_csv(
    beddb_path,
    sep='\t',
    names=colnames,
    comment='#',
    dtype=str,
    low_memory=False,
    na_values=['-']
)

# Extract gene ID without .t1 suffix
beddb["base_gene_id"] = beddb["transcript_id"].str.extract(r'(hanno\.g\d+)')

# Drop rows with missing CDS coordinates
beddb = beddb[pd.notna(beddb['cds_start']) & pd.notna(beddb['cds_end'])].copy()

# Convert CDS start/end to integer for calculation
beddb["cds_start"] = beddb["cds_start"].astype(int)
beddb["cds_end"] = beddb["cds_end"].astype(int)

# Group by gene
grouped = beddb.groupby("base_gene_id", group_keys=False)
total_genes = grouped.ngroups

# Total CDS length across all protein-coding genes
total_coding_bp = grouped[["cds_start", "cds_end"]].apply(
    lambda df: df.drop_duplicates().apply(lambda x: x["cds_end"] - x["cds_start"], axis=1).sum()
).sum()

# Define functional category checks
def has_function(desc): return any(pd.notna(d) and "hypothetical protein" not in d.lower() for d in desc)
def has_pfam(val): return any(pd.notna(x) for x in val)
def has_kog(val): return any(pd.notna(x) for x in val)
def has_kegg(val): return any(pd.notna(x) for x in val)
def has_go(val): return any(pd.notna(x) for x in val)
def has_product(desc): return any(pd.notna(x) for x in desc)

# Category definitions
categories = {
    "Function prediction":     ('description', has_function),
    "Pfam domains":            ('pfam',        has_pfam),
    "KOG assignments":         ('kog',         has_kog),
    "KEGG pathway mapping":    ('KEGG_Pathway',has_kegg),
    "GO term annotation":      ('go_terms',    has_go),
    "Product descriptions":    ('description', has_product)
}

# Process each category
summary_rows = []
for label, (field, check_fn) in categories.items():
    gene_count = grouped[field].apply(check_fn).sum()
    cds_bp = grouped[["cds_start", "cds_end", field]].apply(
        lambda df: (
            check_fn(df[field]) and df[["cds_start", "cds_end"]].drop_duplicates().apply(
                lambda x: x["cds_end"] - x["cds_start"], axis=1
            ).sum()
        )
    ).sum()
    percent = round((cds_bp / genome_size) * 100, 2)
    summary_rows.append([label, gene_count, cds_bp, percent])

# Add overall row
summary_rows.append([
    "All protein-coding genes",
    total_genes,
    total_coding_bp,
    round((total_coding_bp / genome_size) * 100, 2)
])

# Format output
summary = pd.DataFrame(summary_rows, columns=["Category", "Gene Count", "Coding bp", "% of Genome"])

# Print result
print(f"\nGenome size used: {genome_size:,} bp")
print("\nAnnotation Summary by Genome Occupancy (CDS-based):\n")
print(summary.to_string(index=False))
