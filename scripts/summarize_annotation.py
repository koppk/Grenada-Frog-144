#!/usr/bin/env python3
"""
summarize_annotation.py

Summarize genome annotation content from a HANNO `.bedDB` file.

This script computes the number and percentage of protein-coding genes in the annotation that:
- have predicted functional annotations (non-hypothetical)
- have Pfam domains
- are assigned to KOG categories
- are mapped to KEGG pathways
- are annotated with Gene Ontology (GO) terms
- have any product description

Input: A tab-separated `.bedDB` file produced by the HANNO pipeline.
Output: Printed table summarizing annotation completeness by category.

Usage:
    python3 summarize_annotation.py BESTMODELS-FINAL.bedDB

Author: Katharina Kopp
Date: 2025-05-08
"""

import pandas as pd
import sys

# Check for correct usage
if len(sys.argv) != 2:
    print("Usage: python3 summarize_annotation.py <bedDB_file>")
    sys.exit(1)

beddb_path = sys.argv[1]

# Define column names (matching HANNO .bedDB format)
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

# Load the annotation file
beddb = pd.read_csv(beddb_path, sep='\t', names=colnames, comment='#', low_memory=False)

# Extract gene-level identifier by removing transcript suffix (e.g., .t1, .t2)
beddb["base_gene_id"] = beddb["transcript_id"].str.extract(r'(hanno\.g\d+)')

# Group all transcript entries by gene
grouped = beddb.groupby("base_gene_id")
total_genes = grouped.ngroups

# Function: gene has a meaningful product description (not hypothetical)
def has_function(descs):
    return any(pd.notna(d) and "hypothetical protein" not in d.lower() and d.strip() != "-" for d in descs)

# Compute annotation coverage across categories
genes_with_function = grouped['description'].apply(has_function).sum()
genes_with_pfam     = grouped['pfam'].apply(lambda x: any(pd.notna(i) and i.strip() != "-" for i in x)).sum()
genes_with_kog      = grouped['kog'].apply(lambda x: any(pd.notna(i) and i.strip() != "-" for i in x)).sum()
genes_with_kegg     = grouped['KEGG_Pathway'].apply(lambda x: any(pd.notna(i) and i.strip() != "-" for i in x)).sum()
genes_with_go       = grouped['go_terms'].apply(lambda x: any(pd.notna(i) and i.strip() != "-" for i in x)).sum()
genes_with_product  = grouped['description'].apply(lambda x: any(pd.notna(i) and i.strip() != "-" for i in x)).sum()

# Format output
summary = pd.DataFrame({
    "Category": [
        "Total protein-coding genes",
        "Genes with function prediction",
        "Genes with Pfam domains",
        "Genes with KOG assignments",
        "Genes with KEGG pathway mapping",
        "Genes with GO term annotation",
        "Genes with product descriptions"
    ],
    "Count": [
        total_genes,
        genes_with_function,
        genes_with_pfam,
        genes_with_kog,
        genes_with_kegg,
        genes_with_go,
        genes_with_product
    ],
    "Percentage": [
        100.0,
        round((genes_with_function / total_genes) * 100, 1),
        round((genes_with_pfam / total_genes) * 100, 1),
        round((genes_with_kog / total_genes) * 100, 1),
        round((genes_with_kegg / total_genes) * 100, 1),
        round(
