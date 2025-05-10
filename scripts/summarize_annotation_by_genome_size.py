#!/usr/bin/env python3

import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python3 summarize_annotation_by_genome_size.py <bedDB_file> <genome_size_bp>")
    sys.exit(1)

beddb_path = sys.argv[1]
try:
    genome_size = int(sys.argv[2])
except ValueError:
    print("ERROR: genome_size_bp must be an integer (e.g., 1748533034)")
    sys.exit(1)

# Load and preprocess
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

bed = pd.read_csv(
    beddb_path, sep="\t", names=colnames, comment="#",
    dtype=str, na_values=["-"], low_memory=False
)

bed["base_gene_id"] = bed["transcript_id"].str.extract(r"(hanno\.g\d+)")
bed["cds_start"] = pd.to_numeric(bed["cds_start"], errors="coerce")
bed["cds_end"] = pd.to_numeric(bed["cds_end"], errors="coerce")
bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
bed["end"] = pd.to_numeric(bed["end"], errors="coerce")
bed["num_exons"] = pd.to_numeric(bed["num_exons"], errors="coerce")

summary_rows = []

# Row 1: Gene spans
gene_spans = bed.groupby("base_gene_id")[["start", "end"]].agg(["min", "max"]).dropna()
gene_lengths = (gene_spans["end"]["max"] - gene_spans["start"]["min"]).sum()
summary_rows.append(["Protein coding genes", gene_spans.shape[0], gene_lengths, round(gene_lengths / genome_size * 100, 2)])

# Row 2: mRNA
mrna_spans = bed.dropna(subset=["transcript_id", "start", "end"])[["transcript_id", "start", "end"]].drop_duplicates()
mrna_spans["length"] = mrna_spans["end"] - mrna_spans["start"]
summary_rows.append(["mRNA", mrna_spans.shape[0], mrna_spans["length"].sum(), round(mrna_spans["length"].sum() / genome_size * 100, 2)])

# Row 3: CDS
cds_blocks = bed.dropna(subset=["cds_start", "cds_end"])[["cds_start", "cds_end"]].drop_duplicates()
cds_blocks["length"] = cds_blocks["cds_end"] - cds_blocks["cds_start"]
summary_rows.append(["CDS", cds_blocks.shape[0], cds_blocks["length"].sum(), round(cds_blocks["length"].sum() / genome_size * 100, 2)])

# Row 4: Exons
bed["exon_length_total"] = bed["exon_sizes"].apply(lambda val: sum(int(x) for x in str(val).strip(",").split(",") if x) if pd.notna(val) else 0)
summary_rows.append(["Exons", int(bed["num_exons"].dropna().sum()), bed["exon_length_total"].sum(), round(bed["exon_length_total"].sum() / genome_size * 100, 2)])

# Functional annotations
grouped = bed.groupby("base_gene_id", group_keys=False)

def has_nonhypo(descs): return any(pd.notna(d) and "hypothetical protein" not in d.lower() for d in descs)
def has_any(col): return any(pd.notna(col))

categories = {
    "Function prediction": ("description", has_nonhypo),
    "Pfam domains": ("pfam", has_any),
    "KOG assignments": ("kog", has_any),
    "KEGG pathway mapping": ("KEGG_Pathway", has_any),
    "GO term annotation": ("go_terms", has_any),
    "Product descriptions": ("description", has_any)
}

cds_subset = bed.dropna(subset=["cds_start", "cds_end"]).copy()
cds_subset["cds_length"] = cds_subset["cds_end"] - cds_subset["cds_start"]

for label, (field, checker) in categories.items():
    passing = grouped[field].apply(checker)
    ids = passing[passing].index
    subset = cds_subset[cds_subset["base_gene_id"].isin(ids)]
    span = subset[["transcript_id", "cds_start", "cds_end"]].drop_duplicates()
    span["length"] = span["cds_end"] - span["cds_start"]
    total_bp = span["length"].sum()
    summary_rows.append([label, len(ids), total_bp, round(total_bp / genome_size * 100, 2)])

summary = pd.DataFrame(summary_rows, columns=["Category", "Feature Count", "Total Length (bp)", "% of Genome"])
print(f"\nGenome size used: {genome_size:,} bp")
print("\nAnnotation Summary by Genome Occupancy:\n")
print(summary.to_string(index=False))


