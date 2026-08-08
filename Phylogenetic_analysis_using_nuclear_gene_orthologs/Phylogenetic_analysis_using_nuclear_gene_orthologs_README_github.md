# Phylogenetic analysis using nuclear gene orthologs

Scripts for the phylogenetic placement of *Pristimantis euphronides*
within the anuran superfamily Hyloidea using nuclear gene orthologs.
The pipeline covers CDS extraction from the HANNO annotation,
Hyloidea homolog retrieval from GenBank, *Rana temporaria* outgroup
extraction, alignment, quality filtering, taxonomy assignment,
supermatrix construction, and IQ-TREE inference.

Corresponds to Additional file 2 (Supplementary Methods) and
Additional file 3 (Supplementary Results), section
**"Phylogenetic analysis using nuclear gene orthologs"**.

Input data, intermediate outputs, IQ-TREE results, and full-size
tree figure PDFs (Additional file 11: Figures 1–12; Additional
file 3: Figures SR30, SR31) are deposited at Zenodo
(doi: [10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Scripts

```
scripts/
├── extract_selected_gene_coding_sequences.sh
├── fetch_hyloidea_sequences.sh
├── extract_cds_gene_from_RanTemp_RefGenome.sh
├── concat_gene_alignments.sh
├── full_alignment_pipeline.sh
├── final_refined_mafft.sh
├── rename_headers_for_supermatrix.sh
├── get_taxids.sh
├── create_taxid_family_genus_species.sh
├── extract_complete_taxa_supermatrix.sh
├── build_5gene_supermatrix.sh
├── build_supermatrix_clean.sh
├── iqtree3_batch_gene_trees.sh
├── run_iqtree_5gene_supermatrix.curated.sh
├── run_iqtree_12gene_supermatrix.curated.sh
└── make_table.sh
```

| # | Script | AF2 ¶ | Description |
|---|--------|-------|-------------|
| 1 | `extract_selected_gene_coding_sequences.sh` | 1 | Extract 18 *P. euphronides* CDS from HANNO BESTMODELS with gene name synonym mapping |
| 2 | `fetch_hyloidea_sequences.sh` | 2 | Download Hyloidea homologs from NCBI GenBank per gene |
| 3 | `extract_cds_gene_from_RanTemp_RefGenome.sh` | 2 | Extract *R. temporaria* outgroup CDS from GCF_905171775.1 |
| 4 | `concat_gene_alignments.sh` | 3 | Concatenate *P. euphronides* + Hyloidea + *R. temporaria* per gene |
| 5 | `full_alignment_pipeline.sh` | 3 | MAFFT alignment, terminal gap trimming, MAFFT re-alignment |
| 6 | `final_refined_mafft.sh` | 3 | Final high-accuracy MAFFT re-alignment after Jalview curation |
| 7 | `rename_headers_for_supermatrix.sh` | 3 | Standardise headers to Species__ACC__|_Genus_|_Family format |
| 8 | `get_taxids.sh` | 3 | Retrieve NCBI TaxIDs for GenBank accessions |
| 9 | `create_taxid_family_genus_species.sh` | 3 | Assign taxonomy using taxonkit v0.19.0 (NCBI taxdump, 17 April 2025) |
| 10 | `extract_complete_taxa_supermatrix.sh` | 3 | Filter to species with all 12 genes for 12-gene supermatrix |
| 11 | `build_5gene_supermatrix.sh` | 3 | Construct 5-gene supermatrix with partition file |
| 12 | `build_supermatrix_clean.sh` | 3 | Construct 12-gene supermatrix with partition file |
| 13 | `iqtree3_batch_gene_trees.sh` | 4 | IQ-TREE on 12 individual gene alignments |
| 14 | `run_iqtree_5gene_supermatrix.curated.sh` | 4 | IQ-TREE on 5-gene supermatrix (MFP+MERGE, partitioned) |
| 15 | `run_iqtree_12gene_supermatrix.curated.sh` | 4 | IQ-TREE on 12-gene supermatrix (MFP+MERGE, partitioned) |
| 16 | `make_table.sh` | 3 | Fetch TaxID, scientific name, genus, and family for GenBank accessions using NCBI Entrez utilities |

## Prerequisites

- SeqKit
- bash 4+ (for associative arrays)
- MAFFT v7.475
- TrimAl v1.4
- IQ-TREE v2.2.0.3
- NCBI Entrez utilities
- taxonkit v0.19.0 with NCBI Taxonomy database (taxdump)
- Jalview v2.11.4.1 (manual curation steps)

## Author

Kopp K, Pristimantis euphronides genome project
