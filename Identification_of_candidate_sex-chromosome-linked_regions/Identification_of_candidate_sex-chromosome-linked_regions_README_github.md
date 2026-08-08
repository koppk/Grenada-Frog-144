# Identification of candidate sex-chromosome-linked regions

Scripts for the identification of candidate Z- and W-sex-chromosome-linked
regions in the *Pristimantis euphronides* genome assembly from a single
female individual.

Corresponds to Additional file 2, section **"Identification of candidate
sex-chromosome-linked regions"** and Additional file 3 (Supplementary
Results), same section.

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Directory structure

Scripts are collected in a flat `scripts/` directory. For the full
directory structure with subdirectories by analysis stage, see the
Zenodo README
([doi:10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Scripts

Listed in AF2 paragraph order.

### Z-candidate identification (AF2 paras 3–8)

| Script | AF2 | Description |
|--------|-----|-------------|
| `z_candidate_workflow.sh` | 3–4 | Master workflow (layers 1–4) |
| `identify_z_candidate_regions.py` | 4 | Scaffold screening (binomial test, Bonferroni) |
| `detect_boundaries.py` | 4 | Boundary detection (binary segmentation, BIC) |
| `haplotype_identity.sh` | 5 | Inter-haplotype concordance (minimap2 asm5) |
| `mann_whitney_identity.py` | 5 | Exact Mann-Whitney U test |
| `scaffold_repeat_analysis.sh` | 6 | Scaffold repeat characterisation (WM + RM) |
| `parse_scaffold_repeats.py` | 6 | Parse repeats into per-Mb bins |
| `mann_whitney_repeats.py` | 6 | Mann-Whitney U tests for repeat density |
| `scaffold_gene_density.py` | 7 | Per-Mb gene density + gene body length |
| `classify_gene_coverage.sh` | 8 | Per-gene coverage classification (4 classes) |

### W-candidate identification (AF2 paras 2, 9–16)

| Script | AF2 | Description |
|--------|-----|-------------|
| `prepare_contig_sets.sh` | 2 | Partition contigs from AGP; SeqKit stats; WM mk_counts |
| `map_reads_primary_assembly_coverage.sh` | 9 | Per-contig coverage; coverage classification |
| `run_windowmasker_two_sets.sh` | 10 | WindowMasker on placed and unplaced contig sets |
| `parse_windowmasker_results.py` | 10 | Parse WM per-contig fractions |
| `run_repeatmasker_2sets.sh` | 10 | RepeatMasker on placed and unplaced contig sets |
| `parse_repeatmasker_two_sets.sh` | 10 | Parse RM per-contig breakdown |
| `primary_assembly_bam_qc.sh` | 10 | Per-contig QC metrics |
| `run_assembly_qc_pipeline.sh` | 10 | Wrapper calling map_reads + bam_qc |
| `run_miniprot_both_sets.sh` | 11 | Miniprot annotation |
| `gametolog_discovery_hanno7.sh` | 11–13 | HANNO gene partitioning + gametolog tiering |
| `gametolog_blastn.sh` | 13 | BLASTn pairwise alignment of pairs |
| `plot_z_regions_banding.py` | 14 | Scaffold-level compositional banding + segment detection |
| `compute_gene_coverage_by_zone.py` | 14 | Gene coverage classification by compositional zone |
| `plot_tier1a_gametolog_banding.py` | 14 | Tier 1a contig pair banding |
| `plot_tier1a_gametolog_synteny.py` | 14 | Tier 1a synteny diagrams |
| `plot_w_genic_banding.py` | 14 | W-genic contig banding |
| `plot_w_het_banding.py` | 14 | W-heterochromatin contig banding |
| `build_contig_master_table.py` | 15–16 | Master table + W-classification |

### W-exclusive gene search (AF2 paras 18–19)

| Script | AF2 | Description |
|--------|-----|-------------|
| `find_W_exclusive_genes.sh` | 18 | Name cross-reference + BLASTn vs placed scaffolds |
| `screen_pfam_unnamed_w_exclusive.sh` | 19 | Pfam domain screen |
| `blast_unnamed_w_exclusive_vs_nt.sh` | 19 | BLASTn vs NCBI nt |

### DM domain search (AF2 para 17)

| Script | AF2 | Description |
|--------|-----|-------------|
| `search_DM_domain.sh` | 17 | tBLASTn of Dm-W DM domain against full assembly |

### ZW ideogram (AF2 para 1)

| Script | AF2 | Description |
|--------|-----|-------------|
| `w_ideogram_schmid.py` | 1 | ZW ideogram from Schmid et al. (2002) cytogenetic data (Figure SR21) |

## Prerequisites

- Python 3.8+ with: matplotlib
- minimap2 v2.28, samtools, mosdepth v0.3.11, bedtools
- BLAST+ (makeblastdb, blastn, tblastn)
- WindowMasker, RepeatMasker v4.2.3 with Dfam 3.9 (Anura)
- miniprot, SeqKit
- HANNO v0.4 annotation outputs (used as input; pipeline not re-run here)

## Author

Kopp K, Pristimantis euphronides genome project
