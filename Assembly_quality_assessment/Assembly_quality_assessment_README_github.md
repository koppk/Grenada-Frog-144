# Assembly quality assessment

Scripts for the quality assessment of the *Pristimantis euphronides*
genome assembly.

Corresponds to Additional file 2, section **"Assembly quality assessment"**
and Additional file 3 (Supplementary Results), same section.

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Directory structure

On the server, the analysis components resided at
`/data/GrenadaFrog144/Inspector/` (Inspector),
`/data/GrenadaFrog144/Merqury/` (Merqury),
`/data/GrenadaFrog144/assembly_qc/` (per-contig QC),
`/data/software/HANNO/` (gene-space completeness), and
`/data/GrenadaFrog144/Hyloidea_Proteins/` (HANNO evidence).
For GitHub they are grouped under `Assembly_quality_assessment/`:

```
Assembly_quality_assessment/
├── Inspector/                   Inspector v1.2 (1 script)
├── Merqury/                     Merqury v1.3 (1 script)
├── HANNO/                       Gene-space completeness (2 pipeline scripts)
├── Hyloidea_Proteins/           HANNO evidence inputs (1 script)
└── Inter_workflow_concordance_pre-scaffolding/
                                 D-Genies v1.5.0 statistics (1 script)
```

## Scripts

### Inspector

| Script | Description |
|--------|-------------|
| `run_inspector.sh` | Reference-free assembly evaluation (--datatype nanopore) |

### Merqury

| Script | Description |
|--------|-------------|
| `run_merqury.sh` | K-mer-based consensus QV and completeness (k=21, non-trio) |

### HANNO

| Script | Description |
|--------|-------------|
| `WithmRNA_...final_medaka_polished_assembly_consensus.fasta.sh` | HANNO v0.4 for Workflow 1 (Medaka-polished Flye) |
| `WithmRNA_...grenada-frog-HK.ctg.fa.gz.sh` | HANNO v0.4 for Workflow 2 (Wtdbg v1.0) |

### HANNO evidence preparation

| Script | Description |
|--------|-------------|
| `download_and_combine_hyloidea_evidence.sh` | Download from NCBI RefSeq FTP, prefix GTF identifiers, concatenate |

### Inter-workflow concordance (pre-scaffolding)

| Script | Description |
|--------|-------------|
| `compute_inter_workflow_dgenies_stats.sh` | Compute concordance statistics from D-Genies v1.5.0 output (PAF, association table, index files). Reproduces pre-scaffolding numbers in AF3 |

Pre-scaffolding assemblies from Workflow 1 (Flye/Medaka) and Workflow 2
(Wtdbg) were aligned against each other in both directions using D-Genies
v1.5.0 with embedded Minimap2 and default parameters. No custom alignment
script; performed via the D-Genies web interface. D-Genies output files
(PAF, dotplots, association tables) are deposited at Zenodo.

### Cross-section scripts

Per-contig QC outputs in `assembly_qc/` are produced by scripts in
`Identification_of_candidate_sex-chromosome-linked_regions/W-chr-workflow/scripts/`:

| Script | Description |
|--------|-------------|
| `map_reads_primary_assembly_coverage.sh` | Read mapping + per-contig coverage |
| `primary_assembly_bam_qc.sh` | Coverage CV, MAPQ, base quality, variant density, soft-clipping, alignment identity |
| `run_assembly_qc_pipeline.sh` | Wrapper calling both scripts above |

## Prerequisites

- Inspector v1.2
- Merqury v1.3, meryl v1.4.1
- HANNO v0.4 with: miniprot v0.13, minimap2 v2.28, StringTie v2.2.3, TACO v0.6.2, TransDecoder v5.7.1, LAST v1595, BUSCO v3.0.2, eggNOG-mapper
- mosdepth v0.3.11, samtools v1.19.2, bcftools, GNU parallel, SeqKit, NanoStat v1.6.0
- D-Genies v1.5.0 (web interface; Cabanettes & Klopp, 2018)

## Author

Kopp K, Pristimantis euphronides genome project
