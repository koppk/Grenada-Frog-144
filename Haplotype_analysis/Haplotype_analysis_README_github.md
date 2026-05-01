# Haplotype analysis

Haplotype dual assembly, haplotype alignment and structural
comparison, reference genome similarity assessment, reference-guided
scaffolding of haplotype assemblies, and scaffolded haplotype alignment
and structural comparison for the *Pristimantis euphronides* genome.

Corresponds to Additional file 4 (Supplementary Methods — haplotypes)
and Additional file 5 (Supplementary Results — haplotypes).

HapDup v0.12 dual haplotype assemblies, WindowMasker soft-masking of
*Eleutherodactylus coqui* (*E. coqui*) haplotype 2 (GCA_035609135.1), RagTag v2.1.0 scaffolding,
D-Genies v1.5.0 pre- and post-scaffolding comparisons, scaffolded
haplotype assemblies, HapDup intermediate files, and computed
statistics are deposited at
[Zenodo](https://doi.org/10.5281/zenodo.15298547).

## Directory structure

```
Haplotype_analysis/
└── scripts/
    ├── run_minimap2_HapDup_prep.sh
    ├── run_docker_HapDup.sh
    ├── compute_hapdup_phasing_stats.sh
    ├── compute_hapdup_bam_stats.sh
    ├── compute_assembly_stats.sh
    └── compute_haplotype_alignment_stats.sh
```

## Scripts

| Script | Description |
|--------|-------------|
| `run_minimap2_HapDup_prep.sh` | Map ONT HAC reads to assembly for HapDup input (Minimap2 v2.28, SAMtools v1.19.2) |
| `run_docker_HapDup.sh` | Run HapDup v0.12 via Docker (`mkolmogo/hapdup:0.12`) |
| `compute_hapdup_phasing_stats.sh` | Compute variant, phasing, and phase block statistics from PEPPER/Margin output (AF5 Section 1) |
| `compute_hapdup_bam_stats.sh` | Compute read mapping, haplotagging, and structural variant statistics from HapDup BAM files (AF5 Section 1) |
| `compute_assembly_stats.sh` | Compute assembly statistics (SeqKit) and AGP-derived gap/unplaced counts (AF5 Table HR1) |
| `compute_haplotype_alignment_stats.sh` | Compute alignment statistics from D-Genies PAF and association table (AF5 Table HR2) |

## Prerequisites

- Minimap2 v2.28
- SAMtools v1.19.2
- HapDup v0.12 (Docker: `mkolmogo/hapdup:0.12`)
- WindowMasker v1.0.0
- RagTag v2.1.0
- SeqKit v2.10.0
- D-Genies v1.5.0

## Author

Kopp K, Pristimantis euphronides genome project
