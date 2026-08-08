# Evaluation of reference-guided scaffolding

Scripts for the three-tier evaluation of the *E. coqui* reference-guided
RagTag scaffolding of the *Pristimantis euphronides* genome assembly, and
inter-workflow scaffold concordance assessment.

Corresponds to Additional file 2, section **"Evaluation of reference-guided
scaffolding"** and Additional file 3 (Supplementary Results), same section.

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Scripts

| Script | Description | Tier |
|--------|-------------|------|
| `01_four_species_synteny.py` | Compound scaffold identification, breakpoint detection, cross-validation, heatmaps, positional synteny plots | tier (ii) |
| `02_contig_allegiance.py` | Breakpoint classification as predicted fusions or predicted scaffolding artifacts | tier (iii) |
| `03_karyotype_schematic.py` | Karyotype reconstruction figure from tier (ii) + (iii) outputs | Figure |
| `compute_inter_workflow_dgenies_stats.sh` | Compute concordance statistics from D-Genies v1.5.0 output. Reproduces post-scaffolding inter-workflow numbers in AF3 | Inter-workflow |

Tier (i) (pre- and post-scaffolding dotplots vs *E. coqui*) was performed
with D-Genies v1.5.0; PAF files at Zenodo.

Inter-workflow scaffold concordance (post-scaffolding) was assessed by
aligning Workflow 2 (Ragout v2.3) against Workflow 1 (RagTag v2.1.0) in
both directions using D-Genies v1.5.0 with embedded Minimap2. No custom
alignment script; performed via the D-Genies web interface. D-Genies output
files (PAF, dotplots, association tables) are deposited at Zenodo.

### Prerequisites

- Python 3.8+ with: matplotlib
- D-Genies v1.5.0 ([Cabanettes & Klopp, 2018](https://doi.org/10.7717/peerj.4958))
  with Minimap2 v2.28 ([Li, 2018](https://doi.org/10.1093/bioinformatics/bty191))

### Execution

```bash
cd Evaluation_of_reference-guided_scaffolding

python3 scripts/01_four_species_synteny.py
python3 scripts/02_contig_allegiance.py
python3 scripts/03_karyotype_schematic.py

# Inter-workflow concordance statistics
bash scripts/compute_inter_workflow_dgenies_stats.sh \
  output/Inter_workflow_concordance_post-scaffolding/D-genies_post-scaffolding_Workflow_2_vs_Workflow_1

bash scripts/compute_inter_workflow_dgenies_stats.sh \
  output/Inter_workflow_concordance_post-scaffolding/D-genies_post-scaffolding_Workflow_1_vs_Workflow_2
```

## Output files

| File | Script | AF3 reference |
|------|--------|---------------|
| `Table_compound_chromosomes.tsv` | 01 | Table SR7 |
| `Table_karyotype_reconstruction.tsv` | 01 | Table SR8 |
| `Table_small_scaffolds.tsv` | 01 | Additional file 10 |
| `Fig_breakpoint_scaffold_{1,3,4,6}.*` | 01 | Figures SR9–SR12 |
| `Fig_heatmap_*.*` | 01 | Figures SR13–SR17 |
| `allegiance_breakpoint_results.tsv` | 02 | Table SR9 |
| `allegiance_full_scaffolds.*` | 02 | Figure SR18 |
| `allegiance_breakpoint_zooms.*` | 02 | Figure SR19 |
| `Fig_karyotype_schematic.*` | 03 | Figure SR8 |
| `Inter_workflow_concordance_post-scaffolding/*/stats_summary.txt` | compute_inter_workflow_dgenies_stats.sh | Figure SR20 paragraph |

## Author

Kopp K, Pristimantis euphronides genome project
