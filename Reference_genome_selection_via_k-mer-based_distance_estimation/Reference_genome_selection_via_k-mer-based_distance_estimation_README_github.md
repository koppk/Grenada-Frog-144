# Reference genome selection via k-mer-based distance estimation

Pairwise genomic distance estimation between *Pristimantis euphronides*
assemblies and 21 publicly available Hyloidea genome assemblies using
Mash v2.3.

Corresponds to Additional file 2, section **"Reference genome selection
via k-mer-based distance estimation"** and Additional file 3
(Supplementary Results), same section.

Input data, output files, and figures are deposited at Zenodo
([doi:10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Directory structure

```
Reference_genome_selection_via_k-mer-based_distance_estimation/
├── input/
│   └── label_map.tsv
└── scripts/
    └── plot_MashDistance_Heatmap.py
```

## Workflow

```bash
mash sketch *.fna.gz *.fasta.gz
mash triangle *.msh > ref_triangle.tab
python3 plot_MashDistance_Heatmap.py ref_triangle.tab label_map.tsv
```

## Scripts

| Script | Description |
|--------|-------------|
| `plot_MashDistance_Heatmap.py` | 24x24 lower-triangle heatmap (Figure SR3) from `mash triangle` Phylip output |

## Prerequisites

- Mash v2.3
- Python 3.x with numpy, matplotlib, seaborn

## Author

Kopp K, Pristimantis euphronides genome project
