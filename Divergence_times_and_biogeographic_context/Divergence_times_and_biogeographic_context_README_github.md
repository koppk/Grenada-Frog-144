# Divergence times and biogeographic context

Divergence time extraction and biogeographic analysis of *Pristimantis
euphronides* using the Portik et al. (2023) time-calibrated anuran phylogeny
(5,242 tips).

Corresponds to Additional file 2, section **"Divergence times and
biogeographic context"** and Additional file 3 (Supplementary Results),
same section.


## Server location

On the project server (koppsrv@sd-182851), this analysis resides under
`/data/GrenadaFrog144/DivTimeGeo/analysis/`.


## Overview

This analysis quantifies the phylogenetic isolation of *P. euphronides* by
computing tMRCA and topological distance to all anuran species in the Portik
et al. (2023) tree, enriching the dataset with geographic occurrence data from
AmphibiaWeb, ranking available chromosome-level reference genomes by
phylogenetic proximity, and computing per-species type-locality distances to
Pointe Salines (Grenada) for the 12 closest mainland *Pristimantis* congeners.

Key findings:
- The nearest relative with a chromosome-level genome is *Eleutherodactylus
  coqui* (tMRCA 51.85 Ma), while the closest relative overall is *P. shrevei*
  (tMRCA 3.15 Ma) — a 16.5-fold gap.
- All five anuran species co-occurring on Grenada diverged from *P. euphronides*
  more than 50 Ma ago.
- Five of six Caribbean *Pristimantis* are single-island endemics.
- The 12 closest mainland *Pristimantis* congeners (tMRCA ~18.83 million years)
  occur in the Andean highlands of Colombia and Venezuela, with type-locality
  distances to Pointe Salines ranging from 1,267 to 1,666 km (Additional
  file 3: Table SR36).


## Quick start

```bash
cd Divergence_times_and_biogeographic_context/
bash main_divtime_geo_genomes.sh
```

The script runs four Python steps in sequence (tMRCA computation → geographic
enrichment → reference selection summary → type-locality distance computation)
and writes all results to `output/`.


## Folder structure

```
Divergence_times_and_biogeographic_context/
├── input/                                       # External reference data
│   ├── Supplementary_File_S3_time_tree.tre      # Portik et al. 2023 phylogeny
│   ├── Anura_chromosome-level_genomes_unique.tsv # NCBI Genome query result
│   ├── amphib_names.txt                         # AmphibiaWeb species export
│   ├── countries.csv                            # Country centroid coordinates
│   ├── 12_Pristimantis_amphibiansoftheworld.txt # ASW type-locality strings
│   └── 12_Pristimantis_amphibiansoftheworld_WithGoogleCoords.txt # ASW + Google Maps URLs
├── scripts/
│   ├── build_tree_metrics_chromlevel.py          # tMRCA + topological distance
│   ├── add_isocc_and_distances.py                # Geographic enrichment
│   ├── create_reference_selection_summary.py     # Ranked summary + gap stats
│   └── compute_typelocality_distances.py         # Per-species type-locality distances
├── output/                                      # Generated results
│   ├── tmrca_topo_chromlevel.tsv
│   ├── tmrca_topo_chromlevel_with_geo.tsv
│   ├── chromosome_level_ranked_by_tMRCA.tsv
│   ├── reference_selection_summary.txt
│   └── 12_Pristimantis_typelocality_distances.tsv
├── main_divtime_geo_genomes.sh                  # Entry point
└── Divergence_times_and_biogeographic_context_README_github.md
```


## Dependencies

- Python 3.11
- Biopython (Bio.Phylo)
- pandas


## Data availability

Input data and pre-computed output files are archived on Zenodo:
[doi:10.5281/zenodo.15298546](https://doi.org/10.5281/15298546)

The Portik et al. (2023) time tree is available as Supplementary File S3 of:
Portik DM, Streicher JW, Wiens JJ. Frog phylogeny: a time-calibrated,
species-level tree based on hundreds of loci and 5,242 species. Mol Phylogenet
Evol. 2023;188:107907. Available from:
https://github.com/nhm-herpetology/frog-phylogeny/blob/main/Supplementary_File_S3_time_tree.tre.
Accessed January 2026.

## Author

Kopp K, Pristimantis euphronides genome project
