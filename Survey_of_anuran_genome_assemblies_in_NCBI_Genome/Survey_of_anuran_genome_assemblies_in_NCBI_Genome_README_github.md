# Survey of anuran genome assemblies in NCBI Genome

Scripts for a systematic survey of all publicly available anuran genome
assemblies in NCBI Genome, classifying sequencing technology, tissue
source, and taxonomic representation within Terrarana.

Corresponds to **Additional file 12** (Methods and Results) and the
main manuscript sections "Survey of anuran genome assemblies" in
Methods, Results, and Discussion.

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Scripts

Scripts are numbered by execution order (steps 1–4). Each step depends
on output from preceding steps.

| Script | Description | AF12 section |
|--------|-------------|--------------|
| `scrub_anuran_genomes_ncbi.sh` | Query NCBI Genome for all Anura (taxid 8342) assemblies via NCBI Datasets CLI | Methods para 1 |
| `analyze_anuran_genomes.sh` | Deduplicate GCA/GCF pairs, classify sequencing technology, assign family and superfamily via taxonkit | Methods paras 1–3 |
| `fetch_biosample_tissue.sh` | Fetch BioSample metadata via efetch in batches of 50; classify tissue following Taberlet et al. (1999) | Methods para 4 |
| `narrative_funnel.sh` | Generate narrative summary report with all counts and tables | Results, Tables 1–5 |

### Prerequisites

- NCBI Datasets CLI v16 (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
- NCBI Entrez Direct (esearch, efetch) ([Sayers et al., 2022](https://doi.org/10.1093/nar/gkab1112))
- taxonkit v0.19.0 ([Shen & Ren, 2021](https://doi.org/10.1016/j.jgg.2021.03.006))
- NCBI Taxonomy database (taxdump, April 2026)
- Python 3.11+ (xml.etree.ElementTree, standard library only)
- GNU awk, GNU coreutils

### Execution

```bash
# Step 1: Query NCBI for all Anura genome assemblies (~2 min)
bash scripts/scrub_anuran_genomes_ncbi.sh -o output

# Step 2: Deduplicate, classify tech, add taxonomy
bash scripts/analyze_anuran_genomes.sh -i output/anuran_genomes_ncbi.tsv -o output

# Step 3: Fetch BioSample tissue metadata and classify (~5–10 min)
bash scripts/fetch_biosample_tissue.sh -i output/anuran_genomes_ncbi.tsv -o output

# Step 4: Narrative funnel report
bash scripts/narrative_funnel.sh -o output
```

## Author

Kopp K, Pristimantis euphronides genome project
