# Taxonomic classification of reads

Script for taxonomic classification of FAST basecalled Nanopore reads to
detect potential non-host sequences and assess taxonomic composition using
Kraken2 v2.1.2 with a custom database at a confidence threshold of 0.2.

Corresponds to Additional file 2, section **"Taxonomic classification of
reads"** and Additional file 3 (Supplementary Results), same section.

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)).

## Custom database

`DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core` was built by combining
15 chromosome-level Hyloidea genome assemblies (nine species, five families;
downloaded from NCBI Genome database, December 12, 2024) with the Kraken2
PlusPF reference database (release date: Dec 28, 2024;
https://benlangmead.github.io/aws-indexes/k2).

### Hyloidea assemblies included in the database

| Species | Family | Accession |
|---|---|---|
| *Eleutherodactylus coqui* | Eleutherodactylidae | GCF_035609145.1 |
| *Eleutherodactylus coqui* | Eleutherodactylidae | GCA_035609135.1 |
| *Eleutherodactylus coqui* | Eleutherodactylidae | GCA_019857665.1 |
| *Engystomops pustulosus* | Leptodactylidae | GCA_040894005.1 |
| *Engystomops pustulosus* | Leptodactylidae | GCA_040894015.1 |
| *Engystomops pustulosus* | Leptodactylidae | GCA_019512145.1 |
| *Leptodactylus fuscus* | Leptodactylidae | GCA_031893025.1 |
| *Leptodactylus fuscus* | Leptodactylidae | GCA_031893055.1 |
| *Dendropsophus ebraccatus* | Hylidae | GCF_027789765.1 |
| *Hyla sarda* | Hylidae | GCF_029499605.1 |
| *Bufo bufo* | Bufonidae | GCF_905171765.1 |
| *Bufo gargarizans* | Bufonidae | GCF_014858855.1 |
| *Bufotes viridis* | Bufonidae | GCA_033119425.1 |
| *Bufotes viridis* | Bufonidae | GCA_037900795.1 |
| *Ranitomeya imitator* | Dendrobatidae | GCF_032444005.1 |

## Commands

### 1. Database build

```bash
# Add each Hyloidea assembly to the database library
for i in *.fna
do
  kraken2-build --add-to-library $i \
    --db DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core \
    --threads 22
done

# Build the database
kraken2-build --build \
  --db DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core \
  --threads 22

# Inspect the database
kraken2-inspect \
  --db DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core \
  --threads 22 \
  > DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_Inspect_report.txt
```

### 2. Classification

```bash
kraken2 \
  --db DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core \
  --threads 22 \
  --confidence 0.2 \
  --output DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_classification_0.2 \
  --report DB_Hyloidea_arc_bac_plas_vir_fung_prot_UniVec_Core_report_0.2 \
  --use-names \
  --gzip-compressed GrenadaFrog144_ONT_ALL.fastq.gz
```

## Prerequisites

- Kraken2 v2.1.2 ([Wood et al., 2019](https://doi.org/10.1186/s13059-019-1891-0))

## Author

Kopp K, Pristimantis euphronides genome project
