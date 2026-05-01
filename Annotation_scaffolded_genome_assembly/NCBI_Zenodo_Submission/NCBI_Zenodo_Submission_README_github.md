# NCBI Submitter Annotation for *Pristimantis euphronides* GCA_965278355.2

GitHub holds scripts and commands. Full file descriptions, gene numbering,
and Gnomon comparison: see Zenodo
([10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)).


## Assembly

| Field | Value |
|---|---|
| Species | *Pristimantis euphronides* |
| Assembly | aPriEup1.0 |
| GenBank | GCA_965278355.2 |
| BioProject | PRJEB89028 |
| BioSample | SAMEA118241337 |
| Specimen | GrenadaFrog144 (adult female) |
| Sequencing | ONT MinION |
| Pipeline | Flye v2.9.5 → Medaka v2.0.1 → RagTag v2.1.0 |


## Annotation summary

| Metric | Value |
|---|---|
| Protein-coding gene models | 28,071 |
| Genes with gene symbols | 16,891 (60.2%) |
| Locus tag prefix | PRIEUP |
| Locus tag range | PRIEUP\_00001 – PRIEUP\_28071 |


## Directory Structure

```
NCBI_Zenodo_Submission/
├── NCBI_Zenodo_Submission_README_github.md    ← this file
├── generate_seqid_map.py
├── reformat_hanno_to_ncbi.py
└── run_ncbi_reformat.sh
```


## Reformatting from HANNO output

Annotation files generated from HANNO BESTMODELS output using
`reformat_hanno_to_ncbi.py` (https://github.com/koppk/Grenada-Frog-144).
Scaffold-to-accession mapping derived from the NCBI assembly report
for GCA_965278355.2.


## Author

Kopp K, Pristimantis euphronides genome project

## Citation

Kopp K, Lownds S, Harrison B, Kuhl H, Moittie S. From Field to Genome:
Non-Lethal, On-Site Nanopore-Only Assembly of the Endangered Grenada Frog
(*Pristimantis euphronides*). 2026. *BMC Genomics* (submitted for review).

## License

CC BY 4.0
