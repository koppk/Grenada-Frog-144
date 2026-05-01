# Alignment with previously published Pristimantis euphronides and Eleutherodactylus johnstonei sequences

Script for molecular species identification of the sequenced individual
("GrenadaFrog144") by alignment against publicly available *P. euphronides*
GenBank sequences and differential comparison against *E. johnstonei*.

Corresponds to Additional file 2, section **"Alignment with previously
published Pristimantis euphronides and Eleutherodactylus johnstonei
sequences"** and Additional file 3 (Supplementary Results), same section.

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)).

## Workflow

| Step | Description | Method |
|------|-------------|--------|
| 1–5 | Minimap2 alignment to assembly, contig extraction, subsequence extraction | `01_align_published_sequences.sh` (local) |
| 6 | Pairwise BLASTn: contigs vs *P. euphronides* references | NCBI web BLASTn ("Align two sequences") |
| 7 | BLASTn: aligned portions vs *E. johnstonei* (taxid:350008) | NCBI web BLASTn (organism-restricted) |

Steps 6–7 were performed via the NCBI BLASTn web interface
(BLAST+ [Camacho, 2009](https://doi.org/10.1186/1471-2105-10-421)).
Alignment result files are at Zenodo.

## Script

| Script | Description |
|--------|-------------|
| `01_align_published_sequences.sh` | Minimap2 alignment, contig extraction, subsequence extraction (steps 1–5); documents NCBI web BLASTn steps |

### Prerequisites

- minimap2 v2.28 ([Li, 2018](https://doi.org/10.1093/bioinformatics/bty191))
- samtools v1.19.2 ([Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008))
- seqkit ([Shen et al., 2024](https://doi.org/10.1002/imt2.191))

### Execution

```bash
cd Alignment_with_previously_published_Pristimantis_euphronides_and_Eleutherodactylus_johnstonei_sequences
bash scripts/01_align_published_sequences.sh
# Then perform NCBI web BLASTn steps as documented in the script
```

## Author

Kopp K, Pristimantis euphronides genome project
