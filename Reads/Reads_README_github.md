# Reads

Combined Oxford Nanopore Technologies (ONT) sequencing read datasets for
the *Pristimantis euphronides* genome project ("GrenadaFrog144"), a
single female individual sequenced across 12 flow cells using two
MinION Mk1B devices with R10.4.1 flow cells and Kit 14 (SQK-LSK114)
chemistry.

Corresponds to Additional file 2, sections **"Library preparation and
genome sequencing"** and **"Read quality analysis"**, and
Additional file 3, sections **"Read quality analysis"** and
**"Taxonomic classification of reads"**.

## Basecalling

All 12 sequencing runs were basecalled using Dorado Basecall Server
v7.4.12 orchestrated by MinKNOW v24.06.16 on GPU-equipped laptops
(NVIDIA GeForce RTX 3060).

**FAST basecalling** was performed in real-time during each sequencing
run using model `dna_r10.4.1_e8.2_400bps_fast@v4.3.0`. Reads with a
mean quality score below Q8 were filtered out by MinKNOW.

**HAC basecalling** was performed offline after each run using model
`dna_r10.4.1_e8.2_400bps_hac@v4.3.0`. Reads with a mean quality score
below Q9 were filtered out by MinKNOW.

Adapter sequences and ONT-specific artificial sequences were
automatically trimmed by MinKNOW during both basecalling modes.

## Combining reads

The 12 per-run FAST basecalled FASTQ files were concatenated with `cat`
into `GrenadaFrog144_ONT_ALL.fastq.gz`. The 12 per-run HAC basecalled
FASTQ files were concatenated likewise into
`GrenadaFrog144_ONT_HAC_all.fastq.gz`.

## Data availability

Read data were deposited in the European Nucleotide Archive (ENA) under
BioProject PRJEB89028
(https://www.ebi.ac.uk/ena/browser/view/PRJEB89028) and are also
available via NCBI
(https://www.ncbi.nlm.nih.gov/bioproject/PRJEB89028/).
The combined read files and all supplementary
data are archived at Zenodo
(doi: [10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Author

Kopp K, Pristimantis euphronides genome project
