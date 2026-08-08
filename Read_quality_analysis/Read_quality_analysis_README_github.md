# Read quality analysis

Scripts and commands for read quality assessment of Oxford Nanopore
sequencing data from *Pristimantis euphronides* (GrenadaFrog144).

Corresponds to Additional file 2, section **"Read quality analysis"**
and Additional file 3 (Supplementary Results), same section.

## Scripts

| Script | Description |
|--------|-------------|
| `compute_arithmetic_mean_phred.sh` | Per-read arithmetic mean and read-length-weighted mean of Phred scores |

## Commands

```bash
# NanoComp v1.24.2: FAST vs HAC comparison (produces NanoStats.txt and plots)
NanoComp --fastq GrenadaFrog144_ONT_ALL.fastq.gz \
         GrenadaFrog144_ONT_HAC_all.fastq.gz \
         --names GrenadaFrog144_ONT_FAST GrenadaFrog144_ONT_HAC \
         --format png --threads 14

# SeqKit v2.10.0: per-base quality thresholds for HAC reads (Table SM1)
seqkit stats -a GrenadaFrog144_ONT_HAC_all.fastq.gz

# Arithmetic mean of Phred scores
bash compute_arithmetic_mean_phred.sh \
    GrenadaFrog144_ONT_ALL.fastq.gz \
    GrenadaFrog144_ONT_HAC_all.fastq.gz
```

## Prerequisites

- NanoComp v1.24.2 (includes NanoStat v1.6.0 from the NanoPack suite)
- SeqKit v2.10.0
- pigz
- python3

## Data availability

Read data were deposited in the European Nucleotide Archive (ENA) under
BioProject PRJEB89028
(https://www.ebi.ac.uk/ena/browser/view/PRJEB89028) and are also
available via NCBI
(https://www.ncbi.nlm.nih.gov/bioproject/PRJEB89028/).
All supplementary data are archived at Zenodo
(doi: [10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Author

Kopp K, Pristimantis euphronides genome project
