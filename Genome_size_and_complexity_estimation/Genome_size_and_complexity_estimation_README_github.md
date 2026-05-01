# Genome size and complexity estimation

Reference-free genome size estimation, expected sequencing depth
calculation, and read-based G+C content estimation for *Pristimantis euphronides*.

Corresponds to Additional file 2, section **"Genome size and complexity
estimation"** and Additional file 3 (Supplementary Results), same section
(Table SR2; Figure SR1).

## Directory structure

```
Genome_size_and_complexity_estimation/
└── Genome_size_and_complexity_estimation_README_github.md
```

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)).

## Commands

### 1. K-mer counting (Jellyfish v2.3.1)

```bash
pigz -dc -p 16 GrenadaFrog144_ONT_HAC_all.fastq.gz > GrenadaFrog144_ONT_HAC_all.fastq
jellyfish count -C -m 21 -s 1000000000 -t 48 GrenadaFrog144_ONT_HAC_all.fastq -o reads.jf
jellyfish histo -t 10 reads.jf > reads.histo
```

### 2. Genome size estimation (GenomeScope 2.0)

The histogram was submitted to the GenomeScope 2.0 web interface at
http://genomescope.org/genomescope2.0/ with parameters: k-mer length = 21,
ploidy = 2, max k-mer coverage = -1, average k-mer coverage for polyploid
genome = -1. The complete results page is archived in the Zenodo deposit
as `output/GenomeScope.html`.

### 3. Expected sequencing depth

```bash
seqkit stats GrenadaFrog144_ONT_HAC_all.fastq.gz -a
# 55595176615 / 1617259408 = 34.4×
```

### 4. Read-based G+C content (seqtk v1.4)

```bash
seqtk comp GrenadaFrog144_ONT_HAC_all.fastq.gz > Seqtk_Stats.GrenadaFrog144_ONT_HAC_all.fastq.gz.txt
awk -F '\t' '{ sum += $4 } END { print sum }' Seqtk_Stats.GrenadaFrog144_ONT_HAC_all.fastq.gz.txt  # C: 12056034465
awk -F '\t' '{ sum += $5 } END { print sum }' Seqtk_Stats.GrenadaFrog144_ONT_HAC_all.fastq.gz.txt  # G: 11980108195
# (C+G) / total_bases * 100 = 43.23%
```

## Prerequisites

- Jellyfish v2.3.1
- pigz
- SeqKit v2.10.0
- seqtk v1.4
- GenomeScope 2.0 (web interface)

## Author

Kopp K, Pristimantis euphronides genome project
