# Reference-guided scaffolding of genome assemblies

Scripts for the reference-guided scaffolding of the *Pristimantis
euphronides* Medaka-polished Flye assembly against the *Eleutherodactylus
coqui* chromosome-level reference genome using RagTag v2.1.0.

Corresponds to **Additional file 2** (Methods) and **Additional file 3**
(Results), section "Reference-guided scaffolding of genome assemblies",
and the main manuscript section of the same name in Methods, Results,
and Discussion.

Input and output data files are deposited at Zenodo
(doi: [10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)).
The final scaffolded genome
(`Pristimantis_euphronides.genome.fasta.gz`) is also deposited at
NCBI/ENA under accession
[GCA_965278355.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_965278355.2/).

## Directory structure

```
Reference-guided_scaffolding_of_genome_assemblies/
├── Workflow_1/
│   └── scripts/
│       ├── postprocess_ragtag.sh
│       ├── split_contig_length_and_placement.sh
│       └── compute_contig_length_and_placement_stats.py
└── Workflow_2/                         (data on Zenodo only)
```

Workflow 1 inputs (WindowMasker-masked *E. coqui* reference, RagTag
output, Medaka-polished Flye assembly), the renamed and filtered
deposited genome, and Figure SR4 are all on Zenodo. Workflow 2
(`grenada-frog-HK.scf.fa.gz`, Wtdbg v1.0 + Ragout v2.3, produced by
co-author H. Kuhl) is on Zenodo only.

## Scripts

Scripts are numbered by execution order. Each step depends on output
from preceding steps.

| Script | Description | AF2 section |
|--------|-------------|-------------|
| `postprocess_ragtag.sh` | Rename RagTag scaffolds (NC_/NW_ accessions → scaffold_1 to scaffold_31 in *E. coqui* chromosome size order), update AGP, filter contigs < 200 bp from FASTA | Step 3 |
| `split_contig_length_and_placement.sh` | Split contigs from the renamed AGP into three `.len` files (scaffolds 1-13, scaffolds 14-31, unplaced); contigs < 200 bp excluded to match the deposited assembly | Step 4 |
| `compute_contig_length_and_placement_stats.py` | Compute per-group placement summary, length distribution statistics (Table SR6), and the log-scale boxplot (Figure SR4) | Step 4 |

### Prerequisites

- RagTag v2.1.0 ([Alonge et al., 2022](https://doi.org/10.1186/s13059-022-02823-7))
- WindowMasker v1.0.0 ([Morgulis et al., 2006](https://doi.org/10.1093/bioinformatics/bti774))
- SeqKit v2 ([Shen et al., 2016](https://doi.org/10.1371/journal.pone.0163962))
- Seqtk v1.4-r122 (https://github.com/lh3/seqtk)
- Python 3.x with numpy, matplotlib

### Execution

```bash
# Step 1: Soft-mask the E. coqui haplotype 1 reference genome
windowmasker -mk_counts \
    -in GCF_035609145.1_aEleCoq1.hap1_genomic.fna \
    -out GCF_035609145.1_aEleCoq1.hap1_genomic.wmstat

windowmasker -ustat GCF_035609145.1_aEleCoq1.hap1_genomic.wmstat \
    -in GCF_035609145.1_aEleCoq1.hap1_genomic.fna \
    -outfmt fasta \
    -out GCF_035609145.1_aEleCoq1.hap1_genomic.softmasked.fasta

# Step 2: RagTag scaffolding
ragtag.py scaffold \
    GCF_035609145.1_aEleCoq1.hap1_genomic.softmasked.fasta \
    final_medaka_polished_assembly_consensus.fasta \
    -o RagTag_EleCoq1.hap1.softmasked_finalmedaka_polished_assembly \
    -t 24

# Step 3: Postprocess (rename scaffolds, filter < 200 bp)
bash scripts/postprocess_ragtag.sh ragtag_output/ output/

# Step 4: Contig length and placement analysis
bash scripts/split_contig_length_and_placement.sh
python3 scripts/compute_contig_length_and_placement_stats.py
```

## Author

Kopp K, Pristimantis euphronides genome project
