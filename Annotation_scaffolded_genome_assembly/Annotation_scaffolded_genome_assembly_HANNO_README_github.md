# Annotation of Scaffolded Genome Assembly

HANNO v0.4 (Kuhl 2024) annotation of *Pristimantis euphronides* assembly
`aPriEup1.0` with protein and mRNA evidence.

Corresponds to Additional file 2, section **"Annotation of scaffolded
genome assembly"** and Additional file 3 (Supplementary Results), same
section.

GitHub holds scripts and commands. Full data, output files, and detailed
descriptions: see Zenodo
([10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)).


## Directory Structure

```
Annotation_scaffolded_genome_assembly/
├── Annotation_scaffolded_genome_assembly_HANNO_README_github.md    ← this file
│
├── HANNO_WithmRNA_HyloideaProteinSet_Pristimantis_euphronides_genome/
│   │   Workflow 1 (Flye, Medaka, RagTag; ENA/NCBI-submitted genome).
│   ├── HANNO-RUN/                     # HANNO bash script + log
│   └── Follow-up-scripts/             # Post-HANNO processing (§3)
│
├── HANNO_WithmRNA_HyloideaProteinSet_grenada-frog-HK.scf.genome/
│   │   Workflow 2 (wtdbg, Ragout). HANNO log only.
│   └── ...
│
└── NCBI_Zenodo_Submission/            # See NCBI_Zenodo_Submission_README_github.md
```


## 1. HANNO Run

Generated with HANNO v0.4 Perl wrapper (24 CPUs, 252 GB RAM):

```bash
scripts/HANNO.v0.4.pl -t 24 \
  -d WithmRNA_HyloideaProteinSet_HANNO-RUN-Pristimantis_euphronides.genome.fasta \
  -a ../../GrenadaFrog144/Pristimantis_euphronides.genome.fasta \
  -p ../../GrenadaFrog144/Hyloidea_Proteins/combined_Hyloidea_proteins.faa \
  -r ../../GrenadaFrog144/Hyloidea_Proteins/combined_rna.fna.gz \
  -b tetrapoda_odb10 -B 0 -P tetrapoda_odb10/refseq_db.faa -E 0
```

Paths must be relative to `/data/software/HANNO/`.


## 2. Script Modification

Default `stringtie --merge` replaced with `stringtie_merge_safe()` shell
function (4 parallel groups by scaffold size, unique MSTRG ID prefixes).
Reduces peak memory from full genome to largest group.

Run:

```bash
nohup ./WithmRNA_HyloideaProteinSet_HANNO-RUN-Pristimantis_euphronides.genome.sh \
  > ...log 2>&1 &
```


## 3. Post-HANNO Processing

### 3.1 Main pipeline (`run_annotation_pipeline.sh`)

```bash
cd HANNO_WithmRNA_HyloideaProteinSet_Pristimantis_euphronides_genome/Follow-up-scripts/
bash run_annotation_pipeline.sh \
  -i ../WithmRNA_HyloideaProteinSet_HANNO-RUN-Pristimantis_euphronides.genome.fasta \
  -o ../annotation_output \
  -g /data/GrenadaFrog144/Pristimantis_euphronides.genome.fasta \
  -p aPriEup1.0 --genome-size 1748533034
```

| Step | Script | Function |
|---|---|---|
| 1 | `beddb_to_gff3.py` | BED12→GFF3, removes coordinate-identical duplicates |
| 2 | `filter_hanno_fasta.py` | Filters AA/CDS/mRNA FASTAs to `retained_ids.txt` |
| 3 | `parse_hanno_annotation_summary.py` | Structural/functional statistics |
| 4 | — | Copies HANNO GTF |
| 5 | — | Verifies gene counts across files |

### 3.2 Gene counts by scaffold group (`count_genes_by_scaffold_group.sh`)

Single-pass `awk` classification of GFF3 gene features:
- `scaffold_1`–`scaffold_13` (13 longest, one per *E. coqui* hap1 chromosome)
- `scaffold_14`–`scaffold_31` (shorter RagTag scaffolds)
- `contig_*` (unplaced Flye contigs)

```bash
bash count_genes_by_scaffold_group.sh \
  -g ../annotation_output/aPriEup1.0_genomic.gff3 \
  -o ../annotation_output
```


## 4. NCBI Submission

See `NCBI_Zenodo_Submission/NCBI_Zenodo_Submission_README_github.md`.


## Author

Kopp K, Pristimantis euphronides genome project

## Citation

- Kuhl H (2024). HANNO. doi:10.5281/zenodo.11532370.
  https://github.com/HMPNK/HANNO
