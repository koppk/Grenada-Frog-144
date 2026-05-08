# *Pristimantis euphronides* genome project — scripts repository

<p align="center">
  <img src="GrenadaFrogLogo.png" alt="Grenada Frog Logo" width="300">
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15298547.svg)](https://doi.org/10.5281/zenodo.15298547)

This repository contains the computational scripts used in the genome assembly, annotation, and analysis of *Pristimantis euphronides* (Grenada frog, specimen GrenadaFrog144), accompanying the manuscript:

> Kopp K, Lownds S, Harrison B, Kuhl H, Moittie S. From Field to Genome: Non-Lethal, On-Site Nanopore-Only Assembly of the Endangered Grenada Frog (*Pristimantis euphronides*). *BMC Genomics*. In submission.

## Abstract

The Grenada Frog (*Pristimantis euphronides*) is a critically endangered direct-developing anuran endemic to montane rainforest on Grenada in the Lesser Antilles. Reference genomes are foundational tools for the conservation of endangered species, but their production typically depends on destructive sampling, large institutional infrastructure, or the shipment of biomaterials abroad.

We present a 1.75 Gb scaffold-level reference genome of *P. euphronides* assembled from a single wild female sampled non-lethally by toe clip, blood swab, and buccal swab. All wet-laboratory steps, including Oxford Nanopore Technologies (ONT) MinION sequencing, were performed in-country on portable equipment, and downstream analyses ran on laptop-class hardware and a remote server accessed from Grenada. The Flye-Medaka assembly, scaffolded with RagTag against *Eleutherodactylus coqui*, comprised 13 pseudochromosomes with a scaffold N50 of 117.8 Mb, achieved 91.4% tetrapod BUSCO completeness, and supported 28,071 annotated protein-coding genes. Cross-species synteny analysis against three additional Hyloidea genomes computationally predicted an underlying karyotype of n = 16, consistent with prior cytogenetic counts. Coverage analysis identified candidate Z-linked regions on two scaffolds and a set of unplaced contigs with sequence features consistent with the heterochromatic W chromosome described from C-banding, providing a first-pass framework for sex chromosome analyses in *Pristimantis*.

This whole-genome assembly provides a reference for phylogenomic and conservation genomic studies in a hyperdiverse direct-developing frog clade in which previous comparative work had relied on few marker genes. In the absence of chromatin conformation capture methods such as Hi-C, which were not compatible with non-lethal field sampling, reference-guided scaffolding against the closest available chromosome-scale Hyloidea genome, approximately 52 million years divergent, provided a workable alternative. The non-lethal, in-country, ONT-only workflow demonstrates that reference-quality genomes for repeat-rich vertebrates of conservation concern can be produced without destructive sampling or outsourcing to external sequencing or bioinformatics services.

## Scope

This repository contains **scripts and commands only**. All input data, output files, intermediate results, and detailed per-directory documentation (including tool versions, parameters, and reproducibility notes) are deposited on Zenodo (DOI: [10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)). The Zenodo deposit mirrors the directory structure of this repository.

## Repository structure

| Directory | Description |
|-----------|-------------|
| `Reads/` | Sequencing read information |
| `Read_quality_analysis/` | Read quality metrics and Phred score computation |
| `Genome_size_and_complexity_estimation/` | Genome size estimation |
| `Genome_assembly/` | Genome assembly with Flye and polishing with Medaka |
| `Assembly_quality_assessment/` | Inspector, Merqury, BUSCO, and inter-workflow concordance |
| `Taxonomic_classification_of_reads/` | Kraken2 contamination screening |
| `Reference_genome_selection_via_phylogenetic_proximity/` | Phylogenetic reference selection |
| `Reference_genome_selection_via_k-mer-based_distance_estimation/` | Mash distance estimation and heatmap |
| `Reference-guided_scaffolding_of_genome_assemblies/` | RagTag scaffolding against *Eleutherodactylus coqui* |
| `Evaluation_of_reference-guided_scaffolding/` | Synteny analysis and karyotype schematic |
| `Read_mapping_and_genome_coverage_analysis/` | Read mapping, coverage computation, and anomaly detection |
| `Annotation_scaffolded_genome_assembly/` | HANNO gene annotation and NCBI/Zenodo reformatting |
| `Identification_of_candidate_sex-chromosome-linked_regions/` | Z and W candidate regions, gametologs, banding analysis |
| `Phylogenetic_analysis_using_nuclear_gene_orthologs/` | Supermatrix phylogenetics with IQ-TREE |
| `Alignment_with_previously_published_Pristimantis_euphronides_and_Eleutherodactylus_johnstonei_sequences/` | Alignment with published sequences |
| `Haplotype_analysis/` | HapDup haplotype phasing |
| `Survey_of_anuran_genome_assemblies_in_NCBI_Genome/` | Survey of publicly available anuran genomes |
| `Divergence_times_and_biogeographic_context/` | Divergence time estimation and biogeographic analysis |

Each directory contains a chapter-specific README and a `scripts/` subfolder with the analysis scripts. Scripts are written in bash and Python 3.

## Data availability

- **Raw sequencing reads**: ENA/NCBI BioProject [PRJEB89028](https://www.ebi.ac.uk/ena/browser/view/PRJEB89028)
- **Assembly**: ENA [GCA_965278355.2](https://www.ebi.ac.uk/ena/browser/view/GCA_965278355.2), NCBI [GCA_965278355.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_965278355.2/) (aPriEup1.0)
- **Complete analysis deposit (data, results, documentation)**: Zenodo [10.5281/zenodo.15298547](https://doi.org/10.5281/zenodo.15298547)

## License

This repository is licensed under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.
