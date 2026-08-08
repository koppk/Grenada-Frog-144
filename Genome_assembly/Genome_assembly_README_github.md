# Genome assembly

De novo genome assemblies of *Pristimantis euphronides* from two
independent workflows.

Corresponds to Additional file 2, section **"Genome assembly"** and
Additional file 3, section **"Genome assembly statistics"**.

Input and output data files are deposited at Zenodo
([doi:10.5281/zenodo.15298546](https://doi.org/10.5281/zenodo.15298546)).

## Directory structure

On the server, the Flye output resided at
`/data/GrenadaFrog144/flye_nano-hq_HAC_output_GrenadaFrog144/`;
Medaka polishing was run locally on a GPU-equipped sequencing laptop.
For GitHub they are grouped under `Genome_assembly/`:

```
Genome_assembly/
├── Workflow_1/               Flye v2.9.5 + Medaka v2.0.1
└── Workflow_2/               Wtdbg v1.0 (H. Kuhl, CSA v2.6)
```

## Assembly commands

Run commands and parameters are documented in Additional file 2.

| Step | Tool | Key parameters |
|------|------|----------------|
| De novo assembly | Flye v2.9.5 | `--nano-hq`, 24 threads, 3 iterations |
| Consensus polishing | Medaka v2.0.1 | model `r1041_e82_400bps_hac_v4.3.0`, GPU-accelerated, batch size 32 |

## Prerequisites

- Flye v2.9.5
- Medaka v2.0.1 (GPU recommended)

## Author

Kopp K, Pristimantis euphronides genome project
