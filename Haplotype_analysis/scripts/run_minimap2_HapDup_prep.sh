#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Map ONT HAC reads to the Medaka-polished Flye assembly for HapDup input.
# Produces lr_mapping.bam (sorted, indexed) as required by HapDup v0.12.
#
# Input:  final_medaka_polished_assembly_consensus.fasta (in current directory)
#         GrenadaFrog144_ONT_HAC_all.fastq.gz (one directory up)
# Output: lr_mapping.bam, lr_mapping.bam.bai

minimap2 -ax map-ont -t 24 \
    final_medaka_polished_assembly_consensus.fasta \
    ../GrenadaFrog144_ONT_HAC_all.fastq.gz \
    | samtools sort -@ 4 -m 4G > lr_mapping.bam

samtools index -@ 4 lr_mapping.bam
