#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
#
# Run HapDup v0.12 via Docker to generate haplotype-resolved diploid
# assemblies from the Medaka-polished Flye assembly.
#
# Requires lr_mapping.bam (from run_minimap2_HapDup_prep.sh) and
# final_medaka_polished_assembly_consensus.fasta in the current directory.
#
# Output: hapdup/hapdup_dual_1.fasta, hapdup/hapdup_dual_2.fasta
#         hapdup/hapdup_phased_1.fasta, hapdup/hapdup_phased_2.fasta

HD_DIR=$(pwd)

docker run -v "${HD_DIR}:${HD_DIR}" -u "$(id -u):$(id -g)" \
    mkolmogo/hapdup:0.12 \
    hapdup \
    --assembly "${HD_DIR}/final_medaka_polished_assembly_consensus.fasta" \
    --bam "${HD_DIR}/lr_mapping.bam" \
    --out-dir "${HD_DIR}/hapdup" \
    -t 24 \
    --rtype ont
