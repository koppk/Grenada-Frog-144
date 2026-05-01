#!/bin/bash
# search_DM_domain.sh
# Author: Kopp K, Pristimantis euphronides genome project
#
# Search for DM-domain-containing genes in the P. euphronides assembly
# using tBLASTn with Xenopus laevis Dm-W and P. euphronides DMRT1 as queries.
#
# Usage: bash /data/GrenadaFrog144/SexChromosomes/DM_domain_search/scripts/search_DM_domain.sh
#
# Prerequisite: BLAST+ installed (makeblastdb, tblastn)

set -euo pipefail

BASEDIR="/data/GrenadaFrog144/SexChromosomes/DM_domain_search"
QUERYFILE="${BASEDIR}/dm_domain_queries.fasta"
ASSEMBLY="/data/GrenadaFrog144/Pristimantis_euphronides.genome.fasta"
OUTDIR="${BASEDIR}/output_DM_domain_search"
BLASTDB="${OUTDIR}/grenadafrog144_db"

mkdir -p "${OUTDIR}"

# ---------------------------------------------------------------
# Step 1: Verify query file exists
# ---------------------------------------------------------------
if [ ! -f "${QUERYFILE}" ]; then
    echo "ERROR: Query file not found: ${QUERYFILE}"
    exit 1
fi
echo "=== Step 1: Query file found: ${QUERYFILE} ==="

# ---------------------------------------------------------------
# Step 2: Build BLAST database from scaffolded assembly
# ---------------------------------------------------------------
if [ ! -f "${BLASTDB}.ndb" ]; then
    echo "=== Step 2: Building BLAST database ==="
    makeblastdb -in "${ASSEMBLY}" -dbtype nucl -parse_seqids \
        -out "${BLASTDB}"
else
    echo "=== Step 2: BLAST database already exists ==="
fi

# ---------------------------------------------------------------
# Step 3: Run tBLASTn
# ---------------------------------------------------------------
echo "=== Step 3: Running tBLASTn ==="

# Tabular output for parsing
tblastn \
    -query "${QUERYFILE}" \
    -db "${BLASTDB}" \
    -evalue 10 \
    -seg no \
    -word_size 2 \
    -max_target_seqs 100 \
    -num_threads 12 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
    -out "${OUTDIR}/dm_domain_tblastn_tab.txt"

# Full alignments for visual inspection
tblastn \
    -query "${QUERYFILE}" \
    -db "${BLASTDB}" \
    -evalue 10 \
    -seg no \
    -word_size 2 \
    -max_target_seqs 100 \
    -num_threads 12 \
    -outfmt 0 \
    -out "${OUTDIR}/dm_domain_tblastn_alignments.txt"

echo "=== Step 3: tBLASTn complete ==="

# ---------------------------------------------------------------
# Step 4: Parse results
# ---------------------------------------------------------------
echo ""
echo "================================================================"
echo "  ALL HITS"
echo "================================================================"
column -t "${OUTDIR}/dm_domain_tblastn_tab.txt" 2>/dev/null || \
    cat "${OUTDIR}/dm_domain_tblastn_tab.txt"

echo ""
echo "================================================================"
echo "  HITS ON SCAFFOLD_5 (known DMRT1/2/3 cluster — expected)"
echo "================================================================"
grep "scaffold_5" "${OUTDIR}/dm_domain_tblastn_tab.txt" || echo "(none)"

echo ""
echo "================================================================"
echo "  HITS OUTSIDE SCAFFOLD_5 (potential novel DM-domain genes)"
echo "================================================================"
grep -v "scaffold_5" "${OUTDIR}/dm_domain_tblastn_tab.txt" || echo "(none)"

echo ""
echo "================================================================"
echo "  HITS ON UNPLACED CONTIGS (potential W-linked DM-domain genes)"
echo "================================================================"
grep "contig_" "${OUTDIR}/dm_domain_tblastn_tab.txt" || echo "(none)"

echo ""
echo "================================================================"
echo "  INTERPRETATION"
echo "================================================================"
echo ""
echo "  Strong DM domain hit:   >60% identity, >40 aa aligned"
echo "  Plausible diverged hit: 30-60% identity, >30 aa aligned"
echo "  Likely false positive:  <30% identity or <25 aa aligned"
echo ""
echo "  Full alignments: ${OUTDIR}/dm_domain_tblastn_alignments.txt"
echo "  Tabular results: ${OUTDIR}/dm_domain_tblastn_tab.txt"
echo ""
echo "=== DONE ==="
