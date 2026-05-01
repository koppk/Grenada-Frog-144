#!/bin/bash
#
# download_and_combine_hyloidea_evidence.sh
# ==========================================
# Download protein, mRNA, and GTF annotation files for seven Hyloidea
# species from NCBI RefSeq FTP and concatenate them into combined
# evidence datasets for the HANNO v0.4 annotation pipeline.
#
# Replaces the original prefix_clean_concat_gtfs.sh, which handled
# only the GTF concatenation from pre-downloaded files.
#
# Species (7, all with BUSCO-evaluated RefSeq annotations):
#   Bufo gargarizans         GCF_014858855.1  (ASM1485885v1)
#   Bufo bufo                GCF_905171765.1  (aBufBuf1.1)
#   Dendropsophus ebraccatus GCF_027789765.1  (aDenEbr1.pat)
#   Eleutherodactylus coqui  GCF_035609145.1  (aEleCoq1.hap1)
#   Engystomops pustulosus   GCF_040894005.1  (aEngPut4.maternal)
#   Hyla sarda               GCF_029499605.1  (aHylSar1.hap1)
#   Ranitomeya imitator      GCF_032444005.1  (aRanImi1.pri)
#
# Per species, three files are downloaded from NCBI RefSeq FTP:
#   *_protein.faa.gz      Predicted protein sequences
#   *_rna.fna.gz          Predicted mRNA sequences
#   *_genomic.gtf.gz      Gene annotation in GTF format
#
# Output (3 combined files):
#   combined_Hyloidea_proteins.faa.gz   Concatenated protein sequences
#   combined_Hyloidea_rna.fna.gz        Concatenated mRNA sequences
#   combined_Hyloidea_species.gtf.gz    Concatenated GTF with species prefixes
#
# The GTF concatenation adds a six-letter species prefix (e.g. BufGar,
# EleCoq) to all gene_id and transcript_id fields to prevent identifier
# collisions across species.
#
# Requirements: wget, gzip, awk (gawk with gensub support)
#
# Usage:
#   mkdir -p Hyloidea_Proteins && cd Hyloidea_Proteins
#   bash download_and_combine_hyloidea_evidence.sh
#
# Author: Kopp K, Pristimantis euphronides genome project

set -euo pipefail

# ── Species definitions ────────────────────────────────────────────
# Format: accession|assembly_name|species_prefix
declare -A ASSEMBLIES=(
    ["GCF_014858855.1"]="ASM1485885v1"
    ["GCF_905171765.1"]="aBufBuf1.1"
    ["GCF_027789765.1"]="aDenEbr1.pat"
    ["GCF_035609145.1"]="aEleCoq1.hap1"
    ["GCF_040894005.1"]="aEngPut4.maternal"
    ["GCF_029499605.1"]="aHylSar1.hap1"
    ["GCF_032444005.1"]="aRanImi1.pri"
)

declare -A PREFIXES=(
    ["GCF_014858855.1"]="BufGar"
    ["GCF_905171765.1"]="BufBuf"
    ["GCF_027789765.1"]="DenEbr"
    ["GCF_035609145.1"]="EleCoq"
    ["GCF_040894005.1"]="EngPus"
    ["GCF_029499605.1"]="HylSar"
    ["GCF_032444005.1"]="RanImi"
)

declare -A SPECIES_NAMES=(
    ["GCF_014858855.1"]="Bufo gargarizans"
    ["GCF_905171765.1"]="Bufo bufo"
    ["GCF_027789765.1"]="Dendropsophus ebraccatus"
    ["GCF_035609145.1"]="Eleutherodactylus coqui"
    ["GCF_040894005.1"]="Engystomops pustulosus"
    ["GCF_029499605.1"]="Hyla sarda"
    ["GCF_032444005.1"]="Ranitomeya imitator"
)

# Ordered accession list (alphabetical by species)
ACCESSIONS=(
    GCF_014858855.1
    GCF_905171765.1
    GCF_027789765.1
    GCF_035609145.1
    GCF_040894005.1
    GCF_029499605.1
    GCF_032444005.1
)

# ── FTP URL construction ──────────────────────────────────────────
# GCF_014858855.1 -> GCF/014/858/855/GCF_014858855.1_ASM1485885v1
make_ftp_url() {
    local acc="$1"
    local asm="${ASSEMBLIES[$acc]}"
    local digits="${acc#GCF_}"
    digits="${digits%%.*}"
    local d1="${digits:0:3}"
    local d2="${digits:3:3}"
    local d3="${digits:6:3}"
    echo "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${d1}/${d2}/${d3}/${acc}_${asm}"
}

# ── Step 1: Download ──────────────────────────────────────────────
echo "=== Step 1: Downloading protein, mRNA, and GTF files ==="
echo "Date: $(date)"
echo ""

for acc in "${ACCESSIONS[@]}"; do
    asm="${ASSEMBLIES[$acc]}"
    name="${SPECIES_NAMES[$acc]}"
    base_url=$(make_ftp_url "$acc")
    base_name="${acc}_${asm}"

    echo "  ${name} (${acc})"

    for suffix in protein.faa.gz rna.fna.gz genomic.gtf.gz; do
        url="${base_url}/${base_name}_${suffix}"
        outfile="${base_name}_${suffix}"
        if [ -f "$outfile" ]; then
            echo "    ${suffix}: exists, skipping"
        else
            echo "    ${suffix}: downloading"
            wget -q -O "$outfile" "$url"
        fi
    done
    echo ""
done

echo "  Download complete."
echo ""

# ── Step 2: Concatenate protein sequences ─────────────────────────
echo "=== Step 2: Concatenating protein sequences ==="

PROT_OUT="combined_Hyloidea_proteins.faa.gz"
PROT_TMP="combined_Hyloidea_proteins.faa"
> "$PROT_TMP"

for acc in "${ACCESSIONS[@]}"; do
    asm="${ASSEMBLIES[$acc]}"
    name="${SPECIES_NAMES[$acc]}"
    infile="${acc}_${asm}_protein.faa.gz"
    n_seq=$(zcat "$infile" | grep -c '^>')
    echo "  ${name}: ${n_seq} sequences"
    zcat "$infile" >> "$PROT_TMP"
done

n_total=$(grep -c '^>' "$PROT_TMP")
echo "  Total: ${n_total} sequences"
gzip -c "$PROT_TMP" > "$PROT_OUT"
rm "$PROT_TMP"
echo "  Output: ${PROT_OUT} ($(du -h "$PROT_OUT" | cut -f1))"
echo ""

# ── Step 3: Concatenate mRNA sequences ────────────────────────────
echo "=== Step 3: Concatenating mRNA sequences ==="

RNA_OUT="combined_Hyloidea_rna.fna.gz"

cat *_rna.fna.gz > "$RNA_OUT"

echo "  Output: ${RNA_OUT} ($(du -h "$RNA_OUT" | cut -f1))"
echo ""

# ── Step 4: Prefix and concatenate GTF files ──────────────────────
echo "=== Step 4: Concatenating GTF files with species prefixes ==="

GTF_OUT="combined_Hyloidea_species.gtf.gz"
GTF_TMP="combined_Hyloidea_species.gtf"
> "$GTF_TMP"

for acc in "${ACCESSIONS[@]}"; do
    asm="${ASSEMBLIES[$acc]}"
    prefix="${PREFIXES[$acc]}"
    name="${SPECIES_NAMES[$acc]}"
    infile="${acc}_${asm}_genomic.gtf.gz"

    echo "  ${name}: prefix [${prefix}]"

    zcat "$infile" | awk -v pfx="$prefix" '
    BEGIN { OFS="\t" }
    /^#/ { next }
    {
        gsub(/\r/, "")
        $9 = gensub(/gene_id "([^"]+)"/, "gene_id \"" pfx "_\\1\"", "g", $9)
        $9 = gensub(/transcript_id "([^"]+)"/, "transcript_id \"" pfx "_\\1\"", "g", $9)
        print
    }' >> "$GTF_TMP"
done

gzip -c "$GTF_TMP" > "$GTF_OUT"
rm "$GTF_TMP"
echo "  Output: ${GTF_OUT} ($(du -h "$GTF_OUT" | cut -f1))"
echo ""

# ── Summary ────────────────────────────────────────────────────────
echo "=== Output files ==="
echo "  ${PROT_OUT}  ($(du -h "$PROT_OUT" | cut -f1))"
echo "  ${RNA_OUT}  ($(du -h "$RNA_OUT" | cut -f1))"
echo "  ${GTF_OUT}  ($(du -h "$GTF_OUT" | cut -f1))"
echo ""
echo "Per-species downloads retained in working directory."
echo "Done: $(date)"
