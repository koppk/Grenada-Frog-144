#!/bin/bash
set -euo pipefail

# =============================================================================
# fetch_biosample_tissue.sh
#
# For each anuran genome assembly, fetch the BioSample tissue field and
# classify tissue source using the framework of Taberlet et al. (1999):
# destructive (animal killed), nondestructive (animal survives),
# ambiguous (muscle/tongue — biopsy or post-mortem unknown).
#
# Prerequisites:
#   conda install -c bioconda entrez-direct
#
# Usage:
#   bash fetch_biosample_tissue.sh -i anuran_genomes_ncbi.tsv -o output_dir
#
# Author: Kopp K, Pristimantis euphronides genome project
# =============================================================================

usage() {
    echo "Usage: $0 -i INPUT_TSV -o OUTPUT_DIR"
    exit 1
}

INPUT=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTPUT_DIR" ]]; then usage; fi

mkdir -p "$OUTPUT_DIR"

BIOSAMPLE_LIST="$OUTPUT_DIR/biosample_accessions.txt"
BIOSAMPLE_RAW="$OUTPUT_DIR/biosample_raw.jsonl"
BIOSAMPLE_TSV="$OUTPUT_DIR/biosample_tissue.tsv"
MERGED="$OUTPUT_DIR/assemblies_with_tissue.tsv"
REPORT="$OUTPUT_DIR/tissue_analysis.out"

{
echo "============================================================"
echo "BioSample tissue metadata fetch and classification"
echo "============================================================"
echo "Date:    $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Input:   $INPUT"
echo "============================================================"
echo ""

# -------------------------------------------------------------------
# 1. Extract unique BioSample accessions
# -------------------------------------------------------------------
echo ">>> Step 1: Extracting BioSample accessions..."
tail -n +2 "$INPUT" | cut -d$'\t' -f6 | sort -u | grep -v "^$" > "$BIOSAMPLE_LIST"
N_BS=$(wc -l < "$BIOSAMPLE_LIST")
echo "  Unique BioSample accessions: $N_BS"
echo ""

# -------------------------------------------------------------------
# 2. Batch fetch BioSample metadata
# -------------------------------------------------------------------
echo ">>> Step 2: Fetching BioSample metadata via efetch (batches of 50)..."

> "$BIOSAMPLE_RAW"
TOTAL=$(wc -l < "$BIOSAMPLE_LIST")

# Batch accessions into groups of 50 for efetch
BATCH_FILE="$OUTPUT_DIR/_batches.txt"
awk '{
    batch[NR] = $0
}
END {
    bsize = 50
    for (i = 1; i <= NR; i += bsize) {
        s = ""
        for (j = i; j < i + bsize && j <= NR; j++) {
            if (s != "") s = s ","
            s = s batch[j]
        }
        print s
    }
}' "$BIOSAMPLE_LIST" > "$BATCH_FILE"

BNUM=0
NBATCHES=$(wc -l < "$BATCH_FILE")
while IFS= read -r batch; do
    BNUM=$((BNUM + 1))
    efetch -db biosample -id "$batch" -format docsum 2>/dev/null \
        >> "$BIOSAMPLE_RAW" || true
    echo "  Batch $BNUM / $NBATCHES"
    sleep 0.4
done < "$BATCH_FILE"
rm -f "$BATCH_FILE"

echo ""

# -------------------------------------------------------------------
# 3. Extract tissue field from XML
# -------------------------------------------------------------------
echo ">>> Step 3: Extracting tissue fields..."

echo -e "biosample_accession\torganism\ttissue\tisolation_source\tsample_type" > "$BIOSAMPLE_TSV"

PYPARSER="$OUTPUT_DIR/_parse_biosample.py"
cat > "$PYPARSER" << 'PYEOF'
import xml.etree.ElementTree as ET
import sys

infile = sys.argv[1]

try:
    tree = ET.parse(infile)
except ET.ParseError:
    # Multiple XML documents concatenated — wrap in root
    with open(infile) as f:
        content = f.read()
    # Strip XML declarations and wrap
    import re
    content = re.sub(r'<\?xml[^?]*\?>\s*', '', content)
    content = re.sub(r'<!DOCTYPE[^>]*>\s*', '', content)
    content = "<root>" + content + "</root>"
    tree = ET.ElementTree(ET.fromstring(content))

seen = set()
for ds in tree.getroot().iter("DocumentSummary"):
    acc_el = ds.find("Accession")
    if acc_el is None:
        continue
    acc = acc_el.text or ""
    if acc in seen:
        continue
    seen.add(acc)

    org_el = ds.find("Organism")
    org = org_el.text if org_el is not None else ""

    tissue = ""
    isolation_source = ""
    sample_type = ""

    for attr in ds.iter("Attribute"):
        hname = (attr.get("harmonized_name") or "").lower()
        aname = (attr.get("attribute_name") or "").lower()
        val = (attr.text or "").strip()

        if hname == "tissue" or aname == "tissue_type" or aname == "tissue":
            if not tissue:
                tissue = val
        elif hname == "isolation_source" or aname == "isolation_source" or aname == "isolation source":
            isolation_source = val
        elif hname == "sample_type" or aname == "sample_type" or aname == "sample type":
            sample_type = val

    print(f"{acc}\t{org}\t{tissue}\t{isolation_source}\t{sample_type}")

PYEOF

python3 "$PYPARSER" "$BIOSAMPLE_RAW" >> "$BIOSAMPLE_TSV"
rm -f "$PYPARSER"

N_TISSUE=$(tail -n +2 "$BIOSAMPLE_TSV" | wc -l)
echo "  BioSamples with metadata: $N_TISSUE"
echo ""

# -------------------------------------------------------------------
# 4. Classify tissue source (Taberlet et al. 1999 framework)
# -------------------------------------------------------------------
echo ">>> Step 4: Classifying tissue source (Taberlet 1999)..."

CLASSIFIED="$OUTPUT_DIR/tissue_classified.tsv"
echo -e "biosample_accession\torganism\ttissue\ttissue_class" > "$CLASSIFIED"

tail -n +2 "$BIOSAMPLE_TSV" | awk -F'\t' '
BEGIN { OFS="\t" }
{
    acc = $1; org = $2; tissue = $3; iso = $4; stype = $5
    t = tolower(tissue " " iso " " stype)

    has_destructive = 0; has_ambiguous = 0; has_nondestructive = 0; has_unclassified = 0

    # Destructive: internal organs, whole organism, larvae
    if (t ~ /liver|brain|heart|kidney|gonad|ovary|testes|testis|spleen|lung|intestin|pancreas|eye|gut|whole|carcass|entire|euthan|larvae|tadpole/)
        has_destructive = 1

    # Ambiguous: muscle, tongue, bare skin, blood, erythrocytes
    # (could be from live animal or post-mortem collection)
    if (t ~ /muscle|tongue|blood|erythrocyte/) has_ambiguous = 1
    if (t ~ /skin/ && t !~ /skin swab/ && t !~ /skin clip/) has_ambiguous = 1

    # Nondestructive: swab, toe clip, fin clip, tail clip
    # (procedures that inherently require a live animal)
    if (t ~ /swab|buccal|oral swab|toe|toe clip|fin clip|tail clip|saliva|mucus|skin swab|skin clip/)
        has_nondestructive = 1

    # Most harmful present wins
    if (has_destructive) cls = "destructive"
    else if (has_ambiguous) cls = "ambiguous"
    else if (has_nondestructive) cls = "nondestructive"
    else if (t ~ /tissue|dna|genomic|mixed|misc|unknown/) cls = "unclassified"
    else if (tissue ~ /^[ \t]*$/ && iso ~ /^[ \t]*$/ && stype ~ /^[ \t]*$/) cls = "not_reported"
    else cls = "unclassified"

    print acc, org, tissue, cls
}' >> "$CLASSIFIED"

echo "  Classified: $CLASSIFIED"
echo ""

# -------------------------------------------------------------------
# 5. Merge with genome assembly data
# -------------------------------------------------------------------
echo ">>> Step 5: Merging with genome assembly data..."

echo -e "accession\torganism\ttech\tlevel\tsize\tsn50\tcn50\tbiosample\ttissue\ttissue_class" > "$MERGED"
tail -n +2 "$INPUT" | awk -F'\t' '
BEGIN { OFS="\t" }
NR==FNR && FNR==1 { next }
NR==FNR {
    tissue[$1] = $3
    cls[$1] = $4
    next
}
{
    bs = $6
    t = (bs in tissue) ? tissue[bs] : "not_fetched"
    c = (bs in cls) ? cls[bs] : "not_fetched"
    print $1, $2, $5, $4, $7, $8, $9, bs, t, c
}' "$CLASSIFIED" - >> "$MERGED"

echo "  Merged: $MERGED"
echo ""

# -------------------------------------------------------------------
# 6. Summary
# -------------------------------------------------------------------
echo "============================================================"
echo "SUMMARY"
echo "============================================================"
echo ""

echo "--- Tissue classification counts (Taberlet 1999) ---"
tail -n +2 "$MERGED" | cut -d$'\t' -f10 | sort | uniq -c | sort -rn
echo ""

echo "--- Nondestructive assemblies (animal survives) ---"
echo -e "accession\torganism\ttech\ttissue"
tail -n +2 "$MERGED" | awk -F'\t' '$10 == "nondestructive" {print $1"\t"$2"\t"$3"\t"$9}'
echo ""

echo "--- Ambiguous assemblies (muscle/tongue/skin — biopsy, swab or post-mortem) ---"
echo -e "accession\torganism\ttech\ttissue"
tail -n +2 "$MERGED" | awk -F'\t' '$10 == "ambiguous" {print $1"\t"$2"\t"$3"\t"$9}'
echo ""

echo "--- Unclassified / not reported ---"
echo -e "accession\torganism\ttissue_raw\tclass"
tail -n +2 "$MERGED" | awk -F'\t' '$10 == "unclassified" || $10 == "not_reported" {print $1"\t"$2"\t"$9"\t"$10}'
echo ""

echo "============================================================"

} > "$REPORT" 2>&1

echo "Done."
echo "Full report: $REPORT"
echo "Merged data: $MERGED"
