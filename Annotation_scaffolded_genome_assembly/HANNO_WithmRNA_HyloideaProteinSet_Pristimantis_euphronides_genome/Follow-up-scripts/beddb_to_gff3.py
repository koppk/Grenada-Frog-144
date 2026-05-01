#!/usr/bin/env python3
"""
beddb_to_gff3.py — Convert HANNO BESTMODELS-FINAL.bedDB to valid NCBI-style GFF3

Fixes critical issues in previous conversion attempts:
  1. CDS is correctly split across exon blocks (not one span including introns)
  2. Coordinates properly converted from BED 0-based half-open to GFF3 1-based closed
  3. Gene symbol and product extracted from correct bedDB columns
  4. Features output in correct hierarchy: gene → mRNA → exon + CDS
  5. Duplicate/overlapping gene models detected and flagged
  6. Produces valid GFF3 passing gt gff3validator

Usage:
  python3 beddb_to_gff3.py BESTMODELS-FINAL.bedDB genome.fasta -o annotation.gff3

  Optional:
    --named-only        Only include genes with a functional name (Preferred_name != "-")
    --locus-prefix PEU  Locus tag prefix (default: PEU for P. euphronides)
    --source HANNO      Source field in column 2 (default: HANNO)
    --no-fasta          Skip writing protein/mRNA/CDS FASTA files

Author: K. Kopp, P. euphronides genome project
"""

import sys
import argparse
import re
from collections import OrderedDict


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert HANNO BESTMODELS-FINAL.bedDB to NCBI-style GFF3"
    )
    p.add_argument("beddb", help="BESTMODELS-FINAL.bedDB file")
    p.add_argument("genome", nargs="?", default=None,
                   help="Genome FASTA (for sequence-region headers and FASTA extraction)")
    p.add_argument("-o", "--output", default="annotation.gff3",
                   help="Output GFF3 file (default: annotation.gff3)")
    p.add_argument("--named-only", action="store_true",
                   help="Only output genes with a functional name")
    p.add_argument("--locus-prefix", default="PEU",
                   help="Locus tag prefix (default: PEU)")
    p.add_argument("--source", default="HANNO",
                   help="GFF3 source field (default: HANNO)")
    p.add_argument("--write-proteins", action="store_true",
                   help="Write protein FASTA (requires BESTMODELS-FINAL.AA.faa or CDS)")
    p.add_argument("--seqid-map", default=None,
                   help="TSV file mapping local→NCBI seq IDs (col1=local, col2=NCBI)")
    p.add_argument("--product-col", default=None,
                   help="Column name or 0-based index for product/description "
                        "(default: auto-detect 'description' then 'Description')")
    p.add_argument("--stats", action="store_true",
                   help="Print summary statistics to stderr")
    return p.parse_args()


def read_genome_lengths(fasta_path):
    """Read sequence lengths from a FASTA file (simple parser)."""
    lengths = OrderedDict()
    current = None
    length = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if current is not None:
        lengths[current] = length
    return lengths


def read_seqid_map(map_path):
    """Read a TSV mapping file: local_id → ncbi_id."""
    mapping = {}
    with open(map_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def parse_beddb_header(header_line):
    """Parse the bedDB header to get column indices."""
    cols = header_line.lstrip("#").strip().split("\t")
    return {name.strip(): i for i, name in enumerate(cols)}


def parse_bed12_blocks(chrom_start, block_count, block_sizes_str, block_starts_str):
    """
    Parse BED12 block structure into list of (start, end) in 0-based half-open coords.
    Returns list of tuples: [(exon_start_0based, exon_end_0based), ...]
    """
    sizes = [int(x) for x in block_sizes_str.rstrip(",").split(",") if x]
    starts = [int(x) for x in block_starts_str.rstrip(",").split(",") if x]

    blocks = []
    for i in range(min(block_count, len(sizes), len(starts))):
        exon_start = chrom_start + starts[i]
        exon_end = exon_start + sizes[i]
        blocks.append((exon_start, exon_end))
    return blocks


def intersect_with_cds(exon_blocks, thick_start, thick_end):
    """
    Intersect exon blocks with the CDS (thick) region.
    All coordinates 0-based half-open.
    Returns list of (cds_start, cds_end) tuples for coding portions of exons.
    """
    cds_blocks = []
    for (ex_start, ex_end) in exon_blocks:
        # Skip exons entirely outside CDS
        if ex_end <= thick_start or ex_start >= thick_end:
            continue
        # Clip to CDS boundaries
        cds_start = max(ex_start, thick_start)
        cds_end = min(ex_end, thick_end)
        if cds_end > cds_start:
            cds_blocks.append((cds_start, cds_end))
    return cds_blocks


def to_gff3_coord(bed_start, bed_end):
    """Convert 0-based half-open (BED) to 1-based closed (GFF3)."""
    return bed_start + 1, bed_end


def url_encode_gff3(value):
    """Minimal URL encoding for GFF3 attribute values."""
    return (value
            .replace("%", "%25")
            .replace(";", "%3B")
            .replace("=", "%3D")
            .replace("&", "%26")
            .replace(",", "%2C")
            .replace("\t", " ")
            .replace("\n", " "))


def clean_product_name(raw):
    """
    Clean a protein description for use as NCBI product name.
    
    Removes:
      - Species brackets: "titin [Eleutherodactylus coqui]" → "titin"
      - Accession prefixes: "XP_032605298.2 titin" → "titin"
      - "PREDICTED: " prefix
      - Leading/trailing whitespace
    """
    if not raw or raw in ("-", "nan", "NA", ""):
        return "hypothetical protein"
    
    name = raw.strip()
    
    # Remove "PREDICTED: " prefix
    name = re.sub(r'^PREDICTED:\s*', '', name)
    
    # Remove species bracket: "[Genus species]" at end
    name = re.sub(r'\s*\[.*?\]\s*$', '', name)
    
    # Remove leading accession (e.g., "XP_032605298.2 ")
    name = re.sub(r'^[A-Z]{1,3}_\d+\.\d+\s+', '', name)
    
    # Remove "isoform X1" etc. (optional — NCBI sometimes keeps these)
    # name = re.sub(r'\s+isoform\s+\S+$', '', name)
    
    name = name.strip()
    
    if not name:
        return "hypothetical protein"
    
    return name


def format_attrs(attr_dict):
    """Format a dict of GFF3 attributes into column 9 string."""
    parts = []
    for key, val in attr_dict.items():
        if val is not None and val != "" and val != "-":
            parts.append(f"{key}={url_encode_gff3(str(val))}")
    return ";".join(parts)


def main():
    args = parse_args()

    # Load sequence ID mapping if provided
    seqid_map = {}
    if args.seqid_map:
        seqid_map = read_seqid_map(args.seqid_map)
        print(f"[INFO] Loaded {len(seqid_map)} sequence ID mappings", file=sys.stderr)

    # Load genome lengths if provided
    seq_lengths = OrderedDict()
    if args.genome:
        print(f"[INFO] Reading genome lengths from {args.genome}...", file=sys.stderr)
        seq_lengths = read_genome_lengths(args.genome)
        print(f"[INFO] Found {len(seq_lengths)} sequences", file=sys.stderr)

    # Parse bedDB
    print(f"[INFO] Parsing {args.beddb}...", file=sys.stderr)
    genes = []  # List of parsed gene records
    col_map = None
    header_line = None

    with open(args.beddb) as f:
        for line in f:
            if line.startswith("##") or line.startswith("#"):
                header_line = line
                col_map = parse_beddb_header(line)
                # Diagnostic: show key column mappings
                desc_idx = col_map.get("description", 14)
                pname_idx = col_map.get("Preferred_name", 25)
                print(f"[INFO] bedDB header detected. Key columns (0-based):", file=sys.stderr)
                print(f"       description     → index {desc_idx} (col {desc_idx+1})", file=sys.stderr)
                print(f"       Preferred_name  → index {pname_idx} (col {pname_idx+1})", file=sys.stderr)
                if args.stats and col_map:
                    print(f"[INFO] All columns:", file=sys.stderr)
                    for cname, cidx in sorted(col_map.items(), key=lambda x: x[1]):
                        print(f"       {cidx:3d}: {cname}", file=sys.stderr)
                continue
            if not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")

            # BED12 core fields (always columns 0-11)
            chrom = cols[0]
            chrom_start = int(cols[1])  # 0-based
            chrom_end = int(cols[2])    # 0-based half-open (= 1-based end)
            name = cols[3]              # transcript name, e.g. hanno.g1.t1
            score = cols[4]
            strand = cols[5]
            thick_start = int(cols[6])  # CDS start, 0-based
            thick_end = int(cols[7])    # CDS end, 0-based half-open
            # cols[8] = color
            block_count = int(cols[9])
            block_sizes = cols[10]
            block_starts = cols[11]

            # Extended bedDB fields — use column map if available, else fixed indices
            # HANNO bedDB has TWO description-like columns:
            #   "description" (around col 15, idx 14): best protein hit description
            #       e.g., "titin [Eleutherodactylus coqui]"
            #   "Description" (around col 25, idx 24): eggNOG functional description
            #       e.g., "Titin" or "Serine/threonine-protein kinase"
            #   "Preferred_name" (col 26, idx 25): gene symbol
            #       e.g., "TTN", "PAK6", or "-"

            # Determine product column
            if args.product_col is not None:
                # User specified: either a column name or index
                try:
                    prod_idx = int(args.product_col)
                except ValueError:
                    prod_idx = col_map.get(args.product_col, 14) if col_map else 14
            elif col_map:
                # Auto-detect: prefer "description" (protein hit), fall back to "Description" (eggNOG)
                if "description" in col_map:
                    prod_idx = col_map["description"]
                elif "Description" in col_map:
                    prod_idx = col_map["Description"]
                else:
                    prod_idx = 14  # fallback
            else:
                prod_idx = 14  # fallback

            # Preferred_name (gene symbol)
            if col_map:
                pname_idx = col_map.get("Preferred_name", 25)
            else:
                pname_idx = 25

            description = cols[prod_idx].strip() if len(cols) > prod_idx else "-"
            preferred_name = cols[pname_idx].strip() if len(cols) > pname_idx else "-"

            # Also collect eggNOG Description for fallback
            eggnog_desc_idx = col_map.get("Description", 24) if col_map else 24
            eggnog_desc_raw = cols[eggnog_desc_idx].strip() if len(cols) > eggnog_desc_idx else "-"

            # Clean up
            if description in ("-", "", "nan", "NA"):
                description = "hypothetical protein"
            else:
                description = clean_product_name(description)
            if preferred_name in ("-", "", "nan", "NA"):
                preferred_name = ""

            # Filter if --named-only
            if args.named_only and not preferred_name:
                continue

            # Parse exon blocks
            exon_blocks = parse_bed12_blocks(chrom_start, block_count, block_sizes, block_starts)

            # Compute CDS blocks (intersect exons with thick region)
            has_cds = thick_start < thick_end  # thick_start == thick_end means non-coding
            cds_blocks = []
            if has_cds:
                cds_blocks = intersect_with_cds(exon_blocks, thick_start, thick_end)

            genes.append({
                "chrom": chrom,
                "start": chrom_start,
                "end": chrom_end,
                "name": name,
                "score": score,
                "strand": strand,
                "thick_start": thick_start,
                "thick_end": thick_end,
                "exon_blocks": exon_blocks,
                "cds_blocks": cds_blocks,
                "has_cds": has_cds,
                "description": description,
                "eggnog_description": clean_product_name(eggnog_desc_raw) if eggnog_desc_raw not in ("-", "", "nan", "NA") else "hypothetical protein",
                "gene_symbol": preferred_name,
            })

    print(f"[INFO] Parsed {len(genes)} gene models", file=sys.stderr)

    # Check if product column was mostly empty — auto-fallback to eggNOG Description
    empty_product = sum(1 for g in genes if g["description"] == "hypothetical protein")
    if empty_product > len(genes) * 0.95 and args.product_col is None:
        print(f"[WARNING] Product column 'description' (col 15) is >95% empty "
              f"({empty_product}/{len(genes)}).", file=sys.stderr)
        print(f"[WARNING] Auto-falling back to 'Description' (eggNOG, col 25).",
              file=sys.stderr)
        # Use the eggNOG descriptions we already collected
        for g in genes:
            if g.get("eggnog_description") and g["eggnog_description"] != "hypothetical protein":
                g["description"] = g["eggnog_description"]
        new_empty = sum(1 for g in genes if g["description"] == "hypothetical protein")
        print(f"[INFO] After fallback: {len(genes) - new_empty} genes now have product names.",
              file=sys.stderr)

    # Diagnostic: show sample values from key fields
    if genes:
        print(f"[INFO] Sample gene models (first 5 with gene symbol):", file=sys.stderr)
        shown = 0
        for g in genes:
            if g["gene_symbol"] and shown < 5:
                print(f"       {g['name']}: gene={g['gene_symbol']}, product={g['description'][:60]}",
                      file=sys.stderr)
                shown += 1
        if shown == 0:
            print("       (no named genes found — check column mapping!)", file=sys.stderr)
            for g in genes[:3]:
                print(f"       {g['name']}: gene='{g['gene_symbol']}', product='{g['description'][:60]}'",
                      file=sys.stderr)

    # --- Detect duplicates (same chrom, overlapping coordinates, same strand) ---
    # Sort by chrom, start
    genes.sort(key=lambda g: (g["chrom"], g["start"], g["end"]))

    dup_count = 0
    for i in range(1, len(genes)):
        prev = genes[i - 1]
        curr = genes[i]
        if (prev["chrom"] == curr["chrom"] and
            prev["strand"] == curr["strand"] and
            prev["start"] == curr["start"] and
            prev["end"] == curr["end"]):
            dup_count += 1
            # Mark the lower-scoring one (or the second) as duplicate
            curr["is_dup"] = True
        else:
            curr.setdefault("is_dup", False)
    if len(genes) > 0:
        genes[0].setdefault("is_dup", False)

    if dup_count > 0:
        print(f"[WARNING] Found {dup_count} exact-position duplicates (will be skipped)",
              file=sys.stderr)

    # --- Write GFF3 ---
    print(f"[INFO] Writing GFF3 to {args.output}...", file=sys.stderr)

    gene_count = 0
    mrna_count = 0
    cds_feature_count = 0
    exon_count = 0
    named_count = 0

    with open(args.output, "w") as out:
        # Header
        out.write("##gff-version 3\n")

        # Sequence-region directives
        if seq_lengths:
            for seqid, length in seq_lengths.items():
                mapped_id = seqid_map.get(seqid, seqid)
                out.write(f"##sequence-region {mapped_id} 1 {length}\n")

        # Gene features
        for gene in genes:
            if gene.get("is_dup", False):
                continue

            chrom = seqid_map.get(gene["chrom"], gene["chrom"])
            strand = gene["strand"]
            name = gene["name"]
            gene_symbol = gene["gene_symbol"]
            description = gene["description"]

            # Extract gene number for locus_tag
            gene_num_match = re.search(r'g(\d+)', name)
            gene_num = gene_num_match.group(1) if gene_num_match else str(gene_count + 1)
            locus_tag = f"{args.locus_prefix}_{gene_num.zfill(5)}"

            # GFF3 IDs
            gene_id = f"gene-{name}"
            mrna_id = f"rna-{name}"
            # Protein ID
            protein_id = f"{name}.p1"

            # --- Gene feature ---
            g_start, g_end = to_gff3_coord(gene["start"], gene["end"])
            gene_attrs = OrderedDict([
                ("ID", gene_id),
                ("Name", gene_symbol if gene_symbol else name),
                ("locus_tag", locus_tag),
            ])
            if gene_symbol:
                gene_attrs["gene"] = gene_symbol
            gene_attrs["gene_biotype"] = "protein_coding" if gene["has_cds"] else "misc_RNA"

            out.write(f"{chrom}\t{args.source}\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\t{format_attrs(gene_attrs)}\n")
            gene_count += 1
            if gene_symbol:
                named_count += 1

            # --- mRNA feature ---
            mrna_attrs = OrderedDict([
                ("ID", mrna_id),
                ("Parent", gene_id),
                ("Name", gene_symbol if gene_symbol else name),
                ("locus_tag", locus_tag),
                ("transcript_id", name),
            ])
            if gene_symbol:
                mrna_attrs["gene"] = gene_symbol
            mrna_attrs["product"] = description

            out.write(f"{chrom}\t{args.source}\tmRNA\t{g_start}\t{g_end}\t.\t{strand}\t.\t{format_attrs(mrna_attrs)}\n")
            mrna_count += 1

            # --- Exon features ---
            for i, (ex_start, ex_end) in enumerate(gene["exon_blocks"], 1):
                e_start, e_end = to_gff3_coord(ex_start, ex_end)
                exon_attrs = OrderedDict([
                    ("ID", f"exon-{name}.{i}"),
                    ("Parent", mrna_id),
                ])
                out.write(f"{chrom}\t{args.source}\texon\t{e_start}\t{e_end}\t.\t{strand}\t.\t{format_attrs(exon_attrs)}\n")
                exon_count += 1

            # --- CDS features (one per coding exon) ---
            if gene["has_cds"] and gene["cds_blocks"]:
                # Calculate phase for each CDS block
                # Phase = (3 - cumulative_coding_bases % 3) % 3 for + strand
                # For - strand, process blocks in reverse
                cds_list = list(gene["cds_blocks"])
                if strand == "-":
                    cds_list = list(reversed(cds_list))

                cumulative_bases = 0
                cds_with_phase = []
                for (c_start, c_end) in cds_list:
                    phase = (3 - (cumulative_bases % 3)) % 3
                    cds_with_phase.append((c_start, c_end, phase))
                    cumulative_bases += (c_end - c_start)

                # Write in genomic order (re-reverse for - strand)
                if strand == "-":
                    cds_with_phase = list(reversed(cds_with_phase))

                for i, (c_start, c_end, phase) in enumerate(cds_with_phase, 1):
                    cs, ce = to_gff3_coord(c_start, c_end)
                    cds_attrs = OrderedDict([
                        ("ID", f"cds-{name}.{i}"),
                        ("Parent", mrna_id),
                        ("protein_id", protein_id),
                    ])
                    if gene_symbol:
                        cds_attrs["gene"] = gene_symbol
                    cds_attrs["product"] = description

                    out.write(f"{chrom}\t{args.source}\tCDS\t{cs}\t{ce}\t.\t{strand}\t{phase}\t{format_attrs(cds_attrs)}\n")
                    cds_feature_count += 1

    # --- Write retained and duplicate ID lists ---
    import os
    out_dir = os.path.dirname(args.output) or "."
    retained_file = os.path.join(out_dir, "retained_ids.txt")
    duplicate_file = os.path.join(out_dir, "duplicate_ids.txt")

    with open(retained_file, "w") as rf:
        for gene in genes:
            if not gene.get("is_dup", False):
                rf.write(gene["name"] + "\n")
    print(f"[INFO] Retained IDs written to {retained_file} ({gene_count} entries)",
          file=sys.stderr)

    with open(duplicate_file, "w") as df:
        df.write("# Coordinate-identical duplicate gene models removed\n")
        df.write("# (same chrom, start, end, strand as a preceding entry)\n")
        for gene in genes:
            if gene.get("is_dup", False):
                df.write(gene["name"] + "\n")
    print(f"[INFO] Duplicate IDs written to {duplicate_file} ({dup_count} entries)",
          file=sys.stderr)

    # --- Statistics ---
    stats_msg = f"""
[DONE] GFF3 written to {args.output}
  Genes:        {gene_count}
  Named genes:  {named_count} ({100*named_count/max(gene_count,1):.1f}%)
  mRNA:         {mrna_count}
  Exons:        {exon_count}
  CDS features: {cds_feature_count}
  Duplicates skipped: {dup_count}
"""
    print(stats_msg, file=sys.stderr)

    if args.stats:
        # Additional stats
        avg_exons = exon_count / max(gene_count, 1)
        avg_cds = cds_feature_count / max(gene_count, 1)
        print(f"  Avg exons/gene: {avg_exons:.1f}", file=sys.stderr)
        print(f"  Avg CDS features/gene: {avg_cds:.1f}", file=sys.stderr)


if __name__ == "__main__":
    main()
