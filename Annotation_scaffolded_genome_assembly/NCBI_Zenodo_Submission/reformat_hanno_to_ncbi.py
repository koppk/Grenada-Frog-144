#!/usr/bin/env python3
"""
reformat_hanno_to_ncbi.py — Reformat HANNO-based annotation files into
NCBI-compliant format for two distinct purposes:

  OUTPUT MODE 1 — "zenodo" (default):
    NCBI FTP-style files for Zenodo/reviewers.
    - Scaffold names → NCBI INSDC accessions (CBDIFN02######.1)
    - hanno.g#.t# → sequential locus_tag numbering (PRIEUP_00001)
    - Clean FASTA headers matching NCBI FTP conventions
    - NO gnl|dbname| in headers (that's submission-internal only)

  OUTPUT MODE 2 — "table2asn":
    GFF3 ready for table2asn processing and direct NCBI submission.
    - Keeps ORIGINAL scaffold seqids (must match genome FASTA)
    - protein_id = gnl|dbname|locus_tag   (NCBI requirement, gets replaced
      with real accessions like CAN####### during processing)
    - transcript_id = gnl|dbname|mrna.locus_tag

Both modes always produce:
  - hanno_to_ncbi_id_mapping.tsv (traceability: old → new IDs)

Usage:
  # Zenodo files (for reviewers — looks like NCBI FTP)
  python3 reformat_hanno_to_ncbi.py \\
      --gff3 aPriEup1.0_genomic.gff3 \\
      --cds  aPriEup1.0_cds_from_genomic.fna \\
      --protein aPriEup1.0_protein.faa \\
      --mrna aPriEup1.0_rna_from_genomic.fna \\
      --seqid-map scaffold_to_accession.tsv \\
      --locus-prefix PRIEUP \\
      -o zenodo_annotation/

  # table2asn GFF3 (for NCBI submission — with gnl| protein_ids)
  python3 reformat_hanno_to_ncbi.py \\
      --gff3 aPriEup1.0_genomic.gff3 \\
      --seqid-map scaffold_to_accession.tsv \\
      --locus-prefix PRIEUP \\
      --mode table2asn \\
      --dbname KoppSGU \\
      -o ncbi_submission/

Author: K. Kopp, P. euphronides genome project
"""

import argparse
import os
import sys
import re
from collections import OrderedDict
from urllib.parse import unquote


# =============================================================================
# Argument parsing
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Reformat HANNO annotation → NCBI-compliant files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--gff3", required=True, help="Input GFF3 from beddb_to_gff3.py")
    p.add_argument("--cds", default=None, help="Input CDS FASTA")
    p.add_argument("--protein", default=None, help="Input protein FASTA")
    p.add_argument("--mrna", default=None, help="Input mRNA FASTA")
    p.add_argument("--seqid-map", required=True,
                    help="TSV: scaffold_name → NCBI accession")
    p.add_argument("--locus-prefix", default="PRIEUP",
                    help="Registered locus_tag prefix (default: PRIEUP)")
    p.add_argument("--mode", choices=["zenodo", "table2asn", "both"], default="both",
                    help="Output mode (default: both)")
    p.add_argument("--dbname", default="KoppSGU",
                    help="Lab identifier for gnl|dbname| in table2asn mode (default: KoppSGU). "
                         "Only used in table2asn GFF3; never appears in Zenodo files.")
    p.add_argument("--source", default="HANNO",
                    help="GFF3 source field, column 2 (default: HANNO)")
    p.add_argument("--assembly-prefix", default="aPriEup1.0",
                    help="Assembly name for output filenames (default: aPriEup1.0)")
    p.add_argument("-o", "--outdir", default="ncbi_annotation",
                    help="Output directory (default: ncbi_annotation)")
    return p.parse_args()


# =============================================================================
# Helpers
# =============================================================================
def load_seqid_map(path):
    """Load scaffold→NCBI accession mapping (TSV)."""
    mapping = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def parse_gff3_attrs(attr_str):
    """Parse GFF3 column 9 into OrderedDict."""
    attrs = OrderedDict()
    for pair in attr_str.split(";"):
        if "=" in pair:
            key, val = pair.split("=", 1)
            attrs[key.strip()] = val.strip()
    return attrs


def format_gff3_attrs(attrs):
    """Format dict → GFF3 column 9."""
    return ";".join(f"{k}={v}" for k, v in attrs.items() if v)


def build_hanno_to_locus_map(gff3_path, locus_prefix):
    """
    Parse GFF3 for gene features, assign sequential locus_tags in
    genomic order:
      1. scaffold_1 .. scaffold_N  (numerically, = by decreasing length)
      2. contigs ordered by decreasing length
      3. within each sequence, genes ordered by start position

    Returns: hanno_to_locus dict, list of gene tuples
    """
    # First pass: collect sequence lengths from ##sequence-region headers
    seq_lengths = {}
    genes = []

    with open(gff3_path) as f:
        for line in f:
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                if len(parts) == 4:
                    seq_lengths[parts[1]] = int(parts[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9 or parts[2] != "gene":
                continue
            attrs = parse_gff3_attrs(parts[8])
            gene_id = attrs.get("ID", "")
            hanno_name = gene_id.replace("gene-", "")
            gene_symbol = attrs.get("gene", "")
            genes.append((hanno_name, parts[0], int(parts[3]), int(parts[4]),
                          parts[6], gene_symbol))

    def seqid_sort_key(seqid):
        """
        Sort key: scaffolds first (by number), then contigs (by length desc).
        Returns tuple: (priority, sub_sort_value)
          - scaffold_N  → (0, N)           — scaffolds first, numerically
          - contig_N    → (1, -length)      — contigs next, longest first
          - anything else → (2, -length)
        """
        m_scaf = re.match(r'^scaffold_(\d+)$', seqid)
        if m_scaf:
            return (0, int(m_scaf.group(1)))

        # For contigs and everything else: sort by length descending
        length = seq_lengths.get(seqid, 0)
        return (1, -length)

    genes.sort(key=lambda g: (seqid_sort_key(g[1]), g[2]))

    hanno_to_locus = OrderedDict()
    for i, (hanno_name, *_rest) in enumerate(genes, 1):
        hanno_to_locus[hanno_name] = f"{locus_prefix}_{str(i).zfill(5)}"

    # Report ordering
    import sys
    first_contig_idx = None
    for i, (hanno_name, seqid, *_) in enumerate(genes, 1):
        if not seqid.startswith("scaffold_") and first_contig_idx is None:
            first_contig_idx = i
            break
    if first_contig_idx:
        print(f"    Numbering: scaffold genes PRIEUP_00001–PRIEUP_{str(first_contig_idx-1).zfill(5)}, "
              f"then contigs PRIEUP_{str(first_contig_idx).zfill(5)}–PRIEUP_{str(len(genes)).zfill(5)}",
              file=sys.stderr)

    return hanno_to_locus, genes


def build_genbank_location(locations, strand):
    """Build GenBank-style location string from sorted (start, end) tuples."""
    if not locations:
        return ""
    locs = sorted(locations)
    parts = [f"{s}..{e}" for s, e in locs]
    loc = parts[0] if len(parts) == 1 else f"join({','.join(parts)})"
    if strand == "-":
        loc = f"complement({loc})"
    return loc


# =============================================================================
# GFF3 reformatting — Zenodo mode (clean, NCBI FTP-like)
# =============================================================================
def reformat_gff3_zenodo(input_path, output_path, seqid_map, hanno_to_locus, source):
    """
    Produce NCBI FTP-style GFF3 with:
      - INSDC accessions as seqids
      - Clean locus_tag-based IDs (no gnl|)
      - gbkey attributes for NCBI compatibility
    """
    gene_metadata = {}

    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith("##gff-version"):
                fout.write(line)
                fout.write("#!processor HANNO annotation pipeline (Kopp 2025)\n")
                fout.write("#!genome-build aPriEup1.0\n")
                fout.write("#!genome-build-accession GCA_965278355.2\n")
                continue
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                if len(parts) == 4:
                    new_id = seqid_map.get(parts[1], parts[1])
                    fout.write(f"##sequence-region {new_id} {parts[2]} {parts[3]}\n")
                else:
                    fout.write(line)
                continue
            if line.startswith("#"):
                fout.write(line)
                continue
            if not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue

            # Map seqid
            parts[0] = seqid_map.get(parts[0], parts[0])
            parts[1] = source
            attrs = parse_gff3_attrs(parts[8])
            ftype = parts[2]

            if ftype == "gene":
                hanno = attrs.get("ID", "").replace("gene-", "")
                lt = hanno_to_locus.get(hanno, hanno)
                new = OrderedDict([
                    ("ID", f"gene-{lt}"),
                ])
                if attrs.get("gene"):
                    new["gene"] = attrs["gene"]
                new["locus_tag"] = lt
                new["gene_biotype"] = attrs.get("gene_biotype", "protein_coding")
                new["gbkey"] = "Gene"
                parts[8] = format_gff3_attrs(new)

                gene_metadata[lt] = {
                    "seqid": parts[0], "start": int(parts[3]), "end": int(parts[4]),
                    "strand": parts[6], "gene": attrs.get("gene", ""),
                    "product": "", "hanno_name": hanno,
                    "exon_locations": [], "cds_locations": [],
                }

            elif ftype == "mRNA":
                hanno = attrs.get("ID", "").replace("rna-", "")
                lt = hanno_to_locus.get(hanno, hanno)
                new = OrderedDict([
                    ("ID", f"rna-{lt}"),
                    ("Parent", f"gene-{lt}"),
                ])
                if attrs.get("gene"):
                    new["gene"] = attrs["gene"]
                new["locus_tag"] = lt
                new["product"] = attrs.get("product", "hypothetical protein")
                new["gbkey"] = "mRNA"
                parts[8] = format_gff3_attrs(new)

                if lt in gene_metadata:
                    gene_metadata[lt]["product"] = unquote(attrs.get("product", "hypothetical protein"))
                    if attrs.get("gene"):
                        gene_metadata[lt]["gene"] = attrs["gene"]

            elif ftype == "exon":
                parent_hanno = attrs.get("Parent", "").replace("rna-", "")
                lt = hanno_to_locus.get(parent_hanno, parent_hanno)
                old_id = attrs.get("ID", "")
                m = re.search(r'\.(\d+)$', old_id)
                enum = m.group(1) if m else "1"
                new = OrderedDict([
                    ("ID", f"exon-{lt}.{enum}"),
                    ("Parent", f"rna-{lt}"),
                    ("gbkey", "mRNA"),
                ])
                parts[8] = format_gff3_attrs(new)

                if lt in gene_metadata:
                    gene_metadata[lt]["exon_locations"].append(
                        (int(parts[3]), int(parts[4])))

            elif ftype == "CDS":
                parent_hanno = attrs.get("Parent", "").replace("rna-", "")
                lt = hanno_to_locus.get(parent_hanno, parent_hanno)
                new = OrderedDict([
                    ("ID", f"cds-{lt}"),
                    ("Parent", f"rna-{lt}"),
                ])
                if attrs.get("gene"):
                    new["gene"] = attrs["gene"]
                new["locus_tag"] = lt
                new["product"] = attrs.get("product", "hypothetical protein")
                new["gbkey"] = "CDS"
                parts[8] = format_gff3_attrs(new)

                if lt in gene_metadata:
                    gene_metadata[lt]["cds_locations"].append(
                        (int(parts[3]), int(parts[4])))

            fout.write("\t".join(parts) + "\n")

    return gene_metadata


# =============================================================================
# GFF3 reformatting — table2asn mode (for NCBI submission)
# =============================================================================
def reformat_gff3_table2asn(input_path, output_path, hanno_to_locus, dbname, source):
    """
    Produce table2asn-ready GFF3 with:
      - ORIGINAL scaffold seqids (must match genome FASTA submitted to NCBI)
      - protein_id = gnl|dbname|locus_tag  (NCBI requirement, gets replaced
        with real accessions like CAN####### during processing)
      - transcript_id = gnl|dbname|mrna.locus_tag
    """
    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            if not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue

            # Keep original seqid — must match genome FASTA
            parts[1] = source
            attrs = parse_gff3_attrs(parts[8])
            ftype = parts[2]

            if ftype == "gene":
                hanno = attrs.get("ID", "").replace("gene-", "")
                lt = hanno_to_locus.get(hanno, hanno)
                new = OrderedDict([
                    ("ID", f"gene-{lt}"),
                ])
                if attrs.get("gene"):
                    new["gene"] = attrs["gene"]
                new["locus_tag"] = lt
                new["gene_biotype"] = attrs.get("gene_biotype", "protein_coding")
                parts[8] = format_gff3_attrs(new)

            elif ftype == "mRNA":
                hanno = attrs.get("ID", "").replace("rna-", "")
                lt = hanno_to_locus.get(hanno, hanno)
                new = OrderedDict([
                    ("ID", f"rna-{lt}"),
                    ("Parent", f"gene-{lt}"),
                ])
                if attrs.get("gene"):
                    new["gene"] = attrs["gene"]
                new["locus_tag"] = lt
                new["product"] = attrs.get("product", "hypothetical protein")
                new["transcript_id"] = f"gnl|{dbname}|mrna.{lt}"
                parts[8] = format_gff3_attrs(new)

            elif ftype == "exon":
                parent_hanno = attrs.get("Parent", "").replace("rna-", "")
                lt = hanno_to_locus.get(parent_hanno, parent_hanno)
                old_id = attrs.get("ID", "")
                m = re.search(r'\.(\d+)$', old_id)
                enum = m.group(1) if m else "1"
                new = OrderedDict([
                    ("ID", f"exon-{lt}.{enum}"),
                    ("Parent", f"rna-{lt}"),
                ])
                parts[8] = format_gff3_attrs(new)

            elif ftype == "CDS":
                parent_hanno = attrs.get("Parent", "").replace("rna-", "")
                lt = hanno_to_locus.get(parent_hanno, parent_hanno)
                new = OrderedDict([
                    ("ID", f"cds-{lt}"),
                    ("Parent", f"rna-{lt}"),
                ])
                if attrs.get("gene"):
                    new["gene"] = attrs["gene"]
                new["product"] = attrs.get("product", "hypothetical protein")
                new["protein_id"] = f"gnl|{dbname}|{lt}"
                new["transcript_id"] = f"gnl|{dbname}|mrna.{lt}"
                parts[8] = format_gff3_attrs(new)

            fout.write("\t".join(parts) + "\n")


# =============================================================================
# FASTA reformatting (Zenodo mode only — clean headers)
# =============================================================================
def reformat_fasta(input_path, output_path, hanno_to_locus, gene_metadata, fasta_type):
    """
    Rewrite HANNO FASTA with NCBI FTP-style headers (no gnl|).

    fasta_type: "cds", "prot", or "rna"

    Output header format (matches real NCBI FTP files):
      >lcl|CBDIFN020000001.1_cds_PRIEUP_00001_1 [gene=SLC16A13]
       [locus_tag=PRIEUP_00001] [protein=solute carrier family 16 member 13]
       [location=join(15843..16041,24124..24267)] [gbkey=CDS]
    """
    seq_counter = {}  # per-seqid counter
    total = written = 0

    with open(input_path) as fin, open(output_path, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                total += 1
                raw_id = line[1:].split()[0]
                if raw_id.endswith("(+)") or raw_id.endswith("(-)"):
                    raw_id = raw_id[:-3]

                lt = hanno_to_locus.get(raw_id)
                if not lt or lt not in gene_metadata:
                    keep = False
                    continue

                meta = gene_metadata[lt]
                seqid = meta["seqid"]
                gene_sym = meta.get("gene", "")
                product = meta.get("product", "hypothetical protein")
                strand = meta.get("strand", "+")

                # Locations
                if fasta_type in ("cds", "prot"):
                    locs = meta.get("cds_locations", [])
                    gbkey = "CDS"
                else:
                    locs = meta.get("exon_locations", [])
                    gbkey = "mRNA"

                location = build_genbank_location(locs, strand)

                # Per-scaffold sequence counter
                seq_counter[seqid] = seq_counter.get(seqid, 0) + 1
                seq_num = seq_counter[seqid]

                # Build header — clean, no gnl|
                header = f">lcl|{seqid}_{fasta_type}_{lt}_{seq_num}"
                quals = []
                if gene_sym:
                    quals.append(f"[gene={gene_sym}]")
                quals.append(f"[locus_tag={lt}]")
                quals.append(f"[protein={product}]")
                if location:
                    quals.append(f"[location={location}]")
                quals.append(f"[gbkey={gbkey}]")

                fout.write(header + " " + " ".join(quals) + "\n")
                keep = True
                written += 1
            else:
                if keep:
                    fout.write(line)

    return total, written


# =============================================================================
# Traceability mapping
# =============================================================================
def write_id_mapping(hanno_to_locus, gene_metadata, output_path):
    with open(output_path, "w") as f:
        f.write("# HANNO → NCBI ID mapping for P. euphronides (GCA_965278355.2)\n")
        f.write("# hanno_id\tlocus_tag\tncbi_seqid\tgene_symbol\tproduct\n")
        for hanno_name, lt in hanno_to_locus.items():
            meta = gene_metadata.get(lt, {})
            f.write(f"{hanno_name}\t{lt}\t{meta.get('seqid','')}\t"
                    f"{meta.get('gene','')}\t{meta.get('product','hypothetical protein')}\n")


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.assembly_prefix

    print("=" * 70, file=sys.stderr)
    print("HANNO → NCBI Annotation Reformatter", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print(f"  Mode:           {args.mode}", file=sys.stderr)
    print(f"  Locus prefix:   {args.locus_prefix}", file=sys.stderr)
    if args.mode in ("table2asn", "both"):
        print(f"  DB name:        {args.dbname} (table2asn protein_ids only)", file=sys.stderr)
    print(f"  Output dir:     {args.outdir}", file=sys.stderr)
    print("=" * 70, file=sys.stderr)

    # --- Load mappings ---
    print("[1] Loading scaffold → accession mapping...", file=sys.stderr)
    seqid_map = load_seqid_map(args.seqid_map)
    print(f"    {len(seqid_map)} scaffolds mapped", file=sys.stderr)

    print("[2] Building hanno → locus_tag mapping...", file=sys.stderr)
    hanno_to_locus, genes = build_hanno_to_locus_map(args.gff3, args.locus_prefix)
    print(f"    {len(hanno_to_locus)} genes: {args.locus_prefix}_00001 — "
          f"{args.locus_prefix}_{str(len(hanno_to_locus)).zfill(5)}", file=sys.stderr)

    gene_metadata = {}

    # --- Zenodo output ---
    if args.mode in ("zenodo", "both"):
        zenodo_dir = args.outdir if args.mode == "zenodo" else os.path.join(args.outdir, "zenodo")
        os.makedirs(zenodo_dir, exist_ok=True)

        print("[3] Writing Zenodo GFF3 (NCBI FTP-style, with accessions)...", file=sys.stderr)
        gff_out = os.path.join(zenodo_dir, f"{prefix}_genomic.gff")
        gene_metadata = reformat_gff3_zenodo(
            args.gff3, gff_out, seqid_map, hanno_to_locus, args.source)
        n_named = sum(1 for m in gene_metadata.values() if m.get("gene"))
        print(f"    {len(gene_metadata)} genes ({n_named} with gene symbols) → {gff_out}",
              file=sys.stderr)

        # FASTAs
        step = 4
        for fasta_arg, ftype, suffix, label in [
            (args.cds, "cds", f"{prefix}_cds_from_genomic.fna", "CDS"),
            (args.protein, "prot", f"{prefix}_protein.faa", "Protein"),
            (args.mrna, "rna", f"{prefix}_rna_from_genomic.fna", "mRNA"),
        ]:
            if fasta_arg:
                print(f"[{step}] Writing Zenodo {label} FASTA...", file=sys.stderr)
                out = os.path.join(zenodo_dir, suffix)
                tot, kept = reformat_fasta(fasta_arg, out, hanno_to_locus,
                                           gene_metadata, ftype)
                print(f"    {tot} → {kept} sequences → {out}", file=sys.stderr)
            step += 1

        # ID mapping
        map_out = os.path.join(zenodo_dir, "hanno_to_ncbi_id_mapping.tsv")
        write_id_mapping(hanno_to_locus, gene_metadata, map_out)
        print(f"[{step}] ID mapping → {map_out}", file=sys.stderr)

        # Gzip all Zenodo files (NCBI FTP convention)
        step += 1
        print(f"[{step}] Gzipping Zenodo files (NCBI FTP convention)...", file=sys.stderr)
        import gzip as gzip_mod
        import shutil
        gz_count = 0
        for fname in os.listdir(zenodo_dir):
            fpath = os.path.join(zenodo_dir, fname)
            if os.path.isfile(fpath) and not fname.endswith(".gz"):
                with open(fpath, 'rb') as f_in:
                    with gzip_mod.open(fpath + '.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                gz_count += 1
                print(f"    {fname} → {fname}.gz", file=sys.stderr)

    # --- table2asn output ---
    if args.mode in ("table2asn", "both"):
        t2a_dir = args.outdir if args.mode == "table2asn" else os.path.join(args.outdir, "table2asn")
        os.makedirs(t2a_dir, exist_ok=True)

        print(f"[T] Writing table2asn GFF3 (original seqids, gnl|{args.dbname}| protein_ids)...",
              file=sys.stderr)
        t2a_out = os.path.join(t2a_dir, f"{prefix}_annotation.gff")
        reformat_gff3_table2asn(args.gff3, t2a_out, hanno_to_locus, args.dbname, args.source)
        print(f"    → {t2a_out}", file=sys.stderr)

        # If we didn't run zenodo mode, still need gene_metadata for the mapping
        if not gene_metadata:
            gene_metadata = reformat_gff3_zenodo(
                args.gff3, os.devnull, seqid_map, hanno_to_locus, args.source)

        map_out = os.path.join(t2a_dir, "hanno_to_ncbi_id_mapping.tsv")
        write_id_mapping(hanno_to_locus, gene_metadata, map_out)

    # --- Summary ---
    print("", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print("DONE", file=sys.stderr)
    print("=" * 70, file=sys.stderr)

    if args.mode in ("zenodo", "both"):
        zd = args.outdir if args.mode == "zenodo" else os.path.join(args.outdir, "zenodo")
        print(f"\n  Zenodo / reviewer files (gzipped, NCBI FTP-style):", file=sys.stderr)
        for fname in sorted(os.listdir(zd)):
            if fname.endswith(".gz"):
                fpath = os.path.join(zd, fname)
                size_mb = os.path.getsize(fpath) / (1024 * 1024)
                print(f"    {fpath}  ({size_mb:.1f} MB)", file=sys.stderr)

    if args.mode in ("table2asn", "both"):
        td = args.outdir if args.mode == "table2asn" else os.path.join(args.outdir, "table2asn")
        print(f"\n  table2asn submission GFF3 (gnl| protein_ids, original seqids):", file=sys.stderr)
        print(f"    {td}/{prefix}_annotation.gff", file=sys.stderr)
        print(f"\n  table2asn command:", file=sys.stderr)
        print(f"    table2asn -M n -J -c w -euk \\", file=sys.stderr)
        print(f"      -t template.sbt \\", file=sys.stderr)
        print(f"      -gaps-min 10 \\", file=sys.stderr)
        print(f"      -j \"[organism=Pristimantis euphronides] [geo_loc_name=Grenada]\" \\",
              file=sys.stderr)
        print(f"      -i genome.fsa \\", file=sys.stderr)
        print(f"      -f {td}/{prefix}_annotation.gff \\", file=sys.stderr)
        print(f"      -o output.sqn -Z -V b", file=sys.stderr)

    print("=" * 70, file=sys.stderr)


if __name__ == "__main__":
    main()
