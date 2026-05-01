#!/usr/bin/env python3
"""
generate_seqid_map.py — Generate a scaffold→NCBI accession mapping file.

For genomes submitted to ENA/NCBI, the submitted FASTA headers
(scaffold_1, scaffold_2, etc.) get assigned INSDC accessions
(CBDIFN020000001.1, CBDIFN020000002.1, etc.).

This script creates the mapping by either:
  (A) Matching scaffold names to NCBI accessions via the AGP or assembly report
  (B) Using the known WGS prefix and sequential numbering based on the FASTA order

Usage:
  # Option A: From NCBI assembly report (download from NCBI genome page)
  python3 generate_seqid_map.py \\
      --assembly-report GCA_965278355.2_aPriEup1.0_assembly_report.txt \\
      -o scaffold_to_accession.tsv

  # Option B: From your genome FASTA + WGS prefix
  python3 generate_seqid_map.py \\
      --genome genome.fasta \\
      --wgs-prefix CBDIFN02 \\
      --accession-version 1 \\
      -o scaffold_to_accession.tsv

  # Option C: From a two-column file you already have
  python3 generate_seqid_map.py \\
      --manual-map existing_mapping.tsv \\
      -o scaffold_to_accession.tsv

Author: K. Kopp, P. euphronides genome project
"""

import argparse
import sys
import os


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate scaffold → NCBI accession mapping file"
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--assembly-report",
                       help="NCBI assembly_report.txt (has both submitted and assigned names)")
    group.add_argument("--genome",
                       help="Genome FASTA file (scaffold names from headers)")
    group.add_argument("--manual-map",
                       help="Existing two-column TSV mapping (scaffold<TAB>accession)")

    p.add_argument("--wgs-prefix", default="CBDIFN02",
                   help="WGS accession prefix (default: CBDIFN02)")
    p.add_argument("--accession-version", default="1",
                   help="Accession version suffix (default: 1)")
    p.add_argument("--start-number", type=int, default=1,
                   help="Starting sequence number for WGS accessions (default: 1)")
    p.add_argument("-o", "--output", default="scaffold_to_accession.tsv",
                   help="Output TSV file")
    return p.parse_args()


def extract_scaffold_names_from_fasta(fasta_path):
    """Extract sequence IDs from FASTA file in order."""
    names = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                seqid = line[1:].split()[0].strip()
                names.append(seqid)
    return names


def parse_assembly_report(report_path):
    """
    Parse NCBI assembly_report.txt to get scaffold→accession mapping.

    The assembly report has columns:
      0: Sequence-Name
      1: Sequence-Role
      2: Assigned-Molecule
      3: Assigned-Molecule-loc/type
      4: GenBank-Accn
      5: Relationship
      6: RefSeq-Accn
      7: Assembly-Unit
      8: Sequence-Length
      9: UCSC-style-name
    """
    mapping = {}
    with open(report_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                scaffold_name = parts[0]
                genbank_accn = parts[4]
                if genbank_accn and genbank_accn != "na":
                    mapping[scaffold_name] = genbank_accn
    return mapping


def generate_wgs_mapping(scaffold_names, wgs_prefix, version, start_num):
    """
    Generate sequential WGS accessions for scaffolds.

    WGS accessions: {prefix}0000001.{version}, {prefix}0000002.{version}, ...
    The total digits in the number part depends on the prefix length.
    Standard ENA WGS: 8-char prefix + 7-digit number for 15 total base chars.
    CBDIFN02 (8 chars) + 7 digits = CBDIFN020000001
    """
    mapping = {}
    # Determine number of digits: 15 - len(prefix) = number width
    # But for CBDIFN02 with existing records CBDIFN020000001, it's 7 digits
    n_digits = 15 - len(wgs_prefix)
    if n_digits < 1:
        n_digits = 7  # fallback

    for i, name in enumerate(scaffold_names, start_num):
        acc_num = str(i).zfill(n_digits)
        accession = f"{wgs_prefix}{acc_num}.{version}"
        mapping[name] = accession

    return mapping


def main():
    args = parse_args()

    mapping = {}

    if args.assembly_report:
        print(f"Parsing assembly report: {args.assembly_report}", file=sys.stderr)
        mapping = parse_assembly_report(args.assembly_report)
        print(f"  Found {len(mapping)} scaffold→accession mappings", file=sys.stderr)

    elif args.genome:
        print(f"Reading scaffold names from: {args.genome}", file=sys.stderr)
        scaffolds = extract_scaffold_names_from_fasta(args.genome)
        print(f"  Found {len(scaffolds)} scaffolds", file=sys.stderr)
        print(f"  Generating WGS accessions with prefix {args.wgs_prefix}", file=sys.stderr)
        mapping = generate_wgs_mapping(
            scaffolds, args.wgs_prefix, args.accession_version, args.start_number
        )

    elif args.manual_map:
        print(f"Reading manual mapping from: {args.manual_map}", file=sys.stderr)
        with open(args.manual_map) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]
        print(f"  Loaded {len(mapping)} mappings", file=sys.stderr)

    # Write output
    with open(args.output, "w") as f:
        f.write("# Scaffold to NCBI accession mapping\n")
        f.write("# scaffold_name\tncbi_accession\n")
        for scaffold, accession in mapping.items():
            f.write(f"{scaffold}\t{accession}\n")

    print(f"Wrote {len(mapping)} mappings to {args.output}", file=sys.stderr)

    # Show first few for verification
    items = list(mapping.items())
    print(f"\nFirst 5 mappings:", file=sys.stderr)
    for scaffold, acc in items[:5]:
        print(f"  {scaffold}\t→\t{acc}", file=sys.stderr)
    if len(items) > 5:
        print(f"  ... ({len(items) - 5} more)", file=sys.stderr)


if __name__ == "__main__":
    main()
