#!/usr/bin/env python3
"""
Generate reference genome selection summary for scaffolding.

Compares phylogenetic proximity against chromosome-level genome availability
to identify the best reference and quantify the data gap for Strabomantidae.
Author: Kopp K, Pristimantis euphronides genome project
"""

import argparse
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Summarize reference genome selection based on tMRCA ranking"
    )
    parser.add_argument("--input", required=True,
                        help="TSV with tMRCA and chromosome-level genome flags")
    parser.add_argument("--focal", default="Pristimantis_euphronides")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    
    # validate columns
    required = ["tip", "species", "tMRCA_Ma", "topo_edges", "has_chromosome_level_genome"]
    missing = set(required) - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {missing}")

    df["tMRCA_Ma"] = pd.to_numeric(df["tMRCA_Ma"], errors="coerce")
    df["topo_edges"] = pd.to_numeric(df["topo_edges"], errors="coerce")

    # exclude focal species from comparisons
    df = df[df["tip"] != args.focal].copy()
    
    if df.empty:
        raise SystemExit("No species remaining after excluding focal")

    # split by genome availability
    has_genome = df["has_chromosome_level_genome"].str.lower() == "true"
    with_genome = df[has_genome].sort_values("tMRCA_Ma")
    without_genome = df[~has_genome].sort_values("tMRCA_Ma")

    if with_genome.empty:
        raise SystemExit("No chromosome-level genomes found")
    if without_genome.empty:
        raise SystemExit("All species have chromosome-level genomes (unexpected)")

    # closest in each category
    best_ref = with_genome.iloc[0]
    closest_relative = without_genome.iloc[0]

    # phylogenetic gap metrics
    gap_ratio = best_ref["tMRCA_Ma"] / closest_relative["tMRCA_Ma"]
    
    # how many species are closer than best reference but lack genomes?
    closer_no_genome = without_genome[without_genome["tMRCA_Ma"] < best_ref["tMRCA_Ma"]]
    n_closer = len(closer_no_genome)
    
    # count Pristimantis among those
    n_pristimantis = closer_no_genome["tip"].str.startswith("Pristimantis").sum()

    # write outputs
    os.makedirs(args.outdir, exist_ok=True)

    # text summary
    summary_file = os.path.join(args.outdir, "reference_selection_summary.txt")
    with open(summary_file, "w") as f:
        f.write("Reference Genome Selection Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Focal: {args.focal.replace('_', ' ')}\n")
        f.write(f"Source: Portik et al. (2023) time-calibrated phylogeny\n\n")

        f.write("BEST AVAILABLE REFERENCE (chromosome-level):\n")
        f.write(f"  {best_ref['species']}\n")
        f.write(f"  tMRCA: {best_ref['tMRCA_Ma']:.2f} Ma, {int(best_ref['topo_edges'])} edges\n\n")

        f.write("CLOSEST RELATIVE (no chromosome-level genome):\n")
        f.write(f"  {closest_relative['species']}\n")
        f.write(f"  tMRCA: {closest_relative['tMRCA_Ma']:.2f} Ma, {int(closest_relative['topo_edges'])} edges\n\n")

        f.write("PHYLOGENETIC DATA GAP:\n")
        f.write(f"  Best reference is {gap_ratio:.1f}x more distant than closest relative\n")
        f.write(f"  {n_closer} species closer than best reference lack chromosome-level genomes\n")
        f.write(f"    of which {n_pristimantis} are Pristimantis spp.\n\n")

        f.write("CHROMOSOME-LEVEL GENOMES RANKED BY tMRCA:\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'Species':<35} {'tMRCA (Ma)':>10} {'Edges':>6}\n")
        f.write("-" * 60 + "\n")
        for _, row in with_genome.head(10).iterrows():
            f.write(f"{row['species']:<35} {row['tMRCA_Ma']:>10.2f} {int(row['topo_edges']):>6}\n")

    print(f"Wrote {summary_file}")

    # TSV of ranked chromosome-level genomes
    ranked_file = os.path.join(args.outdir, "chromosome_level_ranked_by_tMRCA.tsv")
    with_genome.to_csv(ranked_file, sep="\t", index=False)
    print(f"Wrote {ranked_file}")

    # stdout summary
    print(f"\nBest reference: {best_ref['species']} ({best_ref['tMRCA_Ma']:.2f} Ma)")
    print(f"Closest relative without genome: {closest_relative['species']} ({closest_relative['tMRCA_Ma']:.2f} Ma)")
    print(f"Gap ratio: {gap_ratio:.1f}x | {n_closer} closer species lack genomes ({n_pristimantis} Pristimantis)")


if __name__ == "__main__":
    main()
