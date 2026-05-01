#!/usr/bin/env python3
"""
01_four_species_synteny.py

Four-species comparative synteny analysis of the P. euphronides
scaffolded assembly. Identifies compound scaffolds containing material
from two or more D. ebraccatus chromosomes, detects breakpoint positions,
and cross-validates against E. pustulosus.

Input:
    Six PAF files from D-Genies (Minimap2) whole-genome alignments
    in input/dgenies_pairwise/.

Output (in output/):
    Tables:
        Table_synteny_all_pairwise.tsv
        Table_compound_chromosomes.tsv
        Table_breakpoints.tsv
        Table_karyotype_reconstruction.tsv
        Table_small_scaffolds.tsv
    Figures (PNG + PDF + TIFF, 300 dpi):
        Fig_breakpoint_scaffold_{1,3,4,6}
        Fig_heatmap_PriEup_vs_DenEbr
        Fig_heatmap_PriEup_vs_EngPus
        Fig_heatmap_EleCoq_vs_DenEbr
        Fig_heatmap_EleCoq_vs_EngPus
        Fig_heatmap_EngPus_vs_DenEbr

Dependencies: Python 3.8+, matplotlib

Author: Kopp K, Pristimantis euphronides genome project
"""


import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgba
import matplotlib.gridspec as gridspec


# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

SPECIES = {
    "PriEup": {
        "full": "Pristimantis euphronides",
        "n": 16, "2n": 32,
        "assembly": "GCA_965278355.2",
        "chr_prefix": "scaffold_",
        "note": "RagTag-scaffolded against E. coqui"
    },
    "EleCoq": {
        "full": "Eleutherodactylus coqui",
        "n": 13, "2n": 26,
        "assembly": "GCF_035609145.1",
        "chr_prefix": "NC_0898",
        "note": "scaffolding reference"
    },
    "EngPus": {
        "full": "Engystomops pustulosus",
        "n": 11, "2n": 22,
        "assembly": "GCF_040894005.1",
        "chr_prefix": "NC_0924",
        "note": "independent outgroup (Leptodactylidae)"
    },
    "DenEbr": {
        "full": "Dendropsophus ebraccatus",
        "n": 15, "2n": 30,
        "assembly": "GCF_027789765.1",
        "chr_prefix": "NC_0914",
        "note": "independent outgroup, highest n (Hylidae)"
    },
}

# Colour palette for DenEbr chromosomes (D. ebraccatus units)
DENEBR_COLOURS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
    "#dcbeff", "#9A6324", "#800000", "#aaffc3", "#808000",
]


# ─────────────────────────────────────────────────────────────────────────────
# PAF parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_paf(paf_file):
    """
    Parse a PAF file into a list of record dicts.

    PAF columns:
        0: query name,  1: query length,  2: query start,  3: query end,
        4: strand,       5: target name,   6: target length, 7: target start,
        8: target end,   9: matching bases, 10: alignment block length,
        11: mapping quality
    """
    records = []
    with open(paf_file) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            records.append({
                "qname": cols[0],  "qlen": int(cols[1]),
                "qstart": int(cols[2]), "qend": int(cols[3]),
                "strand": cols[4],
                "tname": cols[5],  "tlen": int(cols[6]),
                "tstart": int(cols[7]), "tend": int(cols[8]),
                "matches": int(cols[9]),
                "block_len": int(cols[10]),
                "mapq": int(cols[11]),
            })
    return records


def aggregate_chr_pairs(records, q_prefixes=None, t_prefixes=None):
    """
    Aggregate alignment block lengths per (query_seq, target_seq) pair.
    Optionally filter by accession prefixes.

    Returns:
        dict: {(qname, tname): {"total_aln": int, "n_blocks": int,
                                  "qlen": int, "tlen": int}}
    """
    agg = defaultdict(lambda: {"total_aln": 0, "n_blocks": 0,
                                "qlen": 0, "tlen": 0})
    for r in records:
        qn, tn = r["qname"], r["tname"]
        if q_prefixes and not any(qn.startswith(p) for p in q_prefixes):
            continue
        if t_prefixes and not any(tn.startswith(p) for p in t_prefixes):
            continue
        key = (qn, tn)
        agg[key]["total_aln"] += r["block_len"]
        agg[key]["n_blocks"] += 1
        agg[key]["qlen"] = r["qlen"]
        agg[key]["tlen"] = r["tlen"]
    return dict(agg)


# ─────────────────────────────────────────────────────────────────────────────
# Chromosome naming
# ─────────────────────────────────────────────────────────────────────────────

def build_chr_names(agg, prefix, species_code):
    """
    Build sorted chromosome name mapping: {accession: "Sp1", "Sp2", ...}
    sorted by chromosome size descending.
    """
    sizes = {}
    for (q, t), v in agg.items():
        if q.startswith(prefix):
            sizes[q] = max(sizes.get(q, 0), v["qlen"])
        if t.startswith(prefix):
            sizes[t] = max(sizes.get(t, 0), v["tlen"])
    ordered = sorted(sizes.keys(), key=lambda x: -sizes[x])
    return {acc: f"{species_code}{i+1}" for i, acc in enumerate(ordered)}, sizes


# ─────────────────────────────────────────────────────────────────────────────
# Table writers
# ─────────────────────────────────────────────────────────────────────────────

def write_tsv(path, header, rows):
    """Write a TSV file with header and rows."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"  Written: {path} ({len(rows)} rows)")


# ─────────────────────────────────────────────────────────────────────────────
# Breakpoint analysis
# ─────────────────────────────────────────────────────────────────────────────

def compute_binned_synteny(records, scaffold, target_prefix, bin_size,
                           min_block=10000):
    """
    Bin alignment blocks along a scaffold and tally target chromosome
    coverage per bin.

    Returns:
        list of (bin_start, bin_end, {target_chr: total_block_len})
    """
    scf_len = 0
    blocks = []
    for r in records:
        if r["qname"] != scaffold:
            continue
        scf_len = r["qlen"]
        if not r["tname"].startswith(target_prefix):
            continue
        if r["block_len"] < min_block:
            continue
        blocks.append((r["qstart"], r["qend"], r["tname"], r["block_len"]))

    if scf_len == 0:
        return []

    n_bins = int(scf_len / bin_size) + 1
    bins = []
    for b in range(n_bins):
        bstart = b * bin_size
        bend = min((b + 1) * bin_size, scf_len)
        targets = defaultdict(float)
        for qs, qe, tn, bl in blocks:
            mid = (qs + qe) / 2
            if bstart <= mid < bend:
                targets[tn] += bl
        bins.append((bstart, bend, dict(targets)))
    return bins


def identify_breakpoints(bins, de_names, min_dominance_kb=200,
                         merge_distance_bins=3):
    """
    Walk along binned synteny and identify positions where the dominant
    DenEbr chromosome changes. Applies two filters:

    1. Minimum dominance: a bin's dominant target must have ≥min_dominance_kb
       of aligned sequence to count (filters noise from sparse bins).
    2. Merge nearby transitions: if the dominant target oscillates within
       merge_distance_bins, keep only the major transition.

    Returns:
        list of (bin_index, position_bp, from_chr, to_chr)
    """
    min_dom = min_dominance_kb * 1000  # convert to bytes

    # First pass: assign dominant target per bin (with min dominance)
    dom_per_bin = []
    for i, (bstart, bend, targets) in enumerate(bins):
        if not targets:
            dom_per_bin.append((i, bstart, None))
            continue
        dom_chr = max(targets, key=targets.get)
        if targets[dom_chr] < min_dom:
            dom_per_bin.append((i, bstart, None))
        else:
            dom_per_bin.append((i, bstart, dom_chr))

    # Second pass: smooth by majority vote in sliding window
    smoothed = []
    for i, (idx, bstart, dom) in enumerate(dom_per_bin):
        if dom is None:
            smoothed.append((idx, bstart, dom))
            continue
        # Look at window of ±merge_distance_bins
        window_start = max(0, i - merge_distance_bins)
        window_end = min(len(dom_per_bin), i + merge_distance_bins + 1)
        votes = defaultdict(float)
        for j in range(window_start, window_end):
            _, _, d = dom_per_bin[j]
            if d is not None:
                # Weight central bins more
                weight = merge_distance_bins + 1 - abs(j - i)
                votes[d] += weight
        majority = max(votes, key=votes.get) if votes else dom
        smoothed.append((idx, bstart, majority))

    # Third pass: identify transitions in smoothed series
    raw_breakpoints = []
    prev_dom = None
    for idx, bstart, dom in smoothed:
        if dom is None:
            continue
        if prev_dom is not None and dom != prev_dom:
            raw_breakpoints.append((idx, bstart, prev_dom, dom))
        prev_dom = dom

    # Merge breakpoints that are within merge_distance_bins of each other
    # (keeps the first occurrence of a genuine transition)
    if not raw_breakpoints:
        return []

    merged = [raw_breakpoints[0]]
    for bp in raw_breakpoints[1:]:
        last = merged[-1]
        if bp[0] - last[0] <= merge_distance_bins:
            # Replace if this represents the same overall transition direction
            if bp[3] != last[1]:  # overall from→to is different
                continue
        else:
            merged.append(bp)

    return merged


# ─────────────────────────────────────────────────────────────────────────────
# Figure: breakpoint/positional synteny along a compound scaffold
# ─────────────────────────────────────────────────────────────────────────────

def plot_breakpoint_figure(bins_de, bins_ep, scaffold, scf_len,
                           de_names, ep_names, de_colours,
                           breakpoints, outpath, bin_size):
    """
    Create a two-panel figure showing DenEbr and EngPus synteny along
    a compound scaffold, with breakpoint annotations.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.69, 4.0), sharex=True)

    for ax, bins, names, colours, ref_label in [
        (ax1, bins_de, de_names, de_colours, "D. ebraccatus"),
        (ax2, bins_ep, ep_names, None, "E. pustulosus"),
    ]:
        # Collect all targets in this alignment
        all_targets = set()
        for _, _, tgts in bins:
            all_targets.update(tgts.keys())

        # Sort targets by total alignment descending
        target_totals = defaultdict(float)
        for _, _, tgts in bins:
            for t, v in tgts.items():
                target_totals[t] += v
        sorted_targets = sorted(all_targets, key=lambda x: -target_totals[x])

        # Assign colours
        if colours:
            tgt_colours = {t: colours.get(t, "#cccccc") for t in sorted_targets}
        else:
            cmap = plt.cm.Set2
            tgt_colours = {t: cmap(i % 8) for i, t in enumerate(sorted_targets)}

        # Stacked bar chart
        bottoms = [0.0] * len(bins)
        for tgt in sorted_targets:
            values = []
            for _, _, tgts in bins:
                values.append(tgts.get(tgt, 0) / 1e6)
            positions = [b[0] / 1e6 for b in bins]
            widths = [(b[1] - b[0]) / 1e6 for b in bins]
            label = names.get(tgt, tgt[-8:])
            ax.bar(positions, values, width=widths, bottom=bottoms,
                   color=tgt_colours[tgt], label=label, align="edge",
                   edgecolor="none", linewidth=0)
            bottoms = [b + v for b, v in zip(bottoms, values)]

        ax.set_ylabel("Aligned (Mb)", fontsize=10)
        ax.set_title(f"vs {ref_label}", fontsize=11, fontstyle="italic", loc="left")
        ax.legend(fontsize=8, loc="upper right", ncol=min(len(sorted_targets), 4),
                  framealpha=0.8)
        ax.set_xlim(0, scf_len / 1e6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Mark breakpoints on both panels
    for bp_i, bp_pos, bp_from, bp_to in breakpoints:
        for ax in (ax1, ax2):
            ax.axvline(bp_pos / 1e6, color="red", linewidth=1.5,
                       linestyle="--", alpha=0.7, zorder=10)
        ax1.annotate(f"BP: {bp_pos/1e6:.0f} Mb\n{de_names.get(bp_from, '?')}→"
                     f"{de_names.get(bp_to, '?')}",
                     xy=(bp_pos / 1e6, ax1.get_ylim()[1]),
                     fontsize=8, color="red", ha="center", va="bottom")

    ax2.set_xlabel(f"Position along {scaffold} (Mb)", fontsize=10)
    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    # PDF
    fig.savefig(outpath.replace('.png', '.pdf'), bbox_inches="tight")
    # TIFF
    fig.savefig(outpath.replace('.png', '.tiff'), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {outpath} (.png, .pdf, .tiff)")


# ─────────────────────────────────────────────────────────────────────────────
# Figure: chromosome correspondence heatmap
# ─────────────────────────────────────────────────────────────────────────────

def plot_correspondence_matrix(agg_data, q_order, t_order, q_names, t_names,
                               q_label, t_label, title, outpath,
                               min_aln_mb=1.0):
    """
    Plot a heatmap of alignment block length between two species'
    chromosomes, annotated with Mb values.
    """
    n_q = len(q_order)
    n_t = len(t_order)

    matrix = [[0.0] * n_t for _ in range(n_q)]
    for i, q in enumerate(q_order):
        for j, t in enumerate(t_order):
            val = agg_data.get((q, t), {}).get("total_aln", 0) / 1e6
            matrix[i][j] = val

    fig, ax = plt.subplots(figsize=(6.69,
                                     max(4.0, n_q * 0.4 + 1.5)))

    im = ax.imshow(matrix, cmap="YlOrRd", aspect="auto",
                   interpolation="nearest")

    # Annotate cells with Mb values
    for i in range(n_q):
        for j in range(n_t):
            val = matrix[i][j]
            if val >= min_aln_mb:
                colour = "white" if val > max(max(r) for r in matrix) * 0.6 else "black"
                ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                        fontsize=7, color=colour, fontweight="bold")

    ax.set_xticks(range(n_t))
    ax.set_xticklabels([t_names.get(t, t[-6:]) for t in t_order],
                       rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n_q))
    ax.set_yticklabels([q_names.get(q, q) for q in q_order], fontsize=9)
    ax.set_xlabel(t_label, fontsize=11, fontstyle="italic")
    ax.set_ylabel(q_label, fontsize=11, fontstyle="italic")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Alignment block length (Mb)", fontsize=10)

    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    # PDF
    fig.savefig(outpath.replace('.png', '.pdf'), bbox_inches="tight")
    # TIFF
    fig.savefig(outpath.replace('.png', '.tiff'), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {outpath} (.png, .pdf, .tiff)")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

# ── Hardcoded input paths ──────────────────────────────────────────
BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYNTENY = os.path.join(BASEDIR, 'input', 'dgenies_pairwise')

PAF_PE_EC = os.path.join(SYNTENY, 'PriEup_EleCoq/map_Pristimantis_euphronides.genome_to_GCF_035609145.1_aEleCoq1.hap1_genomic.paf')
PAF_PE_EP = os.path.join(SYNTENY, 'PriEup_EngPus/map_Pristimantis_euphronides.genome_to_GCF_040894005.1_aEngPut4.maternal_genomic.paf')
PAF_PE_DE = os.path.join(SYNTENY, 'PriEup_DenEbr/map_Pristimantis_euphronides.genome_to_GCF_027789765.1_aDenEbr1.pat_genomic.paf')
PAF_EC_EP = os.path.join(SYNTENY, 'EleCoq_EngPus/map_GCF_035609145.1_aEleCoq1.hap1_genomic_to_GCF_040894005.1_aEngPut4.maternal_genomic.paf')
PAF_EC_DE = os.path.join(SYNTENY, 'EleCoq_DenEbr/map_GCF_035609145.1_aEleCoq1.hap1_genomic_to_GCF_027789765.1_aDenEbr1.pat_genomic.paf')
PAF_EP_DE = os.path.join(SYNTENY, 'EngPus_DenEbr/map_GCF_040894005.1_aEngPut4.maternal_genomic_to_GCF_027789765.1_aDenEbr1.pat_genomic.paf')

OUTDIR = os.path.join(BASEDIR, 'output')
MIN_ALN_MB = 5.0
BIN_SIZE_MB = 5.0


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    min_aln = MIN_ALN_MB * 1e6
    bin_size = int(BIN_SIZE_MB * 1e6)

    # ==================================================================
    # Step 1: Parse all 6 PAF files
    # ==================================================================
    print("Parsing PAF files...")
    paf_keys = {
        "pe_ec": PAF_PE_EC,
        "pe_ep": PAF_PE_EP,
        "pe_de": PAF_PE_DE,
        "ec_ep": PAF_EC_EP,
        "ec_de": PAF_EC_DE,
        "ep_de": PAF_EP_DE,
    }
    recs = {}
    for key, path in paf_keys.items():
        print(f"  {key}: {path}")
        recs[key] = parse_paf(path)
        print(f"    → {len(recs[key]):,} alignment records")

    # ==================================================================
    # Step 2: Aggregate and build chromosome name maps
    # ==================================================================
    print("\nAggregating chromosome-level alignments...")

    # Aggregate all pairs
    agg = {}
    for key in paf_keys:
        agg[key] = aggregate_chr_pairs(recs[key])

    # Build chromosome name maps from DenEbr targets (highest n = reference)
    de_names, de_sizes = build_chr_names(agg["ec_de"], "NC_0914", "De")
    ec_names, ec_sizes = build_chr_names(agg["ec_de"], "NC_0898", "Ec")
    ep_names, ep_sizes = build_chr_names(agg["ec_ep"], "NC_0924", "Ep")

    # PriEup scaffolds — name by size
    pe_sizes = {}
    for (q, _), v in agg["pe_ec"].items():
        if q.startswith("scaffold_"):
            pe_sizes[q] = max(pe_sizes.get(q, 0), v["qlen"])
    pe_names = {s: s for s in pe_sizes}  # keep original names

    # Sorted chromosome orders
    de_order = sorted([a for a in de_sizes if a.startswith("NC_")],
                      key=lambda x: -de_sizes[x])
    ec_order = sorted([a for a in ec_sizes if a.startswith("NC_")],
                      key=lambda x: -ec_sizes[x])
    ep_order = sorted([a for a in ep_sizes if a.startswith("NC_")],
                      key=lambda x: -ep_sizes[x])
    pe_order = sorted([s for s in pe_sizes if pe_sizes[s] > 1e6],
                      key=lambda x: -pe_sizes[x])  # large scaffolds only

    # DenEbr colour map
    de_colours = {}
    for i, acc in enumerate(de_order):
        de_colours[acc] = DENEBR_COLOURS[i % len(DENEBR_COLOURS)]

    # Print inventories
    print(f"\n  D. ebraccatus: {len(de_order)} chromosomes")
    for acc in de_order:
        print(f"    {de_names[acc]:>4s}: {acc}  {de_sizes[acc]/1e6:>7.1f} Mb")
    print(f"\n  E. coqui: {len(ec_order)} chromosomes")
    for acc in ec_order:
        print(f"    {ec_names[acc]:>4s}: {acc}  {ec_sizes[acc]/1e6:>7.1f} Mb")
    print(f"\n  E. pustulosus: {len(ep_order)} chromosomes")
    for acc in ep_order:
        print(f"    {ep_names[acc]:>4s}: {acc}  {ep_sizes[acc]/1e6:>7.1f} Mb")
    print(f"\n  P. euphronides: {len(pe_order)} large scaffolds (>{1} Mb)")

    # ==================================================================
    # Step 3: Pairwise alignment summary table (Supplementary)
    # ==================================================================
    print("\nWriting pairwise alignment summaries...")

    comparison_labels = {
        "ec_de": ("E. coqui", "D. ebraccatus", ec_order, de_order,
                  ec_names, de_names, "NC_0898", "NC_0914"),
        "ec_ep": ("E. coqui", "E. pustulosus", ec_order, ep_order,
                  ec_names, ep_names, "NC_0898", "NC_0924"),
        "ep_de": ("E. pustulosus", "D. ebraccatus", ep_order, de_order,
                  ep_names, de_names, "NC_0924", "NC_0914"),
        "pe_de": ("P. euphronides", "D. ebraccatus", pe_order, de_order,
                  pe_names, de_names, "scaffold_", "NC_0914"),
        "pe_ep": ("P. euphronides", "E. pustulosus", pe_order, ep_order,
                  pe_names, ep_names, "scaffold_", "NC_0924"),
        "pe_ec": ("P. euphronides", "E. coqui", pe_order, ec_order,
                  pe_names, ec_names, "scaffold_", "NC_0898"),
    }

    all_rows = []
    for key, (qlabel, tlabel, qord, tord, qn, tn, qpfx, tpfx) in \
            comparison_labels.items():
        for q in qord:
            for t in tord:
                pair = agg[key].get((q, t))
                if pair and pair["total_aln"] / 1e6 >= 1.0:
                    all_rows.append([
                        f"{qlabel} vs {tlabel}",
                        qn.get(q, q), f"{pair['qlen']/1e6:.1f}",
                        tn.get(t, t), f"{pair['tlen']/1e6:.1f}",
                        f"{pair['total_aln']/1e6:.1f}",
                        pair["n_blocks"],
                    ])

    write_tsv(
        os.path.join(OUTDIR,
                     "Table_synteny_all_pairwise.tsv"),
        ["comparison", "query_chr", "query_len_Mb",
         "target_chr", "target_len_Mb",
         "aligned_Mb", "n_blocks"],
        all_rows,
    )

    # ==================================================================
    # Step 4: Identify compound chromosomes
    # ==================================================================
    print("\nIdentifying compound chromosomes...")

    # 4a. E. coqui compounds (vs DenEbr)
    print("\n  E. coqui → D. ebraccatus:")
    ec_compounds = {}
    for ec in ec_order:
        hits = []
        for (q, t), v in agg["ec_de"].items():
            if q == ec and t.startswith("NC_0914") and v["total_aln"] >= min_aln:
                hits.append((t, v["total_aln"]))
        hits.sort(key=lambda x: -x[1])
        ec_compounds[ec] = hits
        de_str = " + ".join(f"{de_names[t]}({v/1e6:.0f}Mb)"
                            for t, v in hits)
        status = "COMPOUND" if len(hits) >= 2 else "simple"
        print(f"    {ec_names[ec]:>4s} ({ec_sizes[ec]/1e6:.0f} Mb): "
              f"{de_str}  [{status}]")

    # 4b. E. pustulosus compounds (vs DenEbr)
    print("\n  E. pustulosus → D. ebraccatus:")
    ep_compounds = {}
    for ep in ep_order:
        hits = []
        for (q, t), v in agg["ep_de"].items():
            if q == ep and t.startswith("NC_0914") and v["total_aln"] >= min_aln:
                hits.append((t, v["total_aln"]))
        hits.sort(key=lambda x: -x[1])
        ep_compounds[ep] = hits
        de_str = " + ".join(f"{de_names[t]}({v/1e6:.0f}Mb)"
                            for t, v in hits)
        status = "COMPOUND" if len(hits) >= 2 else "simple"
        print(f"    {ep_names[ep]:>4s} ({ep_sizes[ep]/1e6:.0f} Mb): "
              f"{de_str}  [{status}]")

    # 4c. P. euphronides scaffolds (vs DenEbr — independent of scaffolding)
    print("\n  P. euphronides → D. ebraccatus:")
    pe_compounds = {}
    compound_scaffolds = []
    for scf in pe_order:
        hits_de = []
        for (q, t), v in agg["pe_de"].items():
            if q == scf and t.startswith("NC_0914") and v["total_aln"] >= min_aln:
                hits_de.append((t, v["total_aln"]))
        hits_de.sort(key=lambda x: -x[1])
        pe_compounds[scf] = hits_de

        # Also get EngPus hits for cross-validation
        hits_ep = []
        for (q, t), v in agg["pe_ep"].items():
            if q == scf and t.startswith("NC_0924") and v["total_aln"] >= min_aln:
                hits_ep.append((t, v["total_aln"]))
        hits_ep.sort(key=lambda x: -x[1])

        # EleCoq (circular)
        hits_ec = []
        for (q, t), v in agg["pe_ec"].items():
            if q == scf and t.startswith("NC_0898") and v["total_aln"] >= min_aln:
                hits_ec.append((t, v["total_aln"]))
        hits_ec.sort(key=lambda x: -x[1])

        is_compound = len(hits_de) >= 2
        status = "COMPOUND" if is_compound else "simple"

        de_str = " + ".join(f"{de_names[t]}({v/1e6:.0f}Mb)"
                            for t, v in hits_de)
        ep_str = " + ".join(f"{ep_names[t]}({v/1e6:.0f}Mb)"
                            for t, v in hits_ep)
        ec_str = " + ".join(f"{ec_names[t]}({v/1e6:.0f}Mb)"
                            for t, v in hits_ec)

        print(f"    {scf:>14s} ({pe_sizes[scf]/1e6:.0f} Mb): [{status}]")
        print(f"      DenEbr: {de_str}")
        print(f"      EngPus: {ep_str}")
        print(f"      EleCoq: {ec_str}  (circular)")

        if is_compound:
            compound_scaffolds.append({
                "scaffold": scf,
                "size": pe_sizes[scf],
                "de_hits": hits_de,
                "ep_hits": hits_ep,
                "ec_hits": hits_ec,
            })

    # Write compound chromosome table
    compound_rows = []
    for cs in compound_scaffolds:
        scf = cs["scaffold"]
        de_str = " + ".join(de_names[t] for t, _ in cs["de_hits"])
        ep_str = " + ".join(ep_names[t] for t, _ in cs["ep_hits"])
        ec_str = " + ".join(ec_names[t] for t, _ in cs["ec_hits"])
        n_units = len(cs["de_hits"])

        for t, aln in cs["de_hits"]:
            compound_rows.append([
                scf, f"{cs['size']/1e6:.1f}",
                ec_str, de_names[t], t,
                f"{de_sizes.get(t, 0)/1e6:.0f}",
                f"{aln/1e6:.1f}",
                n_units,
                ep_str,
            ])

    write_tsv(
        os.path.join(OUTDIR,
                     "Table_compound_chromosomes.tsv"),
        ["PriEup_scaffold", "scaffold_Mb", "EleCoq_chr_circular",
         "DenEbr_unit", "DenEbr_accession", "DenEbr_size_Mb",
         "aligned_Mb", "n_D. ebraccatus_units",
         "EngPus_cross_validation"],
        compound_rows,
    )

    # ==================================================================
    # Step 5: Breakpoint analysis within compound scaffolds
    # ==================================================================
    print("\nAnalysing breakpoints in compound scaffolds...")

    breakpoint_rows = []
    for cs in compound_scaffolds:
        scf = cs["scaffold"]
        scf_len = cs["size"]

        # Binned synteny vs DenEbr
        bins_de = compute_binned_synteny(
            recs["pe_de"], scf, "NC_0914", bin_size)

        # Binned synteny vs EngPus
        bins_ep = compute_binned_synteny(
            recs["pe_ep"], scf, "NC_0924", bin_size)

        # Identify breakpoints from DenEbr bins
        bps = identify_breakpoints(bins_de, de_names)

        for bp_i, bp_pos, bp_from, bp_to in bps:
            # Check if EngPus also shows a transition nearby
            ep_bps = identify_breakpoints(bins_ep, ep_names)
            ep_confirmed = "no"
            for _, ep_pos, _, _ in ep_bps:
                if abs(ep_pos - bp_pos) < bin_size * 3:
                    ep_confirmed = "yes"
                    break

            breakpoint_rows.append([
                scf, f"{scf_len/1e6:.0f}",
                f"{bp_pos/1e6:.0f}",
                de_names.get(bp_from, bp_from),
                de_names.get(bp_to, bp_to),
                ep_confirmed,
            ])

            print(f"  {scf}: breakpoint at ~{bp_pos/1e6:.0f} Mb "
                  f"({de_names.get(bp_from, '?')} → "
                  f"{de_names.get(bp_to, '?')}), "
                  f"EngPus confirmed: {ep_confirmed}")

        # Plot breakpoint figure
        plot_breakpoint_figure(
            bins_de, bins_ep, scf, scf_len,
            de_names, ep_names, de_colours,
            bps,
            os.path.join(OUTDIR,
                         f"Fig_breakpoint_{scf}.png"),
            bin_size,
        )

    write_tsv(
        os.path.join(OUTDIR, "Table_breakpoints.tsv"),
        ["scaffold", "scaffold_Mb", "breakpoint_Mb",
         "DenEbr_unit_proximal", "DenEbr_unit_distal",
         "EngPus_cross_validated"],
        breakpoint_rows,
    )

    # ==================================================================
    # Step 6: Karyotype reconstruction table
    # ==================================================================
    print("\nBuilding karyotype reconstruction...")

    karyo_rows = []
    for de_acc in de_order:
        de_label = de_names[de_acc]
        de_sz = de_sizes.get(de_acc, 0)

        # Which PriEup scaffolds carry this DenEbr unit?
        scf_hits = []
        for (q, t), v in agg["pe_de"].items():
            if t == de_acc and q.startswith("scaffold_") and \
                    v["total_aln"] >= min_aln:
                scf_hits.append((q, v["total_aln"]))
        scf_hits.sort(key=lambda x: -x[1])

        # Which EleCoq chr carries this?
        ec_hits = []
        for (q, t), v in agg["ec_de"].items():
            if t == de_acc and q.startswith("NC_0898") and \
                    v["total_aln"] >= min_aln:
                ec_hits.append((q, v["total_aln"]))
        ec_hits.sort(key=lambda x: -x[1])

        # Which EngPus chr?
        ep_hits = []
        for (q, t), v in agg["ep_de"].items():
            if t == de_acc and q.startswith("NC_0924") and \
                    v["total_aln"] >= min_aln:
                ep_hits.append((q, v["total_aln"]))
        ep_hits.sort(key=lambda x: -x[1])

        status = "split" if len(scf_hits) >= 2 else "intact"
        scf_str = "; ".join(f"{q} ({v/1e6:.0f} Mb)" for q, v in scf_hits)
        ec_str = "; ".join(f"{ec_names[q]} ({v/1e6:.0f} Mb)"
                           for q, v in ec_hits)
        ep_str = "; ".join(f"{ep_names[q]} ({v/1e6:.0f} Mb)"
                           for q, v in ep_hits)

        karyo_rows.append([
            de_label, de_acc, f"{de_sz/1e6:.0f}",
            scf_str or "—",
            len(scf_hits),
            status,
            ec_str or "—",
            ep_str or "—",
        ])

    write_tsv(
        os.path.join(OUTDIR,
                     "Table_karyotype_reconstruction.tsv"),
        ["DenEbr_unit", "DenEbr_accession", "DenEbr_Mb",
         "PriEup_scaffolds", "n_PriEup_scaffolds", "status",
         "EleCoq_chr", "EngPus_chr"],
        karyo_rows,
    )

    # ==================================================================
    # Step 7: Small scaffolds table
    # ==================================================================
    print("\nAnalysing small scaffolds (14–31)...")

    small_rows = []
    all_scf = sorted(
        [q for q in pe_sizes if q.startswith("scaffold_")],
        key=lambda x: int(x.split("_")[1])
    )
    for scf in all_scf:
        scf_num = int(scf.split("_")[1])
        if scf_num < 14:
            continue
        scf_sz = pe_sizes[scf]

        # Collect all chr hits across 3 references
        for ref_key, ref_names, ref_prefix in [
            ("pe_ec", ec_names, "NC_0898"),
            ("pe_ep", ep_names, "NC_0924"),
            ("pe_de", de_names, "NC_0914"),
        ]:
            hits = []
            for (q, t), v in agg[ref_key].items():
                if q == scf and t.startswith(ref_prefix):
                    hits.append((t, v["total_aln"]))
            hits.sort(key=lambda x: -x[1])
            for t, aln_val in hits[:3]:  # top 3 per reference
                ref_label = {"pe_ec": "EleCoq", "pe_ep": "EngPus",
                             "pe_de": "DenEbr"}[ref_key]
                small_rows.append([
                    scf, f"{scf_sz/1e3:.1f}",
                    ref_label, ref_names.get(t, t), t,
                    f"{aln_val/1e3:.1f}",
                ])

    write_tsv(
        os.path.join(OUTDIR, "Table_small_scaffolds.tsv"),
        ["scaffold", "scaffold_kb", "reference", "ref_chr",
         "ref_accession", "aligned_kb"],
        small_rows,
    )

    # ==================================================================
    # Step 8: Correspondence heatmap figures
    # ==================================================================
    print("\nGenerating chromosome correspondence figures...")

    # EleCoq vs DenEbr (main figure)
    agg_ec_de_chr = {(q, t): v for (q, t), v in agg["ec_de"].items()
                     if q.startswith("NC_0898") and t.startswith("NC_0914")}
    plot_correspondence_matrix(
        agg_ec_de_chr, ec_order, de_order, ec_names, de_names,
        "E. coqui", "D. ebraccatus",
        "Chromosome correspondence: E. coqui vs D. ebraccatus",
        os.path.join(OUTDIR,
                     "Fig_heatmap_EleCoq_vs_DenEbr.png"),
    )

    # EngPus vs DenEbr
    agg_ep_de_chr = {(q, t): v for (q, t), v in agg["ep_de"].items()
                     if q.startswith("NC_0924") and t.startswith("NC_0914")}
    plot_correspondence_matrix(
        agg_ep_de_chr, ep_order, de_order, ep_names, de_names,
        "E. pustulosus", "D. ebraccatus",
        "Chromosome correspondence: E. pustulosus vs D. ebraccatus",
        os.path.join(OUTDIR,
                     "Fig_heatmap_EngPus_vs_DenEbr.png"),
    )

    # PriEup vs DenEbr (main figure — non-circular)
    agg_pe_de_chr = {(q, t): v for (q, t), v in agg["pe_de"].items()
                     if q.startswith("scaffold_") and t.startswith("NC_0914")
                     and pe_sizes.get(q, 0) > 1e6}
    plot_correspondence_matrix(
        agg_pe_de_chr, pe_order, de_order, pe_names, de_names,
        "P. euphronides", "D. ebraccatus",
        "Chromosome correspondence: P. euphronides vs D. ebraccatus",
        os.path.join(OUTDIR,
                     "Fig_heatmap_PriEup_vs_DenEbr.png"),
    )

    # PriEup vs EngPus
    agg_pe_ep_chr = {(q, t): v for (q, t), v in agg["pe_ep"].items()
                     if q.startswith("scaffold_") and t.startswith("NC_0924")
                     and pe_sizes.get(q, 0) > 1e6}
    plot_correspondence_matrix(
        agg_pe_ep_chr, pe_order, ep_order, pe_names, ep_names,
        "P. euphronides", "E. pustulosus",
        "Chromosome correspondence: P. euphronides vs E. pustulosus",
        os.path.join(OUTDIR,
                     "Fig_heatmap_PriEup_vs_EngPus.png"),
    )

    # EleCoq vs EngPus
    agg_ec_ep_chr = {(q, t): v for (q, t), v in agg["ec_ep"].items()
                     if q.startswith("NC_0898") and t.startswith("NC_0924")}
    plot_correspondence_matrix(
        agg_ec_ep_chr, ec_order, ep_order, ec_names, ep_names,
        "E. coqui", "E. pustulosus",
        "Chromosome correspondence: E. coqui vs E. pustulosus",
        os.path.join(OUTDIR,
                     "Fig_heatmap_EleCoq_vs_EngPus.png"),
    )

    # ==================================================================
    # Step 9: Summary
    # ==================================================================
    n_ec_compound = sum(1 for h in ec_compounds.values() if len(h) >= 2)
    n_ep_compound = sum(1 for h in ep_compounds.values() if len(h) >= 2)
    n_pe_compound = len(compound_scaffolds)
    total_units_pe = sum(len(cs["de_hits"]) for cs in compound_scaffolds)
    simple_pe = len(pe_order) - n_pe_compound
    max_n = simple_pe + total_units_pe

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nCompound chromosomes (against D. ebraccatus n=15):")
    print(f"  E. coqui (n=13):       {n_ec_compound} compound, "
          f"{13 - n_ec_compound} simple")
    print(f"  E. pustulosus (n=11):  {n_ep_compound} compound, "
          f"{11 - n_ep_compound} simple")
    print(f"  P. euphronides scaffolds: {n_pe_compound} compound, "
          f"{simple_pe} simple")

    print(f"\nKaryotype reconstruction:")
    print(f"  Simple scaffolds:     {simple_pe} → {simple_pe} chromosomes")
    print(f"  Compound scaffolds:   {n_pe_compound} → "
          f"{total_units_pe} D. ebraccatus units")
    print(f"  Maximum predicted n:  {max_n}")
    print(f"  Known n (Schmid):     16")
    if max_n > 16:
        print(f"  Excess:               {max_n - 16} "
              f"(fusions shared with E. coqui)")
    print(f"\nBreakpoints identified: {len(breakpoint_rows)}")
    for row in breakpoint_rows:
        print(f"  {row[0]}: ~{row[2]} Mb ({row[3]} → {row[4]}), "
              f"EngPus confirmed: {row[5]}")

    print(f"\nOutput directory: {OUTDIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
