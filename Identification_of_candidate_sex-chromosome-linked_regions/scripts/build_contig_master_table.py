#!/usr/bin/env python3
"""
build_contig_master_table.py

Joins all contig-level evidence into a single master table.
One row per contig, all evidence columns, for placed + unplaced + satellite contigs.

Input files (all paths configurable at top):
  - Coverage classified: placed + unplaced + unmappable satellites
  - Seqkit per-contig: length + GC%
  - AGP: contig → scaffold mapping
  - WindowMasker: per-contig masked fraction
  - RepeatMasker parsed: per-contig TE breakdown
  - BAM QC Steps 1-4: CV, MAPQ, heterozygosity, clipping
  - Miniprot GFF: gene counts per contig
  - HANNO bedDB: gene counts per contig

Output:
  contig_master_table.tsv — one row per contig with all evidence

Usage:
  python3 build_contig_master_table.py

Author: K. Kopp, P. euphronides genome project
"""

import pandas as pd
import re
import sys
from pathlib import Path
from collections import Counter

# =============================================================================
# PATHS — adjust if needed
# =============================================================================
BASE = "/data/GrenadaFrog144/SexChromosomes/W-chr-workflow/contigComparison"

# Coverage
PLACED_COV   = f"{BASE}/coverage/placed_coverage_classified.tsv"
UNPLACED_COV = f"{BASE}/coverage/unplaced_coverage_classified.tsv"
SATELLITE    = f"{BASE}/coverage/unmappable_satellite_contigs.tsv"

# Seqkit (no headers: contig, length, GC%)
SEQKIT_PLACED   = f"{BASE}/placed_contigs.length_GC_cont.tsv"
SEQKIT_UNPLACED = f"{BASE}/unplaced_contigs.length_GC_cont.tsv"

# AGP
AGP = "/data/GrenadaFrog144/coverage/ragtag.scaffold.renamed.agp"

# WindowMasker
WM_PLACED   = f"{BASE}/windowmasker/wm_placed_per_contig.tsv"
WM_UNPLACED = f"{BASE}/windowmasker/wm_unplaced_per_contig.tsv"

# RepeatMasker parsed
RM_PLACED   = f"{BASE}/RepeatMasker_parsed/placed_TE_per_contig.tsv"
RM_UNPLACED = f"{BASE}/RepeatMasker_parsed/unplaced_TE_per_contig.tsv"

# BAM QC (no headers)
BAMQC_CV   = "/data/GrenadaFrog144/assembly_qc/contig_coverage_uniformity.tsv"
BAMQC_MAPQ = "/data/GrenadaFrog144/assembly_qc/contig_mapq.tsv"
BAMQC_HET  = "/data/GrenadaFrog144/assembly_qc/contig_heterozygosity.tsv"
BAMQC_CLIP = "/data/GrenadaFrog144/assembly_qc/contig_clipping_supplementary.tsv"

# Miniprot GFF
MP_PLACED   = f"{BASE}/placed_contigs_amphibia_miniprot.gff"
MP_UNPLACED = f"{BASE}/unplaced_contigs_amphibia_miniprot.gff"

# HANNO bedDB
HANNO_PLACED   = "/data/software/HANNO/WithmRNA_HANNO-RUN-placed_contigs.fasta/BESTMODELS-FINAL.bedDB"
HANNO_UNPLACED = "/data/software/HANNO/WithmRNA_HANNO-RUN-unplaced_contigs.fasta/BESTMODELS-FINAL.bedDB"

# Output
OUTDIR = f"{BASE}/master_table"
OUTFILE = f"{OUTDIR}/contig_master_table.tsv"

# Z-candidate regions (scaffold coordinates in bp)
Z_REGIONS = {
    "scaffold_2": (0, 42_000_000),
    "scaffold_8": (47_000_000, 64_000_000),
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_coverage():
    """Load placed + unplaced + satellite coverage into one DataFrame."""
    print("  Loading coverage files...")

    # Both files now have consistent headers:
    # contig, length, bases, mean_cov, min, max, coverage_class
    placed = pd.read_csv(PLACED_COV, sep="\t")
    placed = placed[["contig", "mean_cov", "coverage_class"]].copy()
    placed["set"] = "placed"

    unplaced = pd.read_csv(UNPLACED_COV, sep="\t")
    unplaced = unplaced[["contig", "mean_cov", "coverage_class"]].copy()
    unplaced["set"] = "unplaced"

    # Check: satellites should already be in unplaced_classified with coverage_class
    # The unmappable_satellite_contigs.tsv has the 67 zero-coverage contigs
    # They should be in unplaced_classified as coverage_class = "Zero" or similar
    # Let's verify and tag them
    sat = pd.read_csv(SATELLITE, sep="\t")
    sat.rename(columns={"chrom": "contig", "mean": "mean_cov"}, inplace=True)
    sat_contigs = set(sat["contig"].values)

    # Tag satellites in the unplaced set
    unplaced.loc[unplaced["contig"].isin(sat_contigs), "set"] = "satellite"

    # If any satellites are NOT in unplaced_classified, add them
    missing_sat = sat_contigs - set(unplaced["contig"].values)
    if missing_sat:
        print(f"    Adding {len(missing_sat)} satellite contigs not in unplaced_classified")
        sat_missing = sat[sat["contig"].isin(missing_sat)][["contig", "mean_cov"]].copy()
        sat_missing["coverage_class"] = "Zero"
        sat_missing["set"] = "satellite"
        unplaced = pd.concat([unplaced, sat_missing], ignore_index=True)

    df = pd.concat([placed, unplaced], ignore_index=True)
    print(f"    Placed: {(df['set']=='placed').sum()}, "
          f"Unplaced: {(df['set']=='unplaced').sum()}, "
          f"Satellite: {(df['set']=='satellite').sum()}")
    return df


def load_seqkit():
    """Load seqkit per-contig stats (no headers)."""
    print("  Loading seqkit stats...")
    frames = []
    for path in [SEQKIT_PLACED, SEQKIT_UNPLACED]:
        tmp = pd.read_csv(path, sep="\t", header=None,
                          names=["contig", "length", "gc_pct"])
        frames.append(tmp)
    df = pd.concat(frames, ignore_index=True)
    print(f"    {len(df)} contigs with length + GC")
    return df


def load_agp():
    """Parse AGP to get contig → scaffold mapping with positions."""
    print("  Loading AGP...")
    rows = []
    with open(AGP) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[4] != "W":  # skip gaps (U lines)
                continue
            scaffold = parts[0]
            # RagTag AGP includes ALL contigs — unplaced ones appear as
            # self-scaffolds (e.g. contig_15599 -> contig_15599).
            # Only keep entries on actual scaffolds (scaffold_1, scaffold_2, etc.)
            if not scaffold.startswith("scaffold_"):
                continue
            scaff_start = int(parts[1])
            scaff_end = int(parts[2])
            contig = parts[5]
            rows.append({
                "contig": contig,
                "scaffold": scaffold,
                "scaff_start": scaff_start,
                "scaff_end": scaff_end,
            })
    agp = pd.DataFrame(rows)

    # Assign Z-region flag
    def z_flag(row):
        scaffold = row["scaffold"]
        if scaffold not in Z_REGIONS:
            return "no"
        z_start, z_end = Z_REGIONS[scaffold]
        # Overlap: contig overlaps Z-region if scaff_start < z_end AND scaff_end > z_start
        overlap_start = max(row["scaff_start"], z_start)
        overlap_end = min(row["scaff_end"], z_end)
        overlap = max(0, overlap_end - overlap_start)
        contig_len = row["scaff_end"] - row["scaff_start"]
        if contig_len == 0:
            return "no"
        frac = overlap / contig_len
        if frac >= 0.5:
            return "Z-candidate"
        elif frac > 0:
            return "Z-boundary"
        else:
            return "no"

    agp["z_flag"] = agp.apply(z_flag, axis=1)
    n_z = (agp["z_flag"] == "Z-candidate").sum()
    n_b = (agp["z_flag"] == "Z-boundary").sum()
    print(f"    {len(agp)} placed contigs mapped to scaffolds")
    print(f"    Z-candidate: {n_z}, Z-boundary: {n_b}")
    return agp


def load_windowmasker():
    """Load WindowMasker per-contig masked fraction."""
    print("  Loading WindowMasker...")
    frames = []
    for path in [WM_PLACED, WM_UNPLACED]:
        tmp = pd.read_csv(path, sep="\t")
        # Header: contig, length, masked_bp, masked_fraction
        tmp = tmp[["contig", "masked_fraction"]].copy()
        tmp.rename(columns={"masked_fraction": "wm_masked_frac"}, inplace=True)
        frames.append(tmp)
    df = pd.concat(frames, ignore_index=True)
    print(f"    {len(df)} contigs")
    return df


def load_repeatmasker():
    """Load RepeatMasker parsed per-contig TE breakdown."""
    print("  Loading RepeatMasker...")
    te_classes = ["SINE", "LINE", "LTR", "DNA", "RC_Helitron",
                  "Satellite", "Simple_Low", "Unknown", "Other"]
    frames = []
    for path in [RM_PLACED, RM_UNPLACED]:
        # Header: contig, total_masked_bp, SINE, LINE, LTR, DNA, RC_Helitron,
        #         Satellite, Simple_Low, Unknown, Other, contig_length, masked_pct
        # Skip any duplicate header lines
        tmp = pd.read_csv(path, sep="\t")
        # Filter out any rows where 'contig' column literally says 'contig' (double header)
        tmp = tmp[tmp["contig"] != "contig"].copy()
        # Convert numeric columns
        for col in te_classes + ["total_masked_bp", "contig_length", "masked_pct"]:
            tmp[col] = pd.to_numeric(tmp[col], errors="coerce")
        frames.append(tmp)
    df = pd.concat(frames, ignore_index=True)

    # Determine dominant TE class per contig
    def dominant_te(row):
        vals = {c: row[c] for c in te_classes}
        if sum(vals.values()) == 0:
            return "none"
        return max(vals, key=vals.get)

    df["rm_dominant_te"] = df.apply(dominant_te, axis=1)

    # Rename for clarity
    rename_map = {"masked_pct": "rm_masked_pct", "total_masked_bp": "rm_total_masked_bp"}
    for c in te_classes:
        rename_map[c] = f"rm_{c}"
    df.rename(columns=rename_map, inplace=True)

    # Keep only what we need
    keep_cols = (["contig", "rm_total_masked_bp", "rm_masked_pct", "rm_dominant_te"]
                 + [f"rm_{c}" for c in te_classes])
    df = df[keep_cols].copy()
    print(f"    {len(df)} contigs with TE data")
    return df


def load_bamqc():
    """Load all 4 BAM QC steps.
    Detect whether header is present and handle both cases.
    """
    print("  Loading BAM QC...")

    def read_bamqc(path, expected_cols, col_names):
        """Read a BAM QC file, detecting header automatically."""
        with open(path) as fh:
            first_line = fh.readline().strip()
        # If first field starts with "contig" and contains a tab, it's a header
        if first_line.startswith("contig\t") or first_line.startswith("contig "):
            df = pd.read_csv(path, sep="\t")
        else:
            df = pd.read_csv(path, sep="\t", header=None, names=col_names)
        return df

    # Step 1: coverage uniformity — 5 cols
    cv = read_bamqc(BAMQC_CV, 5,
                     ["contig", "bqc_n_windows", "bqc_median_cov",
                      "bqc_mean_cov_windows", "bqc_cv"])
    # Normalise column names if header was present
    cv.columns = ["contig", "bqc_n_windows", "bqc_median_cov",
                   "bqc_mean_cov_windows", "bqc_cv"]
    print(f"    Step1 CV: {len(cv)} contigs")

    # Step 2: MAPQ — 5 cols
    mapq = read_bamqc(BAMQC_MAPQ, 5,
                       ["contig", "bqc_n_alignments", "bqc_mean_mapq",
                        "bqc_frac_mapq0", "bqc_frac_mapq_high"])
    mapq.columns = ["contig", "bqc_n_alignments", "bqc_mean_mapq",
                     "bqc_frac_mapq0", "bqc_frac_mapq_high"]
    print(f"    Step2 MAPQ: {len(mapq)} contigs")

    # Step 3: heterozygosity — 4 cols
    het = read_bamqc(BAMQC_HET, 4,
                      ["contig", "bqc_het_length", "bqc_n_variants",
                       "bqc_variants_per_kb"])
    het.columns = ["contig", "bqc_het_length", "bqc_n_variants",
                    "bqc_variants_per_kb"]
    print(f"    Step3 Het: {len(het)} contigs")

    # Step 4: clipping — 4 cols
    clip = read_bamqc(BAMQC_CLIP, 4,
                       ["contig", "bqc_clip_n_reads", "bqc_soft_clip_frac",
                        "bqc_supp_frac"])
    clip.columns = ["contig", "bqc_clip_n_reads", "bqc_soft_clip_frac",
                     "bqc_supp_frac"]
    print(f"    Step4 Clip: {len(clip)} contigs")

    return cv, mapq, het, clip


def count_miniprot_genes(gff_path):
    """Count mRNA features per contig from miniprot GFF, also extract gene names."""
    gene_counts = Counter()
    gene_names = {}  # contig → set of gene names
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "mRNA":
                continue
            contig = parts[0]
            gene_counts[contig] += 1
            # Extract gene name from Target= field
            attrs = parts[8]
            m = re.search(r'Target=sp\|([^|]+)\|(\S+)', attrs)
            if m:
                gene_id = m.group(2)  # e.g. NMDZ1_XENLA
                gene_name = gene_id.split("_")[0]  # e.g. NMDZ1
                if contig not in gene_names:
                    gene_names[contig] = set()
                gene_names[contig].add(gene_name)
    return gene_counts, gene_names


def load_miniprot():
    """Load miniprot gene counts per contig from both GFF files."""
    print("  Loading miniprot GFF...")
    all_counts = Counter()
    for path in [MP_PLACED, MP_UNPLACED]:
        counts, _ = count_miniprot_genes(path)
        all_counts.update(counts)
    df = pd.DataFrame([
        {"contig": c, "mp_gene_count": n}
        for c, n in all_counts.items()
    ])
    print(f"    {len(df)} contigs with miniprot hits")
    return df


def load_hanno():
    """Load HANNO gene counts per contig from both bedDB files."""
    print("  Loading HANNO bedDB...")
    total_counts = Counter()
    named_counts = Counter()
    for path in [HANNO_PLACED, HANNO_UNPLACED]:
        with open(path) as fh:
            for line in fh:
                if line.startswith("##") or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 26:
                    continue
                contig = parts[0]
                pref_name = parts[25]  # col26 = Preferred_name (0-indexed: 25)
                total_counts[contig] += 1
                if pref_name not in ("-", ""):
                    named_counts[contig] += 1

    # Build DataFrame from all contigs seen
    all_contigs = set(total_counts.keys())
    rows = []
    for c in all_contigs:
        rows.append({
            "contig": c,
            "hanno_total_genes": total_counts[c],
            "hanno_named_genes": named_counts.get(c, 0),
        })
    df = pd.DataFrame(rows)
    print(f"    {len(df)} contigs with HANNO gene models")
    return df


# =============================================================================
# MAIN
# =============================================================================

def main():
    Path(OUTDIR).mkdir(parents=True, exist_ok=True)
    print("Building contig master table...")

    # 1. Base: coverage
    master = load_coverage()

    # 2. Seqkit: length + GC
    seqkit = load_seqkit()
    master = master.merge(seqkit, on="contig", how="left")

    # 3. AGP: scaffold mapping + Z-flag (placed only)
    agp = load_agp()
    master = master.merge(agp, on="contig", how="left")

    # 4. WindowMasker
    wm = load_windowmasker()
    master = master.merge(wm, on="contig", how="left")

    # 5. RepeatMasker
    rm = load_repeatmasker()
    master = master.merge(rm, on="contig", how="left")

    # 6. BAM QC
    cv, mapq, het, clip = load_bamqc()
    master = master.merge(cv, on="contig", how="left")
    master = master.merge(mapq, on="contig", how="left")
    master = master.merge(het, on="contig", how="left")
    master = master.merge(clip, on="contig", how="left")

    # Force numeric on all BAM QC columns (handles non-numeric parse artifacts)
    for col in [c for c in master.columns if c.startswith("bqc_")]:
        master[col] = pd.to_numeric(master[col], errors="coerce")

    # 7. Miniprot
    mp = load_miniprot()
    master = master.merge(mp, on="contig", how="left")
    master["mp_gene_count"] = master["mp_gene_count"].fillna(0).astype(int)

    # 8. HANNO
    hanno = load_hanno()
    master = master.merge(hanno, on="contig", how="left")
    master["hanno_total_genes"] = master["hanno_total_genes"].fillna(0).astype(int)
    master["hanno_named_genes"] = master["hanno_named_genes"].fillna(0).astype(int)

    # Fill missing values for contigs with no TE hits
    for col in [c for c in master.columns if c.startswith("rm_") and c != "rm_dominant_te"]:
        master[col] = master[col].fillna(0)
    master["rm_dominant_te"] = master["rm_dominant_te"].fillna("none")

    # Fill z_flag for unplaced/satellite
    master["z_flag"] = master["z_flag"].fillna("NA")

    # Order columns logically
    col_order = [
        # Identity
        "contig", "set", "length", "gc_pct",
        # Coverage
        "mean_cov", "coverage_class",
        # Scaffold context (placed only)
        "scaffold", "scaff_start", "scaff_end", "z_flag",
        # Masking
        "wm_masked_frac", "rm_masked_pct", "rm_dominant_te",
        "rm_total_masked_bp",
        "rm_SINE", "rm_LINE", "rm_LTR", "rm_DNA", "rm_RC_Helitron",
        "rm_Satellite", "rm_Simple_Low", "rm_Unknown", "rm_Other",
        # BAM QC
        "bqc_cv", "bqc_median_cov", "bqc_mean_cov_windows", "bqc_n_windows",
        "bqc_mean_mapq", "bqc_frac_mapq0", "bqc_frac_mapq_high", "bqc_n_alignments",
        "bqc_variants_per_kb", "bqc_n_variants", "bqc_het_length",
        "bqc_soft_clip_frac", "bqc_supp_frac", "bqc_clip_n_reads",
        # Gene content
        "mp_gene_count", "hanno_total_genes", "hanno_named_genes",
    ]
    # Only include columns that exist
    col_order = [c for c in col_order if c in master.columns]
    master = master[col_order]

    # Sort: placed first (by scaffold, then position), then unplaced by coverage class
    def sort_key(row):
        if row["set"] == "placed":
            return (0, str(row["scaffold"]), int(row.get("scaff_start", 0) or 0))
        elif row["set"] == "unplaced":
            return (1, str(row["coverage_class"]), 0)
        else:  # satellite
            return (2, "", 0)

    master["_sort"] = master.apply(sort_key, axis=1)
    master.sort_values("_sort", inplace=True)
    master.drop("_sort", axis=1, inplace=True)

    # Write
    master.to_csv(OUTFILE, sep="\t", index=False)
    print(f"\nOutput: {OUTFILE}")
    print(f"Total contigs: {len(master)}")

    # Summary stats
    print("\n=== Summary by set ===")
    for s in ["placed", "unplaced", "satellite"]:
        sub = master[master["set"] == s]
        print(f"  {s}: {len(sub)} contigs, "
              f"{sub['length'].sum()/1e6:.1f} Mb, "
              f"median GC={sub['gc_pct'].median():.1f}%")

    print("\n=== Summary by coverage class ===")
    for cls in master["coverage_class"].unique():
        sub = master[master["coverage_class"] == cls]
        print(f"  {cls}: {len(sub)} contigs, "
              f"{sub['length'].sum()/1e6:.1f} Mb, "
              f"mean RM masked={sub['rm_masked_pct'].mean():.1f}%, "
              f"mean WM masked={sub['wm_masked_frac'].mean():.2f}, "
              f"mean het/kb={sub['bqc_variants_per_kb'].mean():.2f}")

    print("\n=== Z-candidate contigs (placed) ===")
    z = master[master["z_flag"] == "Z-candidate"]
    print(f"  {len(z)} contigs, {z['length'].sum()/1e6:.1f} Mb")
    print(f"  Mean coverage: {z['mean_cov'].mean():.2f}")
    print(f"  Mean het/kb: {z['bqc_variants_per_kb'].mean():.2f}")
    print(f"  Mean RM masked: {z['rm_masked_pct'].mean():.1f}%")

    print("\n=== Gene content by set ===")
    for s in ["placed", "unplaced", "satellite"]:
        sub = master[master["set"] == s]
        print(f"  {s}: HANNO={sub['hanno_total_genes'].sum()} "
              f"(named={sub['hanno_named_genes'].sum()}), "
              f"miniprot={sub['mp_gene_count'].sum()}")

    # W-candidate totals for Schmid comparison
    print("\n=== W-candidate totals (Schmid comparison) ===")
    w_genic = master[(master["set"].isin(["unplaced", "satellite"])) &
                     (master["coverage_class"] == "Hemi_0.5x")]
    w_hetero = master[(master["set"].isin(["unplaced", "satellite"])) &
                      (master["coverage_class"].isin(["High_Coverage/Repeat"]))]
    w_zero = master[master["set"] == "satellite"]
    w_auto = master[(master["set"].isin(["unplaced", "satellite"])) &
                    (master["coverage_class"] == "Auto_1.0x")]
    w_low = master[(master["set"].isin(["unplaced", "satellite"])) &
                   (master["coverage_class"] == "Low_Coverage/Other")]

    print(f"  W-genic (Hemi_0.5x unplaced):         {len(w_genic)} contigs, "
          f"{w_genic['length'].sum()/1e6:.1f} Mb")
    print(f"  W-heterochromatin (High_Cov/Repeat):   {len(w_hetero)} contigs, "
          f"{w_hetero['length'].sum()/1e6:.1f} Mb")
    print(f"  W-satellite (Zero coverage):           {len(w_zero)} contigs, "
          f"{w_zero['length'].sum()/1e6:.1f} Mb")
    print(f"  Unscaffolded autosomal (Auto_1.0x):    {len(w_auto)} contigs, "
          f"{w_auto['length'].sum()/1e6:.1f} Mb")
    print(f"  Low/Other:                             {len(w_low)} contigs, "
          f"{w_low['length'].sum()/1e6:.1f} Mb")
    w_total = w_genic["length"].sum() + w_hetero["length"].sum() + w_zero["length"].sum()
    print(f"  ---")
    print(f"  Total W-candidate (genic+hetero+sat):  "
          f"{(w_genic['length'].sum()+w_hetero['length'].sum()+w_zero['length'].sum())/1e6:.1f} Mb")
    print(f"  Schmid et al. 2002 W estimate:         ~500-550 Mb")
    print(f"  Schmid et al. 2002 Z estimate:         ~31 Mb")

    z_total = z["length"].sum() / 1e6 if len(z) > 0 else 0
    print(f"  Z-candidate total (placed):            {z_total:.1f} Mb")

    print("\nDone.")


if __name__ == "__main__":
    main()
