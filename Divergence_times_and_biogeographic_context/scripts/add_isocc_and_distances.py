#!/usr/bin/env python3
"""
Add geographic context to the phylogenetic proximity dataset.

Joins AmphibiaWeb occurrence data (ISO country codes) with the tMRCA table,
computes great-circle distances to Grenada, and flags Caribbean/mainland species.

Geographic distance uses country centroids as a proxy - not ideal for widespread
species but sufficient for identifying broad biogeographic patterns.
Author: Kopp K, Pristimantis euphronides genome project
"""

import argparse
import math
import pandas as pd

# Caribbean island territories (ISO 3166-1 alpha-2)
CARIBBEAN_ISLANDS = {
    "AG", "AI", "AW", "BB", "BL", "BM", "BQ", "BS", "CU", "CW", "DM", "DO",
    "GD", "GP", "HT", "JM", "KN", "KY", "LC", "MF", "MQ", "MS", "PR", "SX",
    "TC", "TT", "UM", "VC", "VG", "VI"
}

# American mainland (excludes Caribbean islands)
AMERICAN_MAINLAND = {
    "AR", "BO", "BR", "BZ", "CA", "CL", "CO", "CR", "EC", "FK", "GF", "GL",
    "GT", "GY", "HN", "MX", "NI", "PA", "PE", "PM", "PY", "SR", "SV", "US",
    "UY", "VE"
}


def haversine(lat1, lon1, lat2, lon2):
    """Great-circle distance in km using the haversine formula (Sinnott 1984)."""
    R = 6371.0088  # Earth mean radius, km
    
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    
    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda/2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def parse_isocc(isocc_str):
    """Split comma-separated ISO codes, normalize to uppercase."""
    if not isinstance(isocc_str, str) or not isocc_str.strip():
        return []
    return [c.strip().upper() for c in isocc_str.split(",") if c.strip()]


def load_amphibiaweb_isocc(path):
    """
    Build species -> isocc mapping from AmphibiaWeb export.
    Handles duplicate entries by taking union of ISO codes.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    
    for col in ["genus", "species", "isocc"]:
        if col not in df.columns:
            raise SystemExit(f"Missing column '{col}' in {path}")
    
    df = df.fillna("")
    df["sci_name"] = (df["genus"].str.strip() + " " + df["species"].str.strip()).str.lower()
    df = df[df["sci_name"].str.len() > 1]  # filter empty names
    
    # aggregate: union of ISO codes for duplicate species entries
    def merge_codes(series):
        all_codes = set()
        for val in series:
            all_codes.update(parse_isocc(val))
        return ",".join(sorted(all_codes))
    
    return df.groupby("sci_name")["isocc"].apply(merge_codes).to_dict()


def load_country_coords(path):
    """Load ISO2 -> (lat, lon) mapping. Returns coords dict and Grenada coords."""
    df = pd.read_csv(path)
    
    for col in ["country", "latitude", "longitude"]:
        if col not in df.columns:
            raise SystemExit(f"Missing column '{col}' in {path}")
    
    df["country"] = df["country"].astype(str).str.upper().str.strip()
    df = df.dropna(subset=["latitude", "longitude"])
    
    coords = {row.country: (row.latitude, row.longitude) 
              for row in df.itertuples(index=False)}
    
    if "GD" not in coords:
        raise SystemExit("Grenada (GD) not found in countries.csv")
    
    return coords, coords["GD"]


def main():
    parser = argparse.ArgumentParser(
        description="Add ISO country codes and geographic distances to tMRCA table"
    )
    parser.add_argument("--tmrca_tsv", required=True)
    parser.add_argument("--amphibia_names", required=True,
                        help="AmphibiaWeb amphib_names.txt export")
    parser.add_argument("--countries_csv", required=True,
                        help="Country centroids with ISO2 codes")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # load input data
    df = pd.read_csv(args.tmrca_tsv, sep="\t", dtype=str)
    name_to_isocc = load_amphibiaweb_isocc(args.amphibia_names)
    iso_coords, grenada_coords = load_country_coords(args.countries_csv)
    gd_lat, gd_lon = grenada_coords

    # match species to ISO codes
    df["isocc"] = df["species"].str.lower().str.strip().map(name_to_isocc).fillna("")

    # compute geographic metrics
    min_dist_list = []
    is_caribbean_list = []
    is_mainland_list = []

    for isocc in df["isocc"]:
        codes = parse_isocc(isocc)
        
        # region classification
        has_caribbean = any(c in CARIBBEAN_ISLANDS for c in codes)
        # mainland only if NOT Caribbean (avoid double-flagging widespread species)
        has_mainland = (not has_caribbean) and any(c in AMERICAN_MAINLAND for c in codes)
        
        is_caribbean_list.append(str(has_caribbean))
        is_mainland_list.append(str(has_mainland))
        
        # minimum distance to Grenada across all occurrence countries
        distances = []
        for code in codes:
            if code in iso_coords:
                lat, lon = iso_coords[code]
                distances.append(haversine(gd_lat, gd_lon, lat, lon))
        
        if distances:
            min_dist_list.append(f"{min(distances):.6f}")
        else:
            min_dist_list.append("")

    df["min_dist_km"] = min_dist_list
    df["is_Caribbean_island"] = is_caribbean_list
    df["is_American_mainland"] = is_mainland_list

    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} records to {args.out}")


if __name__ == "__main__":
    main()
