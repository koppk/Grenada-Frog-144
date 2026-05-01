#!/usr/bin/env python3
# Author: Kopp K, Pristimantis euphronides genome project
"""
Compute great-circle distances from Pristimantis euphronides type localities
to Grenada (Pointe Salines).

Refines the country-centroid distance analysis (see add_isocc_and_distances.py)
for the 12 closest mainland Pristimantis congeners of P. euphronides at the
next phylogenetic node beyond P. shrevei.

Per-species type-locality strings were retrieved from Amphibian Species of the
World v6.2 (https://amphibiansoftheworld.amnh.org). For species whose type
locality contains explicit coordinates (DMS or decimal), those are extracted
directly. For species without explicit coordinates, the most specific named
place was located on Google Maps and the resulting place URL was recorded; the
coordinates are extracted from the URL fragment that follows the place segment.

The reported location for each species is the coordinate pair, sourced either
from the original taxonomic description (where given) or from the Google Maps
URL.
"""

import argparse
import math
import re

import pandas as pd

# Pointe Salines, Grenada: southwesternmost point of the island and the closest
# single point of Grenada to the South American mainland.
POINTE_SALINES_LAT = 12.00
POINTE_SALINES_LON = -61.79

# DMS pattern: e.g. "05° 26′ N, 72° 44′ W" (allows variant unicode marks and
# straight quotes)
DMS_PATTERN = re.compile(
    r"(\d{1,2})\s*[\u00b0\u00ba]\s*(\d{1,2})\s*[\u2032\u2019\u0027\u2018]?\s*([NSns])"
    r"\s*[,;\s]+"
    r"(\d{1,3})\s*[\u00b0\u00ba]\s*(\d{1,2})\s*[\u2032\u2019\u0027\u2018]?\s*([WEwe])"
)

# Decimal pair pattern: e.g. "7.16866-72.2655" or "7.16866, -72.2655"
# Allows comma, hyphen, or en-dash as separator.
DECIMAL_PAIR_PATTERN = re.compile(
    r"(?<![\d.])([+-]?\d{1,2}\.\d+)\s*[,\u2013\-]\s*([+-]?\d{1,3}\.\d+)(?!\d)"
)

# Google Maps place URL: e.g. ".../place/<n>/@7.137,-72.668,17z/..."
GOOGLEMAPS_PLACE_PATTERN = re.compile(
    r"google\.com/maps/place/[^@]*@(-?\d+\.\d+),(-?\d+\.\d+)"
)

# ASW species URL: e.g. ".../Pristimantis/Pristimantis-anolirex"
ASW_PRISTIMANTIS_PATTERN = re.compile(
    r"amphibiansoftheworld\.amnh\.org/.+?/Pristimantis/Pristimantis-([a-z\-]+)"
)


def haversine(lat1, lon1, lat2, lon2):
    """Great-circle distance in km using the haversine formula (Sinnott 1984)."""
    R = 6371.0088  # Earth mean radius, km

    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda/2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def parse_dms(match):
    """Convert a DMS regex match to (lat, lon) decimal degrees."""
    d_lat, m_lat, hem_lat, d_lon, m_lon, hem_lon = match.groups()
    lat = int(d_lat) + int(m_lat) / 60.0
    lon = int(d_lon) + int(m_lon) / 60.0
    if hem_lat.upper() == "S":
        lat = -lat
    if hem_lon.upper() == "W":
        lon = -lon
    return lat, lon


def parse_decimal_pair(match):
    """
    Convert a decimal-pair regex match to (lat, lon) decimal degrees.

    Some original descriptions report coordinates as "lat-lon" without explicit
    hemisphere markers, where the longitude magnitude makes the western
    hemisphere implicit (e.g. "7.16866-72.2655" for 7.16866 N, 72.2655 W in
    northern South America). When the second value is positive and large
    enough to be a longitude magnitude (>30) without a sign, treat it as W.
    """
    v1 = float(match.group(1))
    v2 = float(match.group(2))
    if v2 > 30:
        v2 = -v2
    return v1, v2


def split_into_species_blocks(text):
    """
    Split the ASW-with-Google-coords text file into per-species blocks.
    Each block starts with the AmphibiaSpeciesOfTheWorld species URL.
    """
    blocks = re.split(r"(?=https://amphibiansoftheworld\.amnh\.org/)", text)
    return [b.strip() for b in blocks if b.strip()]


def extract_species_record(block):
    """
    Extract one species record from a block of text.
    Returns a dict with species, latitude_deg, longitude_deg, coord_source,
    and source_url.
    """
    m_species = ASW_PRISTIMANTIS_PATTERN.search(block)
    if not m_species:
        return None
    species = "Pristimantis_" + m_species.group(1)

    m_tl = re.search(r'Type locality:\s*"?([^"\n]+)"?', block)
    type_locality = m_tl.group(1).strip(' "') if m_tl else ""

    # 1. Try DMS in the Type locality string
    m = DMS_PATTERN.search(type_locality)
    if m:
        lat, lon = parse_dms(m)
        return {
            "species": species,
            "latitude_deg": lat,
            "longitude_deg": lon,
            "coord_source": "OriginalDescription_DMS",
            "source_url": "",
        }

    # 2. Try decimal pair in the Type locality string
    m = DECIMAL_PAIR_PATTERN.search(type_locality)
    if m:
        lat, lon = parse_decimal_pair(m)
        return {
            "species": species,
            "latitude_deg": lat,
            "longitude_deg": lon,
            "coord_source": "OriginalDescription_Decimal",
            "source_url": "",
        }

    # 3. Fall back to Google Maps URL within the block
    m = GOOGLEMAPS_PLACE_PATTERN.search(block)
    if m:
        lat = float(m.group(1))
        lon = float(m.group(2))
        # capture the full URL up to the first whitespace for the receipt
        url_match = re.search(r"https://www\.google\.com/maps/place/\S+", block)
        url = url_match.group(0) if url_match else ""
        return {
            "species": species,
            "latitude_deg": lat,
            "longitude_deg": lon,
            "coord_source": "GoogleMaps",
            "source_url": url,
        }

    return {
        "species": species,
        "latitude_deg": None,
        "longitude_deg": None,
        "coord_source": "MISSING",
        "source_url": "",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract type-locality coordinates from ASW + Google Maps text "
                    "and compute haversine distances to Pointe Salines, Grenada"
    )
    parser.add_argument("--asw_with_googlecoords", required=True,
                        help="Combined ASW + Google Maps URLs text file")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.asw_with_googlecoords, "r", encoding="utf-8") as f:
        text = f.read()

    blocks = split_into_species_blocks(text)
    records = []
    for block in blocks:
        record = extract_species_record(block)
        if record is None:
            continue
        if record["latitude_deg"] is None:
            print(f"WARNING: no coordinates found for {record['species']}")
            continue
        record["distance_to_pointesalines_km"] = haversine(
            POINTE_SALINES_LAT, POINTE_SALINES_LON,
            record["latitude_deg"], record["longitude_deg"],
        )
        records.append(record)

    df = pd.DataFrame(records, columns=[
        "species",
        "latitude_deg", "longitude_deg",
        "coord_source", "source_url",
        "distance_to_pointesalines_km",
    ])
    df["latitude_deg"] = df["latitude_deg"].map(lambda x: f"{x:.6f}")
    df["longitude_deg"] = df["longitude_deg"].map(lambda x: f"{x:.6f}")
    df["distance_to_pointesalines_km"] = df["distance_to_pointesalines_km"].map(
        lambda x: f"{x:.6f}"
    )

    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} records to {args.out}")


if __name__ == "__main__":
    main()
