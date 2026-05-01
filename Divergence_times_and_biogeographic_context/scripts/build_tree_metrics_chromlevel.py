#!/usr/bin/env python3
"""
Compute tMRCA and topological distance from a focal species to all tips in a phylogeny.

Uses Portik et al. (2023) time-calibrated tree to calculate:
  - tMRCA_Ma: divergence time to most recent common ancestor
  - topo_edges: number of edges (nodes) separating two tips

Performance note: Bio.Phylo.common_ancestor() is O(n) per call, which becomes
prohibitively slow for 5000+ tips. Instead we pre-compute parent pointers and
node depths, then find LCAs in O(depth) time using the standard "climb to same
depth, then climb together" algorithm (Steel & Penny 1993).
Author: Kopp K, Pristimantis euphronides genome project
"""

import argparse
from Bio import Phylo


def read_species_list(path):
    """Load species names from TSV (first column). Returns lowercase set."""
    species = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # handle TSV: take first column only
            name = line.split("\t")[0].strip()
            if name:
                species.add(name.lower())
    return species


def build_parent_maps(tree):
    """
    Single DFS pass to build lookup dicts for the whole tree.
    
    Returns:
        parent: node -> parent node (root maps to None)
        depth:  node -> number of edges from root
        rdist:  node -> cumulative branch length from root
    """
    root = tree.root
    parent = {root: None}
    depth = {root: 0}
    rdist = {root: 0.0}

    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.clades or []:
            parent[child] = node
            depth[child] = depth[node] + 1
            bl = child.branch_length if child.branch_length else 0.0
            rdist[child] = rdist[node] + bl
            stack.append(child)

    return parent, depth, rdist


def find_lca(node_a, node_b, parent, depth):
    """
    Find lowest common ancestor by climbing parent pointers.
    Standard algorithm: equalize depths first, then climb in tandem.
    """
    da, db = depth[node_a], depth[node_b]
    
    # bring deeper node up to same level
    while da > db:
        node_a = parent[node_a]
        da -= 1
    while db > da:
        node_b = parent[node_b]
        db -= 1
    
    # climb together until we meet
    while node_a is not node_b:
        node_a = parent[node_a]
        node_b = parent[node_b]
    
    return node_a


def main():
    parser = argparse.ArgumentParser(
        description="Calculate tMRCA and topological distance to a focal species"
    )
    parser.add_argument("--tree", required=True,
                        help="Newick time tree (Portik et al. Supplementary_File_S3)")
    parser.add_argument("--chromosome_level_list", required=True,
                        help="TSV with species having chromosome-level assemblies")
    parser.add_argument("--focal", default="Pristimantis_euphronides",
                        help="Focal tip label (default: Pristimantis_euphronides)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    # load chromosome-level species for flagging
    chrom_species = read_species_list(args.chromosome_level_list)

    # parse tree
    tree = Phylo.read(args.tree, "newick")
    terminals = [t for t in tree.get_terminals() if t.name]
    tip_lookup = {t.name: t for t in terminals}

    if args.focal not in tip_lookup:
        raise SystemExit(f"Focal tip '{args.focal}' not found in tree")

    # build traversal structures
    parent, depth, rdist = build_parent_maps(tree)

    # tree height = max root-to-tip distance (ultrametric assumption)
    tree_height = max(rdist[t] for t in terminals)
    focal_node = tip_lookup[args.focal]

    # compute metrics for all tips
    with open(args.out, "w", encoding="utf-8") as out:
        out.write("tip\tspecies\ttMRCA_Ma\ttopo_edges\thas_chromosome_level_genome\n")
        
        for tip_label, tip_node in tip_lookup.items():
            mrca = find_lca(focal_node, tip_node, parent, depth)
            
            # topological distance = sum of edges on path through MRCA
            topo = depth[focal_node] + depth[tip_node] - 2 * depth[mrca]
            
            # tMRCA = tree height minus root-to-MRCA distance
            tmrca = tree_height - rdist[mrca]
            
            species_name = tip_label.replace("_", " ")
            has_chrom = "True" if species_name.lower() in chrom_species else "False"
            
            out.write(f"{tip_label}\t{species_name}\t{tmrca:.8f}\t{topo}\t{has_chrom}\n")


if __name__ == "__main__":
    main()
