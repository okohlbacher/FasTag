#!/usr/bin/env python3
"""Prune an NCBI taxdump to the lineages a FasTag taxdb actually references.

The full dump is ~500 MB across nodes.dmp and names.dmp and describes ~2.5M
taxa. A FasTag index covers a curated handful, and the classifier only ever
walks from those leaves to the root, so everything off those paths is dead
weight in a release.

Kept:
  nodes.dmp  every taxid on a leaf->root path (the ranks the rollup climbs)
  names.dmp  the 'scientific name' row for those taxids -- the only name class
             TaxStats::Taxonomy::load() reads

Usage:
  prune_taxdump.py <nodes.dmp> <names.dmp> <out_dir> <taxid> [taxid ...]
  prune_taxdump.py <nodes.dmp> <names.dmp> <out_dir> --from-proteomes <dir>
"""
import os
import sys


def split_dmp(line):
    # rows are "a\t|\tb\t|\tc\t|\t...\t|\n"
    return [f.strip() for f in line.rstrip("\n").rstrip("|").split("\t|")]


def main(argv):
    if len(argv) < 5:
        print(__doc__, file=sys.stderr)
        return 2
    nodes_p, names_p, out_dir = argv[1], argv[2], argv[3]
    rest = argv[4:]

    if rest[0] == "--from-proteomes":
        d = rest[1]
        leaves = {int(f.split(".")[0]) for f in os.listdir(d) if f.endswith(".fasta")}
    else:
        leaves = {int(x) for x in rest}
    if not leaves:
        print("no taxids given", file=sys.stderr)
        return 2

    parent, rank = {}, {}
    with open(nodes_p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = split_dmp(line)
            if len(f) < 3:
                continue
            try:
                t, p = int(f[0]), int(f[1])
            except ValueError:
                continue
            parent[t], rank[t] = p, f[2]

    # Walk every leaf to the root; the union is what the rollup can ever touch.
    keep = set()
    for leaf in leaves:
        t = leaf
        seen = set()
        while t in parent and t not in seen:
            seen.add(t)
            keep.add(t)
            if parent[t] == t:  # root is its own parent
                break
            t = parent[t]
    keep.update(leaves)

    os.makedirs(out_dir, exist_ok=True)
    n_nodes = 0
    with open(os.path.join(out_dir, "nodes.dmp"), "w", encoding="utf-8") as out:
        for t in sorted(keep):
            if t not in parent:
                continue
            out.write(f"{t}\t|\t{parent[t]}\t|\t{rank[t]}\t|\n")
            n_nodes += 1

    n_names = 0
    with open(names_p, encoding="utf-8", errors="replace") as fh, \
         open(os.path.join(out_dir, "names.dmp"), "w", encoding="utf-8") as out:
        for line in fh:
            f = split_dmp(line)
            if len(f) < 4 or f[3] != "scientific name":
                continue
            try:
                t = int(f[0])
            except ValueError:
                continue
            if t in keep:
                out.write(f"{t}\t|\t{f[1]}\t|\t\t|\tscientific name\t|\n")
                n_names += 1

    print(f"{len(leaves)} leaves -> {n_nodes} nodes, {n_names} names", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
