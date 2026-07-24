#!/usr/bin/env python3
"""Fetch UniProt reference proteomes for the FasTag species reference set.

  fetch_reference_set.py data/taxonomy/reference-set.tsv <out_dir>

Writes <out_dir>/<taxid>.fasta, which is what buildtaxdb expects.

REFERENCE PROTEOMES, not `reviewed:true`. Curation effort tracks how
well-studied an organism is, so a reviewed-only fetch gave 20,431 proteins for
human and 35 for S. pneumoniae -- and that lands directly in the background
model, where a taxon with 35 proteins has an artificially tiny expectation and
anything matching it looks enriched.

Resolution is two-step because UniProt has no single query for "the reference
proteome of taxid X": list the proteomes for the organism, then take the one
whose proteomeType is "Reference proteome" (an organism can also carry
"Excluded" or redundant entries).

Skips files that already exist, so an interrupted run resumes.
"""
import os
import sys
import time
import urllib.parse
import urllib.request

API = "https://rest.uniprot.org"
FTP = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release"
       "/knowledgebase/reference_proteomes")


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read()
        except Exception as e:  # transient 5xx / reset
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))
            print(f"    retry {i+1}: {e}", file=sys.stderr)
    return b""


def reference_proteome(taxid):
    """(UP accession, superkingdom, taxid_actually_fetched) or (None,None,None).

    Two-stage on purpose. `organism_id:` matches the taxid EXACTLY, and most
    bacterial SPECIES carry no reference proteome of their own -- the reference
    sits on a strain (S. aureus 1280 -> USA300 367830). An exact-match-only
    resolver silently returned nothing for 11 of 50 taxa here. Falling back to
    the taxonomy SUBTREE finds the strain, and we then key the file by the
    strain's taxid, which is what the sequences actually are; roll-up to genus
    is the classifier's job and is unaffected.
    """
    import json

    def pick(query):
        q = urllib.parse.quote(query)
        data = get(f"{API}/proteomes/search?query={q}&format=json&size=50")
        rs = json.loads(data).get("results", [])
        for r in rs:
            if r.get("proteomeType") == "Reference proteome":
                return r
        return None

    def kingdom(r):
        lin = [x.get("scientificName", "") if isinstance(x, dict) else str(x)
               for x in r.get("taxonomy", {}).get("lineage", [])]
        for k in ("Bacteria", "Archaea", "Viruses"):
            if k in lin:
                return k
        return "Eukaryota"

    r = pick(f"organism_id:{taxid}") or pick(f"taxonomy_id:{taxid}")
    if not r:
        return None, None, None
    got = r.get("taxonomy", {}).get("taxonId", taxid)
    return r.get("id"), kingdom(r), got


def main(argv):
    if len(argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    manifest, out_dir = argv[1], argv[2]
    os.makedirs(out_dir, exist_ok=True)

    entries = []
    for line in open(manifest, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) >= 2:
            entries.append((f[0], f[1]))

    total_prot = 0
    failures = []
    for i, (taxid, label) in enumerate(entries, 1):
        path = os.path.join(out_dir, f"{taxid}.fasta")
        if os.path.exists(path) and os.path.getsize(path) > 0:
            n = sum(1 for l in open(path, encoding="utf-8", errors="replace") if l.startswith(">"))
            total_prot += n
            print(f"[{i}/{len(entries)}] {taxid:<8} {label:<34} cached {n} proteins")
            continue
        up, kdom, got = reference_proteome(taxid)
        if not up:
            failures.append((taxid, label, "no proteome"))
            print(f"[{i}/{len(entries)}] {taxid:<8} {label:<34} NO PROTEOME", file=sys.stderr)
            continue

        # The CANONICAL one-protein-per-gene file, not the REST proteome query.
        # Querying (proteome:UP...) returns every isoform and unreviewed entry --
        # 147,506 sequences for human against ~2,000 for a bacterium, which is
        # the same unevenness this set exists to remove. The FTP file is 20,652
        # for human: differences now track real proteome size, not curation.
        import gzip
        fasta = b""
        for k in (kdom, "Eukaryota", "Bacteria", "Archaea", "Viruses"):
            try:
                blob = get(f"{FTP}/{k}/{up}/{up}_{got}.fasta.gz", tries=2)
                if blob[:2] == b"\x1f\x8b":
                    fasta = gzip.decompress(blob)
                    break
            except Exception:
                continue
        if not fasta.startswith(b">"):
            failures.append((taxid, label, f"empty fetch ({up})"))
            print(f"[{i}/{len(entries)}] {taxid:<8} {label:<34} EMPTY {up}", file=sys.stderr)
            continue
        # Key the file by the taxid whose sequences these ARE.
        path = os.path.join(out_dir, f"{got}.fasta")
        with open(path, "wb") as fh:
            fh.write(fasta)
        n = fasta.count(b"\n>") + 1
        total_prot += n
        via = "" if str(got) == str(taxid) else f" (via strain {got})"
        print(f"[{i}/{len(entries)}] {taxid:<8} {label:<34} {up} {n} proteins{via}")

    print(f"\n{len(entries)-len(failures)}/{len(entries)} taxa, {total_prot} proteins total")
    for t, l, why in failures:
        print(f"  FAILED {t} {l}: {why}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
