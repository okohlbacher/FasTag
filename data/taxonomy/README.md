# Bundled taxonomy for `-species`

FasTag looks for this directory at `<bin>/../share/FasTag/taxonomy` (a normal
install) or `<bin>/share-FasTag-taxonomy` (the flat release tarball), or wherever
`FASTAG_TAXONOMY_DIR` points. With the index in place, `-species` needs no paths:

```
FasTag -in run.mzML -out run.tags.tsv -species -tag_length 7
```

`-tag_length` must be **at least the index k** (7 for the shipped index) — a
shorter tag cannot be looked up, and FasTag refuses the run up front rather than
producing an empty report.

## What ships where

| file | size | in the release tarball? |
|---|---|---|
| `nodes.dmp` | ~10 KB | yes |
| `names.dmp` | ~19 KB | yes |
| `tax_k7.taxdb` | ~1.06 GB | **no — separate release asset `FasTag-taxonomy-k7.tar.gz`** |

The two `.dmp` files are an NCBI taxdump **pruned to the lineages the index can
reach** (a few hundred taxa instead of ~2.5M). Produced by
`tools/prune_taxdump.py`; identical results to the full ~500 MB dump, and faster
to load.

The **index** is the reduced-set k=7 image over the reference proteomes listed in
`reference-set.tsv` (~50 taxa: model organisms, MS contaminants, common
pathogens, gut microbiome, archaea). It is not in the per-platform tarballs
because it is large and platform-independent. Get it once:

```
# download FasTag-taxonomy-k7.tar.gz from the release, then:
tar xzf FasTag-taxonomy-k7.tar.gz -C "$(dirname "$(command -v FasTag)")/.."   # into share-FasTag-taxonomy/
```

Without it, `-species` fails with a message naming the directory it searched and
pointing at this asset — never silently.

## The index format (v2, `FTX2`)

Memory-mapped and bit-packed, so a ~1 GB index costs **under 1 GB resident** and
loads instantly (nothing is deserialized):

- k-mers packed base-19 (the I/L-folded alphabet is 19 letters) → 4-byte keys;
- postings store a 1-byte taxon **index**, not a 4-byte taxid;
- the file is `mmap`'d and used in place, so its pages are file-backed and
  evictable under memory pressure.

The earlier `FTXI` format (a serialized hash map, ~10× the RAM) still loads for
backward compatibility. `save()` always writes `FTX2`.

## Rebuilding the index

```
python3 tools/fetch_reference_set.py data/taxonomy/reference-set.tsv proteomes
buildtaxdb proteomes tax_k7.taxdb 7                # k = 7; max is 12
python3 tools/prune_taxdump.py <ncbi>/nodes.dmp <ncbi>/names.dmp out/ --from-proteomes proteomes
```

The build is deterministic: the same proteomes give a byte-identical index, so a
release can reproduce it and a checksum verifies a download.

## Scope, honestly

A **reduced reference set**, not a metagenomics database. Reference proteomes
(one protein per gene) keep the taxa evenly weighted — the earlier reviewed-only
set had 20,431 human proteins against 35 for a bacterium, which skewed the
background model. Even so, short peptide 7-mers resolve **genus/family at best**:
on a human sample *Homo* ranks first but *Macaca*, *Bos*, *Sus* follow closely,
because primate and mammalian proteomes share most 7-mers. Read the ranking as
"which lineages are present", not "which species", and treat q as a ranking aid,
not a calibrated FDR.
