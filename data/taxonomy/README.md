# Bundled taxonomy for `-species`

FasTag looks for this directory at `<bin>/../share/FasTag/taxonomy`, or wherever
`FASTAG_TAXONOMY_DIR` points. With it in place, `-species` needs no paths:

```
FasTag -in run.mzML -out run.tags.tsv -species -tag_length 7
```

## What ships here

| file | size | in the release tarball? |
|---|---|---|
| `nodes.dmp` | 3 KB | yes |
| `names.dmp` | 5 KB | yes |
| `tax_k7.taxdb` | 232 MB (88 MB compressed) | **no — separate release asset** |

The two `.dmp` files are an NCBI taxdump **pruned to the lineages the index can
actually reach** (117 taxa instead of ~2.5M). Produced by
`tools/prune_taxdump.py`; output is byte-identical to running against the full
500 MB dump, and loads faster because there is nothing else to parse.

The index itself is not in the per-platform tarballs on purpose: it is
platform-independent, so copying it into all five would add ~440 MB per release
of identical bytes. It is published once per release as
`FasTag-taxonomy.tar.gz`. Install it into this directory:

```
tar xzf FasTag-taxonomy.tar.gz -C "$(dirname "$(command -v FasTag)")/../share/FasTag/taxonomy"
```

Without it, `-species` fails with a message naming the directory it searched —
never silently.

## Rebuilding the index

```
buildtaxdb <proteomes_dir> tax_k7.taxdb 7      # <taxid>.fasta per organism
python3 tools/prune_taxdump.py nodes.dmp names.dmp out/ --from-proteomes <proteomes_dir>
```

## Scope, honestly

17 representative proteomes (Homo, Bos, Sus, Oryctolagus, Saccharomyces + 12
bacteria). It is a **demo-grade reduced set**: good for genus/family-level
discrimination and for spotting the usual contaminants, not a reference for
claiming what is in an unknown sample. Coverage is also uneven — 20,431 human
proteins against 35 for *S. pneumoniae* — which skews the background model, so
compare enrichment across taxa with that in mind.
