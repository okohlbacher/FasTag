# Species / taxonomic detection

FasTag can infer the taxa present in a run directly from its tags — no database
search — by matching tag k-mers against a prebuilt tag→taxon index and ranking
taxa by lowest-common-ancestor evidence. Kraken-2-for-peptide-tags, over a
**reduced** reference set so the resolution is honestly genus/family.

## Use

```bash
# 1. Build the index once, from a directory of <taxid>.fasta files
buildtaxdb proteomes/ ref.taxdb 7        # k=7

# 2. Detect, alongside normal tagging (subsample for speed on big runs)
FasTag -in run.mzML -out tags.tsv -tag_length 7 \
       -taxdb ref.taxdb -taxonomy_nodes nodes.dmp -taxonomy_names names.dmp \
       -species_out taxa.tsv -species_rank genus -subsample_fraction 0.2
```

`species_out` columns: rank, taxid, name, observed (spectra supporting the
taxon), expected (from the taxon's index breadth), enrichment, log_pvalue,
qvalue. Tags shorter than the index `k` are ignored (`-species_min_len`).

## How

- The index maps each I/L-folded k-mer to the SET of taxa whose proteins carry
  it (breadth, not abundance). A **reduced** reference — one representative
  proteome per genus/family — is deliberate: short peptide k-mers cannot resolve
  species, and a curated set keeps the LCA honest and the index small.
- A tag supports a taxon only if EVERY one of its k-mers is in that taxon
  (intersection), so a longer tag is more specific. Each taxon is counted ONCE
  per spectrum, not per tag, so correlated tags cannot manufacture evidence.
- Per-leaf counts and the index breadth are both rolled up the NCBI tree, and
  each node is tested against its subtree breadth (binomial + Benjamini-Hochberg
  in `TaxStats`). Reported per `-species_rank`.

## Validated end-to-end

Two real runs, a reduced 17-taxon reference (target + Enterobacterales
near-neighbours + background bacteria + the spike/contaminant animals):

- **PXD000001** (*Erwinia/Pectobacterium carotovora*, the ECA=*P. atrosepticum*
  SCRI1043 proteome in the index): **Pectobacterium ranks #1** by evidence
  (1807 spectra), and the spikes/contaminants surface — **Bos** (BSA +
  cytochrome C), **Oryctolagus** (rabbit PYGM), **Sus** (trypsin).
- **Human CPTAC (Thermo Lumos)**: **Homo ranks #1** (3326), no bacterial genus
  near the top — the detector identifies the sample organism and does not
  spuriously call the bacterial target. A clean discrimination control.

## Honest limitations

- **Conserved proteins inflate near-neighbours.** Bacterial genera (Escherichia,
  Salmonella, Yersinia) appear on the Erwinia run because conserved bacterial
  proteins share k-mers; the true genus still wins on evidence, but the tail is
  not clean. Rank by `observed`, not `enrichment` (small-DB genera get spuriously
  high enrichment: few proteins -> tiny expected).
- **The q-value is a ranking aid, not a calibrated FDR.** The per-k-mer
  background is a proxy and the chimeric-null problem (F4 in `BACKLOG.md`)
  applies. Trust the top-ranked evidence, not the absolute q.
- **Reduced reference = genus/family only.** Species calls are out of scope by
  construction, and a taxon absent from the reference cannot be called (it will
  present as its nearest represented relative).
- Contaminants (trypsin=Sus, serum=Bos, keratin=Homo) are real signal, not
  noise — they show up because they are in the sample.
