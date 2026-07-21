# mzPeak: what works, what does not, what is next

Reading `.mzpeak` works and produces byte-identical output to the same data as
mzML. It is **not memory-bounded**, and that is upstream rather than here.

## What shipped

`-in run.mzpeak`, when FasTag is built against an OpenMS providing `MzPeakFile`.
Detected at configure time, so the same source builds either way:

```
-- FasTag: mzPeak input available (-in accepts .mzpeak)
```

Reading only. `-out_spectra` still writes mzML, and is refused with mzPeak input.

### How

`MzPeakFile` is push-based (`transform()` calls `consumeSpectrum()` once per
spectrum); FasTag's loop wants a batch to parallelise over. `ChunkingConsumer`
buffers into chunks and hands each to the same per-spectrum work the mzML path
uses — one shared `tag_one`, so the two readers cannot drift.

Two decisions worth keeping:

- **The chunk is bounded by peaks, not spectra.** Spectra in this corpus run from
  ~100 peaks to over 130,000, so "2048 spectra" is anywhere between 200 K and
  270 M peaks. A peak budget is flat whatever the data looks like.
- **MS1 and precursor-less spectra are dropped in the consumer**, not in the
  tagging callback, so they never occupy the buffer. On DDA that is most of the
  input.

### Verified

| | |
|---|---|
| output vs mzML | **byte-identical**, 847,528 tags on the Eclipse DDA file |
| determinism | identical at 1, 4 and 16 threads |
| stock OpenMS | still builds; `.mzpeak` refused with an actionable message |

## The problem: transform() is not streaming

| input | size | peak RSS |
|---|---|---|
| Eclipse DDA, mzpeak | 93 MB | 945 MB |
| **HeLa diaPASEF, mzpeak** | **2.11 GB** | **23,735 MB** |
| same data, mzML | — | 169 MB |

**~11x the file size**, against an mzML path that holds O(threads).

The buffer is not the cause. Varying the chunk budget over a 16x range moved peak
memory under 10%:

| budget | peak RSS |
|---|---|
| 4 M peaks | 1027 MB |
| 1 M peaks | 945 MB |
| 250 K peaks | 986 MB |

If the buffer dominated, memory would track it. It does not, so the term that
matters is inside `MzPeakFile::transform()`, which materialises the run despite
the consumer interface existing precisely to avoid that.

**Guidance: mzPeak is fine for files small relative to RAM. Use mzML for large
runs.**

## Next

1. **Make `MzPeakFile::transform()` stream.** An OpenMS-side fix: read row group
   by row group and emit as it goes, rather than assembling first. This is the
   whole of the problem and everything else is a workaround.
2. **`OnDiscMzPeakExperiment`** — random access over Parquet row groups, matching
   `OnDiscMSExperiment`. Would let the mzPeak path use the *same* pull-based loop
   as mzML instead of a separate chunked one, deleting the consumer entirely.
   The right long-term shape, and an OpenMS contribution rather than a FasTag
   change.
3. **`-out_spectra` for mzPeak**, once writing is worth having. Blocked on
   `MzPeakFile::store()`, whose own documentation says "Run-level metadata and
   precursor facets are not yet emitted" — a converted file would lose the
   precursor, which is exactly what tagging needs.

## Caveat

mzPeak is pre-1.0: *"no stability is guaranteed at this point"*. The API surface
used here is deliberately tiny — one `transform()` call and one consumer — so
churn lands in one class and one dispatch site.
