# mzPeak: what works, what does not, what is next

Reading **and writing** `.mzpeak` work. Reading is not memory-bounded, and that
is upstream rather than here.

## What shipped

`-in run.mzpeak`, when FasTag is built against an OpenMS providing `MzPeakFile`.
Detected at configure time, so the same source builds either way:

```
-- FasTag: mzPeak input available (-in accepts .mzpeak)
```

Writing too: `-out_spectra hits.mzpeak`. All four in/out combinations work, and
the earlier blanket refusal of `-out_spectra` with mzPeak input is gone --
that path only ever needed the run-level `ExperimentalSettings` wired in from
the streaming consumer, which the mzML paths get from `getMetaData()` or the
loaded map.

Writing does NOT go through `FileHandler::storeExperiment()`: it has no mzPeak
branch, so asking it for one silently writes a different format. `MzPeakFile`
is called directly.

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
| mzML -> mzpeak -> tags | **exact round-trip**: 127,035 tags, identical to tagging the source mzML |
| mzpeak -> tags | 42,092 MS2 / 883,939 tags on a 155 MB Lumos file -- same spectrum count as its mzML |
| mzpeak -> mzpeak -> tags | 883,939 tags, unchanged through the write and re-read |

### Timing: mzPeak vs mzML, same data, 16 cores

Same acquisition in both formats -- Thermo Lumos, 42,092 MS2 spectra, 589 MB
as mzML and 155 MB as mzpeak. Apple M-series, 16 cores, 128 GB; page cache
warmed for both files first; two reps per cell, spread under 2%.

| threads | mzML wall | mzpeak wall | mzpeak speedup |
|---|---|---|---|
| 1 | 15.07 s | 6.06 s | 2.5x |
| 2 | 9.23 s | 3.82 s | 2.4x |
| 4 | 6.21 s | 2.63 s | 2.4x |
| 8 | 4.54 s | 1.95 s | 2.3x |
| 16 | **3.83 s** | **1.69 s** | **2.3x** |

Peak RSS tells the other half of the story, and the two shapes are the point:

| threads | mzML | mzpeak |
|---|---|---|
| 1 | 339 MB | 1686 MB |
| 16 | 444 MB | 1686 MB |

**mzPeak buys ~2.3x wall time with ~4x memory.** mzML grows slowly with thread
count (O(threads), one streamed spectrum in flight per worker); mzPeak is
completely flat because `transform()` materialises the whole run up front --
the same non-streaming behaviour documented above, seen from the cost side.

Why it is faster is not mysterious: the mzML path decodes XML, base64 and zlib
per spectrum *inside* the parallel loop, while mzPeak has already paid a
vectorised columnar Parquet decode into RAM and the loop does nothing but tag.
Fitting Amdahl to the two curves puts the serial part at ~3.1 s for mzML
against ~1.4 s for mzpeak, and the parallel part at ~12.0 s against ~4.7 s.

**Two caveats before quoting any of this.**

The OpenMS here does NOT have the fast-reader patch (`doc/OPENMS-FAST-READER.md`),
so the mzML path pays a serial up-front metadata parse that a patched build
skips -- roughly the ~1.7 s difference in the fitted serial terms. A patched
OpenMS would narrow this gap, and nothing here measures by how much.

Tag counts are not identical: 905,760 from mzML against 883,939 from mzpeak,
97.6%. The mzpeak file is a conversion of the same raw data, and its float32
m/z storage against mzML's float64 moves a small number of borderline tags
across the tolerance. Not a correctness difference between the readers -- a
FasTag-written mzpeak round-trips to exactly the source file's tag set, which
is the controlled comparison.

### Two upstream bugs found while wiring this up

Both silent, both fixed in [OpenMS-mzPeakRW PR #1](https://github.com/okohlbacher/OpenMS-mzPeakRW/pull/1),
and either one alone makes mzPeak unusable from an *installed* OpenMS.

1. **Every `dynamic_pointer_cast<arrow::*Array>` returns null** when Arrow is
   built with hidden symbol visibility -- the conda-forge default. `libOpenMS`
   and `libarrow` hold distinct typeinfo, the cross-DSO `dynamic_cast` cannot
   match, and the reader silently drops precursors, CV params and the profile
   `mz_delta_model`. Measured: **0 of 42,092** precursors survived a load, while
   the Parquet held all 42,092. FasTag needs a precursor to tag, so every
   spectrum was rejected and the run reported a confident zero.

   Diagnosed by replicating `readPrecursors_`'s exact logic in a standalone
   Arrow program linked against the *same* libarrow, where it found all 42,092
   rows -- so the difference is purely the shared-library boundary, not the data.

2. **`MzPeakFile.h` was in no CMake list**, so it never installed. Invisible
   when building against an OpenMS *build tree* (which exposes source includes
   directly), fatal for a consumer of an install: `FileTypes.h` ships with
   `MZPEAK` but the class is absent, so FasTag's header-based feature detection
   silently produced a build with no mzPeak support.

Worth remembering as a class: both failures were *green builds that shipped
without the feature*, which is why CI now asserts `MzPeakFile.h` is present in
the installed OpenMS rather than trusting the build to have noticed.

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
