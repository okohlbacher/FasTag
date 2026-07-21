# Design: `OnDiscMzPeakExperiment` — streaming random access for mzPeak

A proposal for OpenMS, mirroring `OnDiscMSExperiment` for mzML. Written from
FasTag's experience of consuming both.

**One-line summary: mzPeak should be *easier* to stream than mzML, not harder,
and the current interface gives that advantage away.**

## The problem, measured

FasTag reads mzML with `OnDiscMSExperiment`: pull-based random access,
`getSpectrum(i)`, one reader per thread. Memory is O(threads):

| input | size | peak RSS |
|---|---|---|
| mzML | 12.2 GB | 1.5 GB |
| mzML | 5.3 GB | 169 MB |

The same tool reading mzPeak through `MzPeakFile::transform()` and a chunking
consumer that buffers a bounded 1 M peaks:

| input | size | peak RSS |
|---|---|---|
| .mzpeak | 93 MB | 945 MB |
| **.mzpeak** | **2.11 GB** | **23,735 MB** |

**~11x the file.** The consumer is not the cause: varying its budget over a 16x
range moved peak memory under 10% (1027 / 945 / 986 MB), which is what shows the
buffer is not the term that matters. `transform()` materialises the run.

## Why pull, not push

Push is not wrong in itself, but it forces two things:

1. **The consumer cannot control what is resident.** It is handed spectra; it
   cannot say "give me spectrum 40,000 now and nothing else". Any bound it
   imposes is downstream of whatever the reader already assembled.
2. **The caller must invent its own parallelism.** FasTag buffers into chunks,
   runs OpenMP over each, then blocks while the next chunk refills
   single-threaded. Decode and compute serialise against each other. The mzML
   path has no such barrier because each thread pulls its own spectrum.

Both dissolve with a pull interface.

## Proposed interface

Deliberately the same shape as `OnDiscMSExperiment`, so a consumer can be written
once against both:

```cpp
class OPENMS_DLLAPI OnDiscMzPeakExperiment
{
public:
  OnDiscMzPeakExperiment() = default;
  OnDiscMzPeakExperiment(const OnDiscMzPeakExperiment&);   // per-thread copies
  // operator= deliberately absent, as in OnDiscMSExperiment

  /// Read the Parquet footers and the metadata table. Does NOT read peak data.
  bool openFile(const String& filename, bool skip_meta_data = false);

  Size getNrSpectra() const;
  Size getNrChromatograms() const;

  /// Decode one spectrum. NOT thread-safe -- give each thread its own copy.
  MSSpectrum getSpectrum(Size id);
  MSChromatogram getChromatogram(Size id);

  boost::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const;

  // --- the part mzML cannot offer -------------------------------------------

  /// Per-spectrum metadata WITHOUT decoding peaks: MS level, RT, precursor m/z,
  /// precursor charge. Backed by spectra_metadata.parquet, so this is a columnar
  /// read of a few narrow columns.
  const SpectrumMetaData& getMetaData(Size id) const;

  /// Row-group boundaries, for callers that want to parallelise over them.
  /// Returns [first_spectrum, last_spectrum) per group.
  const std::vector<std::pair<Size, Size>>& rowGroups() const;
};
```

`getSpectrum` and the copy-constructor rule are identical to
`OnDiscMSExperiment`, so its consumers port by changing a type name.

## Three things Parquet makes better than mzML

**1. The metadata pass becomes nearly free.** This is the big one.

`OnDiscMSExperiment::openFile()` calls `loadMetaData_()`, which parses every
spectrum's metadata before any work starts: measured at ~55 s of a 60 s run on a
5.3 GB mzML, single-threaded, and the reason wall time stopped improving past 16
threads however many cores it was given. We patched `MzMLSpectrumDecoder` to
harvest metadata from the DOM it already builds, which took the same run from
60 s to 9 s (see `doc/OPENMS-FAST-READER.md`).

mzPeak needs none of that. `spectra_metadata.parquet` is one row per spectrum,
already columnar. Reading `ms_level`, `rt`, `precursor_mz`, `precursor_charge`
touches four narrow columns and never decodes a peak. **What cost 55 s in mzML
should cost milliseconds here** — and it is the format's natural layout, not a
workaround bolted onto it.

**2. Filtering before I/O.** With MS level available up front, a caller wanting
only MS2 can skip row groups containing none. On a DDA run that is most of the
MS1 data never read at all. In mzML you must parse a spectrum to learn its level.

**3. Column pruning.** A caller needing only m/z and intensity should not pay for
ion-mobility or other arrays present in the file. Parquet does this natively;
mzML cannot.

## Parallelise over row groups, not spectra

The obvious port of FasTag's loop is `#pragma omp for` over `getSpectrum(i)`. For
Parquet that is wrong: one row group holds many spectra, so per-spectrum random
access re-reads and re-decodes the same group repeatedly, and two threads on
neighbouring spectra contend for it.

Give the caller `rowGroups()` and let it parallelise over those instead:

```cpp
const auto& groups = exp.rowGroups();
#pragma omp parallel
{
  OnDiscMzPeakExperiment reader(exp);          // per-thread copy
#pragma omp for schedule(dynamic, 1)
  for (SignedSize g = 0; g < (SignedSize)groups.size(); ++g)
    for (Size i = groups[g].first; i < groups[g].second; ++i)
      rows[i] = work(reader.getSpectrum(i));   // indexed -> order preserved
}
```

Each thread touches distinct groups, so each is decoded exactly once, there is no
cross-thread contention, and output stays indexed by global spectrum id so order
is independent of scheduling. Internally `getSpectrum` then needs only a
**one-entry per-reader cache** of the currently decoded group; sequential access
within a group hits it every time.

This is strictly better than what mzML can do, where decode granularity is one
spectrum and there is nothing to amortise.

## Memory

Resident = (threads x one decoded row group) + (threads x one spectrum) + caller
state. At ~1 M rows per group and ~20 bytes per row of point-layout peaks that is
~20 MB per thread, so **~320 MB at 16 threads independent of file size**, and
tunable by choosing the row group size at write time — a knob the format already
has.

Against 23.7 GB today.

## What FasTag would delete

`ChunkingConsumer` and its dispatch, ~70 lines, and with them the chunk-refill
barrier. Both input paths would then run the same pull loop over the same shared
per-spectrum function — the arrangement the mzML path already has, and the reason
its readers cannot drift apart.

## Risks and unknowns

- **I have not read mzpeak's internal layout.** This is designed from the OpenMS
  `MzPeakFile` header, the `store()` documentation naming
  `spectra_data.parquet` / `spectra_peaks.parquet` / `spectra_metadata.parquet` /
  `mzpeak_index.json`, and the upstream header *names* `index.h`, `spectra.h`,
  `data/array_index.h`, `util/slice.h` — which suggest a slice/range abstraction
  already exists. If it does, most of this is plumbing rather than new work.
  **Check this first, before anyone commits to the design.**
- **Spectra spanning row groups.** A large spectrum may cross a boundary, so
  `getSpectrum` must stitch two groups. Bounded, but it is the fiddly case.
- **Random access across the ZIP container.** Entries are STORED (uncompressed)
  per the `store()` docs, so byte offsets should be directly seekable. Worth
  confirming rather than assuming.
- **mzPeak is pre-1.0** — "no stability is guaranteed at this point". This should
  land alongside the format stabilising, not ahead of it.
- **Chromatograms** appear only for interface symmetry; FasTag does not use them.

## Smaller fix worth doing first

If `OnDiscMzPeakExperiment` is too large a change for now, make `transform()`
actually stream: read row group by row group and emit as it goes, rather than
assembling first. That alone would take 23.7 GB to something bounded, and every
existing consumer benefits with no API change.

Note `transform()` is documented as two-pass, with `skip_first_pass` to suppress
the first. Two passes explain time, not an 11x memory factor, so the
materialisation is inside a pass rather than a consequence of there being two.
