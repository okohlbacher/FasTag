# FasTag — A fast parallel mass spectrometry tagger and tag-based filter

Partial sequence tags from peptide MS/MS spectra, and a filter that keeps only
the spectra whose tags occur in sequences you supply. A reimplementation of the
[DirecTag](https://doi.org/10.1021/pr800154p) algorithm (Tabb et al.,
*J. Proteome Res.* 2008, 7:3838) as an [OpenMS](https://www.openms.de) TOPP tool.

**632,677 spectra from a 5.3 GB mzML in 9 seconds**, 1.2 GB peak memory — 31x
faster than the reference implementation on the same hardware, returning 5-8%
more tags.

Requires **OpenMS >= 3.5**. CI builds and runs the test suite on **Linux x64 and
arm64, and macOS x64 and arm64**, against bioconda's OpenMS 3.5.0. Also developed
against 3.6.0 (macOS arm64, AppleClang). Windows is not supported — it is the one
platform with no OpenMS conda package, so it must build OpenMS and contrib from
source; see [doc/BACKLOG-ci.md](doc/BACKLOG-ci.md).

**Input is mzML, or [mzPeak](https://github.com/OpenMS/mzpeak) when the OpenMS it
is built against provides `MzPeakFile`** — configure reports which you have.
Output is mzML and TSV. Read [the caveat](#mzpeak-input) before pointing mzPeak
at a large file.

## What it does

Tags are short sequence reads inferred straight from fragment ions — no database
search. Each is scored on peak intensity, m/z fidelity and complementarity, the
three combined by Fisher's method into an E-value. They are useful as a *fast,
sensitive prefilter*: narrow a search space, pick spectra worth a slow method, or
see what is in a run without a database.

FasTag keeps DirecTag's statistical model and changes four things.

**The intensity null is computed, not enumerated.** DirecTag walks every
C(n, k) subset of peak ranks. Measured on the reference binary:

| tag length | subsets | table build |
|---|---|---|
| 3 | 7.9e7 | 0.1 s |
| 5 | 1.7e10 | 15.1 s |
| 6 | 2.0e11 | 173 s |
| 7 | 2.1e12 | ~30 min |

The same quantity is a restricted-partition count, computable in `O(k*n*s_max)`.
FasTag builds every table for lengths 3-16 in **0.15 s**, verified equal to
exhaustive enumeration. DirecTag rebuilds that table on *every run*; this is most
of why FasTag is faster, and what makes tag lengths above four usable.

**Seed length is a minimum.** With `-extension`, each scored seed is walked
further along the spectrum graph and rescored at its realised length.

**A tag can cross one missing peak.** With `-gaps 1`, two peaks separated by a
two-residue sum bridge a hole in the ladder, and every composition matching that
sum is spelled out. Measured against Sage ground truth, this lifts the number of
spectra gaining a correctly placed tag by 44.7%.

**Tags can be restricted to sequences you supply.** `-fasta` reports only tags
occurring in the given proteins; `-out_spectra` writes the spectra carrying them,
giving a tag-filtered mzML. Matching folds I/L (isobaric, so an emitted `L` says
nothing about which residue it was), tries both orientations (tags are stored
N->C under a y-ion assumption, so b-derived tags read reversed), and accepts
isobaric substitutions derived from residue masses at the configured tolerance.
At the default `-isobaric_tolerance 0.04` (20 ppm across two peaks at m/z 1000)
the derived set is `N`=`GG`, `Q`=`GA`/`AG`, `R`=`GV`/`VG`, `W`=`AD`/`GE`/`SV`,
`K`=`GA`. No residue equals a sum of three, so one level of collapse suffices.

A minimum tag length is derived from database size, because a short tag occurs in
almost any protein by chance. Every run prints matches beside the number expected
by chance, so a count is never mistaken for signal:

```
  len         tags      matched    by chance enrichment
    4        17454          316          291       1.1x  <- at chance
    7          212          212            2      86.0x
    9          108          108            0     999.0x
```

## Speed and parallelism

Scaling is near-linear to 16 cores and useful well beyond. Output is identical at
every thread count.

**Laptop, 16 cores**, `tag_length 6`, 20 ppm:

| dataset | spectra | file | 1 thread | 16 threads | speedup |
|---|---|---|---|---|---|
| Orbitrap Astral DDA | 102,236 | 1.7 GB | 32.7 s | **3.2 s** | 10.2x |
| HeLa ddaPASEF | 357,802 | 12.2 GB | 128.7 s | **12.9 s** | 10.0x |

**Server, 128 cores** — 632,677 spectra, 5.3 GB:

| threads | 1 | 16 | 64 | 128 |
|---|---|---|---|---|
| wall | 271 s | 20 s | **9 s** | 7 s |

Throughput tracks **peaks**, not spectra — 46 M peaks/s on the Astral run — which
is what a graph over peak pairs should do.

**Memory is O(threads), not O(file)** — for the default path. Spectra stream from
the index one at a time, so the 12.2 GB ddaPASEF run peaks at 701 MB
single-threaded and 1.5 GB on 16 threads. No mzML is ever fully resident.

`-out_spectra` is the exception: it holds one slot per input spectrum in order to
write them back in input order, so that path needs memory proportional to the
file.

Against the reference implementation, same file, same hardware, run sequentially:

| | FasTag | DirecTag |
|---|---|---|
| 64 threads | **9 s** | 281 s |
| scaling, 1 to 128 threads | **32x** | 2.1x |

DirecTag stops scaling because its ranksum table build is serial and repeated
every run. FasTag's equivalent is a one-off 0.15 s.

### The fast path needs a patched OpenMS

`OnDiscMSExperiment::openFile()` parses every spectrum's metadata serially before
any work begins — ~55 s of a 60 s run, and the reason wall time otherwise stops
improving past 16 threads. A one-file OpenMS patch harvests that metadata from
the spectrum block the reader already parses; see
[doc/OPENMS-FAST-READER.md](doc/OPENMS-FAST-READER.md).

**FasTag works correctly without it**, just slower: it detects the situation at
configure time and again at runtime, warns, and takes the full-load path.

## Build

```bash
cmake -S . -B build \
  -DOpenMS_DIR=/path/to/OpenMS/build \
  -DCMAKE_PREFIX_PATH=/path/to/contrib/build
cmake --build build -j
ctest --test-dir build
```

Configure reports which reader path you get:

```
-- FasTag: using OpenMS 3.6.0 (/path/to/OpenMS/build)
-- FasTag: fast reader available -- this OpenMS reports per-spectrum metadata
```

Two portability notes, both handled automatically but worth knowing:

- **macOS**: curl is pinned to the SDK. CMake searches frameworks first, so a
  third-party `/Library/Frameworks/libcurl.framework` otherwise shadows the SDK
  libcurl that `libOpenMS.dylib` links, and the binary dies at startup with
  `Library not loaded: @rpath/libcurl.framework/...`.
- **Eigen 3.4**: its config-version file rejects the version *range*
  `OpenMSConfig.cmake` asks for (`3.4.0...<6`). Point CMake at the shim if
  configuration fails:
  `-DEigen3_DIR=$PWD/cmake/eigen3-range-shim -DEIGEN3_INCLUDE_DIR=/path/to/contrib/build/include/eigen3`

OpenMS ships a *SameMinorVersion* config-version file, so `find_package(OpenMS 3.5)`
would reject 3.6 — a minimum version cannot be expressed. FasTag asks for any
OpenMS and enforces the floor itself.

## Use

```bash
# tags for every spectrum
FasTag -in run.mzML -out tags.tsv

# high-resolution data, seed 3 extended to at most 15 residues
FasTag -in run.mzML -out tags.tsv -extension 6 -fragment_tolerance 20

# ion-trap MS2 -- set the tolerance to match the analyser
FasTag -in run.mzML -out tags.tsv -fragment_tolerance 0.3 -fragment_tolerance_unit Da

# cleaner spectra: collapse isotopes, allow one gap
FasTag -in run.mzML -out tags.tsv -deisotope -gaps 1

# only tags occurring in a protein of interest, plus the spectra carrying them
FasTag -in run.mzML -out tags.tsv -fasta AGXT.fasta -out_spectra hits.mzML

# mzPeak input, if this build has it
FasTag -in run.mzpeak -out tags.tsv
```

One row per tag: sequence, length, fragment charge, N- and C-terminal flanking
masses, whether it was extended or gapped, its E-value, and which orientation
matched.

**Set the fragment tolerance to match the analyser.** A high-resolution tolerance
on low-resolution data is silent and looks exactly like bad data — on an ion-trap
file, 20 ppm returned 3,007 tags where 0.3 Da returned 824,959. FasTag infers
resolution from peak spacing and warns when the two disagree, but it cannot know
the analyser, so the setting is yours.

### mzPeak input

Available when the OpenMS you build against provides `MzPeakFile`. Configure says
so:

```
-- FasTag: mzPeak input available (-in accepts .mzpeak)
```

Output is byte-identical to the same data as mzML, and identical at every thread
count. Reading only — `-out_spectra` writes mzML and is refused with mzPeak input.

**It is not memory-bounded, and the cause is upstream.** `MzPeakFile::transform()`
materialises the run rather than streaming it:

| input | size | peak RSS |
|---|---|---|
| .mzpeak | 93 MB | 945 MB |
| .mzpeak | **2.11 GB** | **23.7 GB** |
| same data as mzML | — | 169 MB |

~11x the file, against an mzML path that holds O(threads). FasTag's own buffer is
not the term that matters: varying it over a 16x range moved peak memory under
10%. **Prefer mzML for large runs** until the upstream reader streams —
[doc/BACKLOG-mzpeak.md](doc/BACKLOG-mzpeak.md) has the measurements, and
[doc/DESIGN-ondisc-mzpeak.md](doc/DESIGN-ondisc-mzpeak.md) proposes the fix.

## Validation

Against the reference implementation, matched parameters, 717,924 timsTOF
pseudo-DDA spectra:

- **84.5%** of spectra identified by [Sage](https://doi.org/10.1021/acs.jproteome.3c00486)
  at 1% FDR carry a tag that reads the identified peptide (DirecTag 84.2%; the
  original paper reports > 80%)
- tag sets agree **100%** where no peak-count cap applies, ~84% on dense real
  spectra, where the two tools break intensity ties differently at the cut
- tags found by both are 2-4x more likely to be confirmed than either tool's
  unique tags

Against Sage ground truth (14,867 PSMs at 1% FDR), counting spectra that gain a
correctly placed tag. Re-measured on 0.4.0, after the scoring fixes:

| | spectra | vs baseline | correct tags |
|---|---|---|---|
| default | 3,480 | — | 8,054 (7.4%) |
| `-deisotope` | 4,142 | +19.0% | 12,683 (10.6%) |
| `-gaps 1` | 5,055 | +45.3% | 33,236 (7.1%) |
| both | **5,976** | **+71.7%** | 45,919 (8.9%) |

The scoring fixes moved recall by at most 0.4% (3,479 → 3,480 at baseline,
5,035 → 5,055 with gaps), so the earlier figures were sound. They did change
per-tag correctness where gaps are on -- 6.2% → 7.1% -- because a peak can no
longer be its own complement.

Sanity-checked on four datasets spanning two vendors, three instruments, DDA and
DIA-derived input, and 163-1,435 peaks per spectrum: Bruker timsTOF (diaTracer
pseudo-MS2 and real ddaPASEF), Orbitrap Astral, Orbitrap Eclipse (ion trap).
**One acquisition per instrument**, so treat these as coverage of the input
space, not as a statistically sampled benchmark.

`ctest` covers the rank-sum DP against exhaustive enumeration, end-to-end tag
recovery from synthetic spectra, flanking-mass placement, extension, gap ordering
with a negative control, per-row mass closure, determinism, the sequence filter,
and that `-gap_penalty` reorders without ever removing a tag.

### Ranking, and the synthetic benchmark

The Sage measurement above answers "was a correct tag found". It cannot answer
"was it found **first**", because PSM truth is per spectrum, not per tag — and
that blind spot hid a real defect for months.

`bench/benchmark.cpp` closes it: synthetic spectra with exact per-tag truth,
generated from four profiles fitted to the real cornerstones' peak count and
ladder edge density (see [doc/TEST-DATA.md](doc/TEST-DATA.md)). It needs no
external data and runs in seconds.

It found that gapped tags took **95% of rank-1 slots while being 3.6× less likely
to be correct**, so enabling gaps made the single best tag worse even as it
doubled total recall. `-gap_penalty` (default 100) fixes the ordering:

| astral profile, TL=4 | rank-1 | top-5 | total recall |
|---|---|---|---|
| `-gaps 1 -gap_penalty 1` | 17.7% | 47.7% | 71.3% |
| `-gaps 1` (default 100) | **28.6%** | **52.2%** | 71.5% |

Consistent across all four profiles and tag lengths 3–5. On real ddaPASEF the
tag set is byte-identical either way — 4,993,319 tags over 287,382 spectra —
while the rank-1 gapped share falls 96.4% → 76.7%. The penalty **reorders and
never filters**, which is what keeps a sensitive prefilter sensitive.

It also found an uncaught exception: `-deisotope` at an ion-trap tolerance
crashed the tool, a combination reachable straight from the docs.

## Known limitations

- **`-gap_penalty`'s default is calibrated, not derived.** It rests on synthetic
  spectra this project wrote to test itself. The multiplicity argument justifies
  only ~10 of the 100; measurement says the real over-scoring is larger but gives
  no clean optimum. Within the gapped family the E-value still carries almost no
  information about correctness. Real ground truth would settle it —
  [doc/BACKLOG.md](doc/BACKLOG.md).
- **The Sage ground truth is no longer reproducible here.** The PSM table, the
  S23 file and the Sage binary are all gone from the machine those numbers were
  measured on. They stand as recorded; they cannot currently be extended.
- **Defaults are conservative.** `-deisotope`, `-gaps` and `-peaks_per_window`
  are all off. Their measured gains come from one acquisition type, and the peak
  cap in particular was tuned on data whose peak count is clamped flat at 500.
- **No modification support.** Residues are the unmodified 19, so labelled
  samples (TMT and similar) will not match tags spanning a modified residue.
- **mzPeak input is not memory-bounded** and is not covered by CI, which builds
  against a stock OpenMS that lacks the file handler. mzML is the tested path.

## Licence and provenance

MIT — see [LICENSE](LICENSE). Dependencies and licences are in [BOM.md](BOM.md),
with what FasTag deliberately does *not* implement because OpenMS or Boost
already provides it.

The DirecTag algorithm is reimplemented from its publication. The reference
implementation (Apache-2.0) was read to resolve semantics the paper leaves
implicit; **no code was copied**.
