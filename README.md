# FasTag

Partial sequence tags from peptide MS/MS spectra — a reimplementation of the
[DirecTag](https://doi.org/10.1021/pr800154p) algorithm (Tabb et al.,
*J. Proteome Res.* 2008, 7:3838) as an [OpenMS](https://www.openms.de) TOPP tool.

Requires **OpenMS 3.5** or newer. Built and verified against OpenMS 3.5.0
(Linux/gcc-13) and 3.6.0 (macOS arm64/AppleClang).

## Why

DirecTag infers sequence tags directly from fragment ions, scoring each on peak
intensity, m/z fidelity and complementarity, and combining the three by Fisher's
method. FasTag keeps that statistical model and changes three things.

**The rank-sum null is computed by dynamic programming.** DirecTag's
`CalculateIRBins_R` enumerates every C(n, k) subset of peak ranks to build the
intensity null. Measured on the reference binary:

| tag length | subsets enumerated | table build |
|---|---|---|
| 3 | 7.9 × 10⁷ | 0.1 s |
| 5 | 1.7 × 10¹⁰ | 15.1 s |
| 6 | 2.0 × 10¹¹ | 173 s |
| 7 | 2.1 × 10¹² | ~30 min |

The same distribution is a restricted-partition count, computable in
`O(k·n·s_max)`. FasTag builds every table for lengths 3–16 in **0.15 s**, verified
equal to exhaustive enumeration. This is what makes tag lengths above 4 reachable.

**Seed length is a minimum, not a fixed value.** With `-extension`, each scored
seed is greedily extended along the spectrum graph and rescored at its realised
length. On high-resolution data this roughly triples mean tag length for ~1.9×
the cost and ~1 point of accuracy; at 0.5 Da it costs 11 points, so it is off by
default and worth enabling only above ~0.05 Da resolution.

**Tags can be restricted to sequences you supply.** `-fasta` reports only tags
occurring in the given proteins, and `-out_spectra` writes the spectra that carry
them — a tag-filtered spectrum file. Matching folds I/L (isobaric, so an emitted
`L` carries no information about which it was), tries both orientations (tags are
stored N→C under a y-ion assumption, so b-derived tags read reversed), and
optionally accepts isobaric residue/pair substitutions (`N` ≡ `GG`, `Q` ≡ `GA`,
and at wider tolerance `R` ≡ `GV`, `W` ≡ `AD`/`GE`/`SV`, `K` ≡ `GA`).

A minimum tag length is derived from the database size, because a short tag
occurs in almost any protein by chance: against the human proteome the filter is
the identity function below length 6. The tool prints the derived floor and, per
length, how many tags matched beside how many were expected by chance.

## Build

```bash
cmake -S . -B build \
  -DOpenMS_DIR=/path/to/OpenMS/build \
  -DCMAKE_PREFIX_PATH=/path/to/OpenMS/contrib/build
cmake --build build -j
ctest --test-dir build
```

`CMAKE_PREFIX_PATH` must point at the same contrib tree OpenMS was built against,
so `OpenMSConfig.cmake` can resolve Xerces-C, Boost, LIBSVM and friends.

On **macOS** the build resolves curl to the SDK automatically. CMake searches
frameworks first, so a third-party `/Library/Frameworks/libcurl.framework` would
otherwise shadow the SDK libcurl that `libOpenMS.dylib` links, and the binary
would die at startup with `Library not loaded: @rpath/libcurl.framework/...`.
Override with `-DCURL_LIBRARY=...` if you need a different curl.

If configuration fails with *"Could not find a configuration file for package
Eigen3 ... compatible with requested version range 3.4.0...<6"*, you have hit a
bug in Eigen 3.4's own config-version file: it enforces a same-major policy that
rejects the range `OpenMSConfig.cmake` requests, so a perfectly good Eigen 3.4.0
is refused. Point CMake at the shim in this repo:

```bash
  -DEigen3_DIR=$PWD/cmake/eigen3-range-shim \
  -DEIGEN3_INCLUDE_DIR=/path/to/contrib/build/include/eigen3
```

The algorithm core (`include/FasTag/`) is header-only and has no dependencies, so
the tests build and run without OpenMS:

```bash
cmake -S . -B build -DFASTAG_BUILD_TOOL=OFF
cmake --build build -j && ctest --test-dir build
```

## Use

```bash
# tags for every spectrum
FasTag -in run.mzML -out tags.tsv

# seed 3, extended to at most 15 residues, high-resolution data
FasTag -in run.mzML -out tags.tsv -extension 6 -fragment_tolerance 20 -fragment_tolerance_unit ppm

# only tags occurring in a protein of interest, plus the spectra carrying them
FasTag -in run.mzML -out tags.tsv -fasta AGXT.fasta -out_spectra hits.mzML
```

Output is one row per tag: sequence, length, fragment charge, N- and C-terminal
flanking masses, whether it was extended, its E-value, and (when filtering) which
orientation matched.

## Status

The algorithm core is complete and tested: rank-sum DP verified against exhaustive
enumeration, an end-to-end synthetic recovery test, extension, and the sequence
filter (encoding, soundness, completeness, specificity, orientation, degenerate
FASTA input). `ctest` runs all of it without OpenMS.

The TOPP tool builds and runs against OpenMS 3.5.0. It follows TOPP conventions:
`-write_ini` / `-ini` round-trip (21 parameters), standard logging, citations, and
`--help`. It is registered as a non-official tool (`official = false` in the
`TOPPBase` constructor), because it lives outside the OpenMS tree and is therefore
not in `ToolHandler`'s list; everything else behaves as any TOPP tool does.

GitHub CI builds and tests the algorithm core only, since it has no dependencies.
The TOPP tool is not built in CI — that would need an OpenMS installation.

Validation against the reference implementation, on 717,924 timsTOF pseudo-DDA
spectra with matched parameters:

- 7,310 spectra tagged in 2.1 s (97,913 tags); with `-fasta` restricting to a
  single 392-residue protein, 601 tags over 330 spectra, written to a filtered
  mzML that reads back cleanly
- **84.5%** of spectra identified by [Sage](https://doi.org/10.1021/acs.jproteome.3c00486)
  at 1% FDR carry at least one tag that reads the identified peptide
  (DirecTag: 84.2%; the original paper reports >80%)
- tag sets agree **100%** on spectra where no peak-count cap applies, and ~84% on
  dense real spectra, where the two tools break intensity ties differently at the
  100-peak boundary
- where they disagree, tags found by both are 2–4× more likely to be confirmed
  than either tool's unique tags

## Licence

BSD-3-Clause, matching OpenMS. The DirecTag algorithm is reimplemented from its
publication and from reading the reference implementation; no code was copied.
The reference (`freicore`/DirecTag, Vanderbilt University) is Apache-2.0.
