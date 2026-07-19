# FasTag

Partial sequence tags from peptide MS/MS spectra — a reimplementation of the
[DirecTag](https://doi.org/10.1021/pr800154p) algorithm (Tabb et al.,
*J. Proteome Res.* 2008, 7:3838) as an [OpenMS](https://www.openms.de) TOPP tool.

Requires **OpenMS ≥ 3.5**. Built and verified against 3.5.0 (Linux, gcc-13) and
3.6.0 (macOS arm64, AppleClang).

## Why

DirecTag infers sequence tags directly from fragment ions, scoring each on peak
intensity, m/z fidelity and complementarity and combining the three by Fisher's
method. FasTag keeps that model and changes three things.

**The intensity null is computed, not enumerated.** DirecTag's
`CalculateIRBins_R` walks every C(n, k) subset of peak ranks. Measured on the
reference binary:

| tag length | subsets | table build |
|---|---|---|
| 3 | 7.9 × 10⁷ | 0.1 s |
| 5 | 1.7 × 10¹⁰ | 15.1 s |
| 6 | 2.0 × 10¹¹ | 173 s |
| 7 | 2.1 × 10¹² | ~30 min |

The same quantity is a restricted-partition count, computable in `O(k·n·s_max)`.
FasTag builds every table for lengths 3–16 in **0.15 s**, verified equal to
exhaustive enumeration. This is what makes tag lengths above four usable.

**Seed length is a minimum.** With `-extension`, each scored seed is walked
further along the spectrum graph and rescored at its realised length. At high
resolution this roughly triples mean tag length for ~1.9× the cost and about one
point of accuracy; at 0.5 Da it costs eleven points, so it is off by default.

**Tags can be restricted to sequences you supply.** `-fasta` reports only tags
occurring in the given proteins and `-out_spectra` writes the spectra carrying
them — a tag-filtered spectrum file. Matching folds I/L (isobaric, so an emitted
`L` says nothing about which residue it was), tries both orientations (tags are
stored N→C under a y-ion assumption, so b-derived tags read reversed), and
accepts isobaric substitutions derived from the residue masses at the configured
tolerance (`N` ≡ `GG`, `Q` ≡ `GA`, and at wider tolerance `R` ≡ `GV`,
`W` ≡ `AD`/`GE`/`SV`, `K` ≡ `GA`).

A minimum tag length is derived from database size, because a short tag occurs in
almost any protein by chance — against a whole proteome the filter is the
identity function below length 6. Every run prints matches beside the number
expected by chance, so a count is never mistaken for signal:

```
  len         tags      matched    by chance enrichment
    4        17454          316          291       1.1x  <- at chance
    7          212          212            2      86.0x
    9          108          108            0     999.0x
```

## Build

```bash
cmake -S . -B build \
  -DOpenMS_DIR=/path/to/OpenMS/build \
  -DCMAKE_PREFIX_PATH=/path/to/contrib/build
cmake --build build -j
ctest --test-dir build
```

`CMAKE_PREFIX_PATH` must point at the same contrib tree OpenMS was built against,
so `OpenMSConfig.cmake` can resolve Xerces-C, Boost, LIBSVM and the rest.

Two portability notes, both handled automatically but worth knowing:

- **macOS**: curl is pinned to the SDK. CMake searches frameworks first, so a
  third-party `/Library/Frameworks/libcurl.framework` otherwise shadows the SDK
  libcurl that `libOpenMS.dylib` links, and the binary dies at startup with
  `Library not loaded: @rpath/libcurl.framework/...`.
- **Eigen 3.4**: its config-version file rejects the version *range*
  `OpenMSConfig.cmake` asks for (`3.4.0...<6`), so configuration can fail against
  a perfectly good Eigen. Point CMake at the shim if it does:
  `-DEigen3_DIR=$PWD/cmake/eigen3-range-shim -DEIGEN3_INCLUDE_DIR=/path/to/contrib/build/include/eigen3`

OpenMS ships a *SameMinorVersion* config-version file, so `find_package(OpenMS 3.5)`
would reject 3.6 — a minimum version cannot be expressed. FasTag asks for any
OpenMS and enforces the floor itself.

## Use

```bash
# tags for every spectrum
FasTag -in run.mzML -out tags.tsv

# seed 3 extended to at most 15 residues, high-resolution data
FasTag -in run.mzML -out tags.tsv -extension 6 -fragment_tolerance 20

# only tags occurring in a protein of interest, plus the spectra carrying them
FasTag -in run.mzML -out tags.tsv -fasta AGXT.fasta -out_spectra hits.mzML
```

One row per tag: sequence, length, fragment charge, N- and C-terminal flanking
masses, whether it was extended, its E-value, and which orientation matched.

## Validation

On 717,924 timsTOF pseudo-DDA spectra, with matched parameters against the
reference implementation:

- **84.5%** of spectra identified by [Sage](https://doi.org/10.1021/acs.jproteome.3c00486)
  at 1% FDR carry at least one tag that reads the identified peptide
  (DirecTag: 84.2%; the original paper reports >80%)
- tag sets agree **100%** where no peak-count cap applies, and ~84% on dense real
  spectra, where the two tools break intensity ties differently at the cut
- where they disagree, tags found by both are 2–4× more likely to be confirmed
  than either tool's unique tags

`ctest` covers the rank-sum DP against exhaustive enumeration, end-to-end tag
recovery from synthetic spectra, flanking-mass placement, extension,
determinism, and the sequence filter (soundness and completeness against a naive
oracle, orientation, isobaric collapse, the derived length floor, degenerate
input).

## Licence and provenance

MIT — see [LICENSE](LICENSE). Dependencies and their licences are listed in
[BOM.md](BOM.md), together with what FasTag deliberately does *not* implement
because OpenMS or Boost already provides it.

The DirecTag algorithm is reimplemented from its publication; the reference
implementation (Apache-2.0) was read to resolve semantics the paper leaves
implicit, but no code was copied.
