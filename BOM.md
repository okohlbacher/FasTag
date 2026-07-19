# Bill of Materials

FasTag v0.2.0 · MIT · <https://github.com/okohlbacher/FasTag>

## Direct dependencies

| Component | Version | Licence | Why | Linkage |
|---|---|---|---|---|
| [OpenMS](https://www.openms.de) | ≥ 3.5 | BSD-3-Clause | mzML and FASTA I/O, residue chemistry, physical constants, spectrum data structures, TOPP application framework | dynamic |
| [Boost.Math](https://www.boost.org) | ≥ 1.81 | BSL-1.0 | chi-squared and hypergeometric distributions | header-only |

Boost arrives transitively with OpenMS; FasTag adds no dependency OpenMS does
not already require. There is no vendored third-party source in this repository.

## What OpenMS provides, so FasTag does not implement it

| Facility | OpenMS API |
|---|---|
| Residue masses (19, I folded to L) | `ResidueDB::getResidues("Natural19WithoutI")`, `Residue::getMonoWeight(Residue::Internal)` |
| Proton mass | `Constants::PROTON_MASS_U` |
| Water mass | `EmpiricalFormula("H2O").getMonoWeight()` |
| Nearest peak within a tolerance | `MSSpectrum::findNearest(mz, tol)` |
| mzML read and write | `FileHandler::loadExperiment` / `storeExperiment` |
| FASTA parsing | `FASTAFile::load` |
| Spectra, peaks, precursors | `MSExperiment`, `MSSpectrum`, `Peak1D`, `Precursor` |
| CLI, INI files, logging, citations | `TOPPBase` |
| χ² survival function | `boost::math::chi_squared` |
| Hypergeometric tail | `boost::math::hypergeometric_distribution` |

## What FasTag implements, and why it is not borrowed

| Component | Why not from a library |
|---|---|
| Rank-sum null (`ranksum.h`) | Neither OpenMS nor Boost has the exact null of a sum of *distinct* ranks. Boost offers the normal approximation to Mann-Whitney, unusable here because the informative tags sit in the far tail. |
| Top-N peak selection | `NLargest` uses `sortByIntensity()`, an unstable sort with no secondary key. Peaks carry integer counts, so ties at the cut are common — 23% of real spectra have one — and which peak survives would depend on the standard library's sort. Peak selection provably drives tag output, so FasTag uses a total order (intensity desc, m/z desc). |
| Spectrum graph, tag enumeration, scoring, extension | The algorithm itself. `OpenMS::Tagger` enumerates tags but returns bare strings with no scoring and no flanking masses, and expands I/L during recursion, which costs 2^k tags. |
| Isobaric collapse rules | Derived at runtime from `ResidueDB` masses at the configured tolerance. `MassDecompositionAlgorithm` decomposes a mass into compositions, but returns unordered compositions rather than the ordered residue pairs the rules need. |
| k-mer index for the FASTA filter | Application-specific; nothing equivalent in OpenMS. |

## Provenance

The DirecTag algorithm (Tabb et al., *J. Proteome Res.* 2008, 7:3838,
doi:10.1021/pr800154p) is reimplemented from its publication. The reference
implementation (`freicore`/DirecTag, Vanderbilt University, Apache-2.0) was read
to resolve semantics the paper leaves implicit — the y-ion flanking-mass
convention, the fitted rather than observed peak positions, and
`numFragmentChargeStates = max(1, charge − 1)`. **No code was copied**, and
FasTag is not a derivative work of it: the scoring nulls, the tag enumerator and
the data structures are independently written, and several deliberately differ
(exact hypergeometric instead of a tabulated MVH; Fisher DOF from the number of
active subscores rather than from a sum of weights).

## Build-time only

CMake ≥ 3.21, a C++17 compiler. `cmake/eigen3-range-shim/` contains two CMake
files that work around Eigen 3.4's config-version file rejecting the version
range `OpenMSConfig.cmake` requests; they contain no Eigen code.
