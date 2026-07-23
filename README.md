# FasTag — a fast parallel mass spectrometry tagger and tag-based filter

Partial sequence tags from peptide MS/MS spectra, and a filter that keeps only
the spectra whose tags occur in sequences you supply. A reimplementation of the
[DirecTag](https://doi.org/10.1021/pr800154p) algorithm as an
[OpenMS](https://www.openms.de) TOPP tool. The citation — Tabb et al.,
*J. Proteome Res.* 2008, 7:3838 — is for that original DirecTag paper; FasTag
has no separate publication of its own.

Scales across threads with memory independent of file size, and returns 5-8%
more tags than the reference implementation.

## What it does

A *sequence tag* is a short peptide read inferred directly from a fragment-ion
ladder — no database search, no precursor digest. FasTag is a **fast, sensitive
prefilter**: it tells you what is worth looking at before you commit to a slow
method. Typical uses:

- Narrow a huge search space to the spectra that carry a tag matching your
  protein(s) of interest, via `-fasta`.
- Get a sense of what is in a run — species, contaminants, a specific protein —
  without running a full database search.
- Feed tags into a downstream tool that consumes them (open-modification search,
  spectral clustering, de novo seeding).

It is a prefilter, not a final identification method: it is tuned for recall,
not specificity — a tag matching a sequence is a reason to look closer, not a
identification in itself.

## How it works

Each MS2 spectrum's peaks form a graph: an edge connects two peaks whose m/z
difference matches a single amino acid residue's mass within the fragment
tolerance. A **tag** is a walk through that graph — a chain of residue-mass
edges — read as a short partial sequence.

**Seeding.** Every walk of length `-tag_length` (default 3 residues) is a
candidate seed.

**Scoring.** Each candidate is scored on three independent properties, combined
by Fisher's method into a single chi-squared-derived E-value:

- **Intensity** — how unlikely it is that peaks this strong would form a valid
  tag by chance, given the spectrum's overall intensity distribution.
- **m/z fidelity** — how closely each edge matches its residue mass.
- **Complementarity** — whether a tag's flanking masses are corroborated by a
  complementary ion elsewhere in the spectrum (b/y pairs summing to the
  precursor mass).

DirecTag computes the intensity term by enumerating every `C(n, k)` subset of
peak ranks — combinatorial in tag length, and impractical past length 4 or 5.
FasTag computes the same quantity as a restricted-partition count via dynamic
programming, `O(k·n·s_max)` — polynomial instead of combinatorial, verified
identical to exhaustive enumeration. This is most of why FasTag is faster, and
what makes tag lengths above four practical at all.

**Extension** (`-extension`). A seed is a *minimum* length: each scored seed can
be walked further along the graph at either end and rescored at its realised
length, recovering longer tags a fixed seed length alone would miss.

**Gaps** (`-gaps 1`). A tag can cross one missing peak — two peaks separated by
a two-residue combined mass bridge the hole, with every amino-acid composition
matching that sum spelled out as a candidate. Measured against Sage ground
truth, this lifts the number of spectra gaining a correctly placed tag by
44.7%. Because a gapped tag asserts an unobserved split, it is systematically
over-ranked relative to its real correctness (see `-gap_penalty` below).

**Sequence filtering** (`-fasta`). Reported tags can be restricted to ones
occurring in proteins you supply. Matching:

- folds I/L (isobaric — an emitted `L` says nothing about which residue it
  really was),
- tries both orientations (tags are stored N→C under a y-ion assumption, so a
  b-ion-derived tag reads reversed against the protein),
- accepts isobaric residue-pair substitutions within `-isobaric_tolerance`
  (default 0.04 Da): `N`=`GG`, `Q`=`GA`/`AG`, `R`=`GV`/`VG`, `W`=`AD`/`GE`/`SV`,
  `K`=`GA` at the default tolerance.

A short tag occurs in almost any protein by chance, so a minimum tag length for
filtering is derived from database size unless you set `-min_filter_length`.
Every filtered run reports matches beside the number expected by chance, so a
count is never mistaken for signal:

```
  len         tags      matched    by chance enrichment
    4        17454          316          291       1.1x  <- at chance
    7          212          212            2      86.0x
    9          108          108            0     999.0x
```

## Install

Prebuilt binaries for every supported platform, from the
[latest release](https://github.com/okohlbacher/FasTag/releases/latest):

| platform | download |
|---|---|
| Linux x64 | [FasTag-linux-x64.tar.gz](https://github.com/okohlbacher/FasTag/releases/latest/download/FasTag-linux-x64.tar.gz) |
| Linux arm64 | [FasTag-linux-arm64.tar.gz](https://github.com/okohlbacher/FasTag/releases/latest/download/FasTag-linux-arm64.tar.gz) |
| macOS x64 | [FasTag-macos-x64.tar.gz](https://github.com/okohlbacher/FasTag/releases/latest/download/FasTag-macos-x64.tar.gz) |
| macOS arm64 | [FasTag-macos-arm64.tar.gz](https://github.com/okohlbacher/FasTag/releases/latest/download/FasTag-macos-arm64.tar.gz) |
| Windows x64 | [FasTag-windows-x64.zip](https://github.com/okohlbacher/FasTag/releases/latest/download/FasTag-windows-x64.zip) |

Each archive extracts to a `FasTag/` folder — run `FasTag/FasTag` (Linux/macOS)
or `FasTag/FasTag.bat` (Windows); everything else inside is a bundled
dependency the wrapper needs, not something to run directly.

Binaries are not yet code-signed, so macOS Gatekeeper and Windows SmartScreen
will warn on first run; see [doc/BACKLOG-ci.md](doc/BACKLOG-ci.md) if that
matters for your deployment.

Building from source needs OpenMS ≥ 3.5 and CMake; see
[doc/BACKLOG-ci.md](doc/BACKLOG-ci.md) for what each platform requires.

## Usage

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

# mzPeak in, mzPeak out -- any combination of mzML and mzpeak works
FasTag -in run.mzpeak -out tags.tsv
FasTag -in run.mzML   -out tags.tsv -out_spectra hits.mzpeak
```

**Set the fragment tolerance to match the analyser.** A high-resolution
tolerance on low-resolution data is silent and looks exactly like bad data — on
an ion-trap file, 20 ppm returned 3,007 tags where 0.3 Da returned 824,959.
FasTag infers resolution from peak spacing and warns when the two disagree, but
it cannot know the analyser, so the setting is yours.

### mzPeak

[mzPeak](https://github.com/OpenMS/mzpeak) is a Parquet-backed format (Parquet
tables in a ZIP container). FasTag **reads and writes** it: `-in` accepts
`.mzpeak`, `-out_spectra` writes it, and all four in/out combinations work.
Tagging is unaffected by which you use — reading a run as mzpeak gives the same
spectra as reading it as mzML, and a run written to mzpeak and tagged again
reproduces the original tag set exactly.

The released binaries have it. Building it yourself needs an OpenMS with mzPeak
support ([OpenMS-mzPeakRW](https://github.com/okohlbacher/OpenMS-mzPeakRW)) —
stock and bioconda OpenMS do not have it, and configure reports which you got:

```
-- FasTag: mzPeak input available (-in accepts .mzpeak)
```

**Not memory-bounded, upstream** — `MzPeakFile::transform()` currently
materialises the whole run rather than streaming it, unlike the O(threads)
mzML path: a 155 MB `.mzpeak` peaks around 1.7 GB. Prefer mzML for runs large
relative to RAM; see [doc/BACKLOG-mzpeak.md](doc/BACKLOG-mzpeak.md).

## Command-line reference

| Option | Default | Meaning |
|---|---|---|
| `-in <file>` | — | Input spectra: mzML or mzpeak |
| `-out <file>` | — | Tag list, tab-separated |
| `-fasta <file>` | none | Report only tags occurring in these sequences |
| `-out_spectra <file>` | none | Write spectra carrying a reported tag here, as mzML or mzpeak (by extension). Needs memory proportional to the *file*, not the thread count, unlike every other path |
| `-tag_length <n>` | 3 | Seed tag length in residues |
| `-extension <n>` | 0 | Max residues appended per terminus; 0 disables extension |
| `-gaps <n>` | 0 | Allow a tag to cross one missing peak (0 or 1) |
| `-deisotope` | off | Collapse isotope clusters to their monoisotopic peak and move multiply-charged fragments onto the singly-charged scale before peak selection |
| `-fragment_tolerance <value>` | 20 | Fragment mass tolerance |
| `-fragment_tolerance_unit <ppm\|Da>` | ppm | Tolerance unit |
| `-max_peaks <n>` | 400 | Peaks retained per spectrum; 0 uses the internal ceiling of 1024, not unlimited. The ceiling for `-peaks_per_window` |
| `-peaks_per_window <n>` | 10 | Keep this many peaks per 100 Da window instead of the strongest `-max_peaks` overall; 0 disables |
| `-max_tags <n>` | 50 | Tags reported per spectrum; 0 = unlimited |
| `-max_evalue <value>` | 20 | E-value cutoff; 0 disables |
| `-gap_penalty <value>` | 100 | Rank gapped tags as if their E-value were this many times worse. Affects order only, never which tags are reported — gapped tags are otherwise heavily over-ranked (~95% of top-1 slots while ~3.6x less likely to be correct). 1 disables |
| `-isobaric_tolerance <value>` | 0.04 | Isobaric residue-pair substitution tolerance when matching `-fasta`; 0 requires exact strings |
| `-min_filter_length <n>` | 0 | Ignore tags shorter than this when matching `-fasta`; 0 derives a floor from database size |
| `-orientation <both\|forward>` | both | Also match a tag reversed (b-ion reading), or only as written |

Plus the standard OpenMS TOPP options: `-threads <n>` (parallelism — the
core performance lever), `-ini <file>` / `-write_ini <file>` (parameter files),
`-log <file>`, `-no_progress`, `--help` / `--helphelp`.

### Output

One TSV row per tag: `spectrum` (native ID), `tag` (sequence), `length`,
`charge` (fragment charge), `nterm_mass` / `cterm_mass` (flanking masses,
4-decimal precision), `extended` / `gapped` (flags), `evalue`, and `fasta_hit`
(`fwd` / `rev` / `-`, present only with `-fasta`).

## Demo data

FasTag was validated against public runs from
[ProteomeXchange](https://www.proteomexchange.org)/PRIDE and MassIVE, spanning
two vendors, three instruments, and both high- and low-resolution MS2:

| dataset | instrument | accession | file | convert with |
|---|---|---|---|---|
| HeLa, real ddaPASEF | Bruker timsTOF Pro | [PXD054073](https://www.ebi.ac.uk/pride/archive/projects/PXD054073) | `HeLa_T2_DDA_E_A_nLC_timspro_R02.d.7z` (3.17 GB) | `brfp convert -i <.d> -f 2 -o out.mzML` |
| Label-free DDA | Orbitrap Astral | [PXD076528](https://www.ebi.ac.uk/pride/archive/projects/PXD076528) | `7673_YD12_P116861_S00_R1.raw` (2.4 GB) | `ThermoRawFileParser -i in.raw -b out.mzML -f 2` |
| DDA-TMT | Orbitrap Eclipse (ion trap MS2) | [MSV000096674](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000096674) (MassIVE, no PXD) | `ec04479_qy_4cell_SanJose_A1.mzpeak` | already mzpeak |

The Eclipse file's own metadata says `Orbitrap Eclipse` even though the deposit
is catalogued as Astral — read the instrument out of the file itself before
picking a tolerance. Run it with `-fragment_tolerance 0.3
-fragment_tolerance_unit Da`; the other two are high-resolution and want the
20 ppm default. See [doc/TEST-DATA.md](doc/TEST-DATA.md) for the full
provenance and what each dataset is useful for comparing.

## Performance

Scaling is near-linear to 16 cores and useful well beyond; output is identical
at every thread count, and memory is O(threads) rather than O(file) on the
default mzML path (`-out_spectra` is the one exception — see above). Exact
wall-clock time and peak memory depend heavily on the machine, so they aren't
quoted here — measure on your own hardware and data.

Against the reference implementation, run sequentially on the same file and
hardware, FasTag keeps scaling well past where DirecTag flattens out: DirecTag
rebuilds its ranksum table serially on every run, while FasTag's equivalent is
a one-off cost paid once regardless of thread count.

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
correctly placed tag:

| | spectra | vs baseline | correct tags |
|---|---|---|---|
| default | 3,480 | — | 8,054 (7.4%) |
| `-deisotope` | 4,142 | +19.0% | 12,683 (10.6%) |
| `-gaps 1` | 5,055 | +45.3% | 33,236 (7.1%) |
| both | **5,976** | **+71.7%** | 45,919 (8.9%) |

`ctest` covers the rank-sum DP against exhaustive enumeration, end-to-end tag
recovery from synthetic spectra, flanking-mass placement, extension, gap
ordering with a negative control, per-row mass closure, determinism, the
sequence filter, and that `-gap_penalty` reorders without ever removing a tag.

### Ranking, and the synthetic benchmark

Whether a correct tag exists is not the same question as whether it ranks
first, and PSM ground truth only answers the first one. `bench/benchmark.cpp`
generates synthetic spectra with exact per-tag truth, fitted to four real
acquisitions' peak count and ladder edge density (see
[doc/TEST-DATA.md](doc/TEST-DATA.md)), to answer the second.

It found that gapped tags took **95% of rank-1 slots while being 3.6x less
likely to be correct** — enabling gaps made the single best tag worse even as
it doubled total recall. `-gap_penalty` (default 100) fixes the ordering
without filtering anything out:

| astral profile, TL=4 | rank-1 | top-5 | total recall |
|---|---|---|---|
| `-gaps 1 -gap_penalty 1` | 17.7% | 47.7% | 71.3% |
| `-gaps 1` (default 100) | **28.6%** | **52.2%** | 71.5% |

## Known limitations

- **Some defaults are conservative.** `-deisotope` and `-gaps` are off. The peak
  budget defaults to `-peaks_per_window 10 -max_peaks 400`, validated on real
  ground truth (PXD000001) as neutral-to-positive versus the old flat 100-peak
  cap at ~1.35x the runtime; revert with `-peaks_per_window 0 -max_peaks 100`
  for the old speed. See [doc/BACKLOG.md](doc/BACKLOG.md).
- **No modification support.** Residues are the unmodified 19, so labelled
  samples (TMT and similar) will not match tags spanning a modified residue.
- **mzPeak is not memory-bounded on read** (see above), an upstream property of `MzPeakFile::transform()` rather than of FasTag.

## Licence and provenance

MIT — see [LICENSE](LICENSE). Dependencies and licences are in [BOM.md](BOM.md).

The DirecTag algorithm is reimplemented from its publication. The reference
implementation (Apache-2.0) was read to resolve semantics the paper leaves
implicit; **no code was copied**.
