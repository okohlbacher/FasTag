# Backlog

## mzPeak input — designed, not built

**Decision: the consumer + work-queue adapter, when it is built. Not now.**

FasTag reads and writes mzML only. [mzPeak](https://github.com/OpenMS/mzpeak) is
a Parquet-backed format whose columnar layout would make the metadata problem
that forced the OpenMS patch (see `OPENMS-FAST-READER.md`) disappear structurally
rather than by patching: MS level and precursor live in their own columns and can
be read without touching peak data.

### What already works

`OpenMS::MzPeakFile` (branch `feature/mzpeak-file-handler`) reads and writes the
fields FasTag needs. Verified in its source, not assumed:

```cpp
auto charge = ...->GetFieldByName("MS_1000041_charge_state");
if (charge && !charge->IsNull(r)) { pd.has_charge = true; pd.charge = charge->Value(r); }
```

An earlier note here claimed charge does not survive mzPeak conversion. That was
wrong — it was a property of the *files* tested, which came from a converter that
does not populate the field, not of the format or of OpenMS.

`MzPeakFile` also offers `transform(filename, IMSDataConsumer*)`, so a large file
need not be resident.

### The one real obstacle

Not I/O — the parallel model. FasTag's OpenMP loop is **pull-based random
access**: each worker copy-constructs its own reader and calls `getSpectrum(i)`.
`MzPeakFile` offers **push-based streaming** (`transform` feeds a consumer
sequentially) or a whole-file `load`. Neither fits directly.

### Chosen approach: consumer + bounded work queue

Write an `IMSDataConsumer` whose `consumeSpectrum()` pushes onto a bounded queue
that the existing workers drain. Memory stays bounded, no OpenMS change is
needed, and the tagging code is untouched.

The cost is output ordering. Rows currently come out in spectrum order because
each worker writes to `rows[i]` under an index it already has; a queue loses that
index unless it is carried explicitly. Carry it — determinism of the output is
worth more than the simplicity of a bare queue, and `ctest` asserts it.

Rejected alternatives:

- **`OnDiscMzPeakExperiment`** — mirror `OnDiscMSExperiment`'s random-access API
  over Parquet row groups. Architecturally the right answer, and FasTag would
  change by about five lines. But it is an OpenMS contribution of roughly a week,
  and it belongs upstream where every TOPP tool benefits, not here.
- **Reading Parquet directly in FasTag** — adds Arrow and Parquet to a tool whose
  dependency discipline is documented in `BOM.md`, and duplicates what OpenMS
  already does.

### Why it is deferred

- mzPeak is explicitly pre-1.0: *"no stability is guaranteed at this point. Our
  current goal is to stabilize the API by the end of the summer (2026)."*
- `MzPeakFile` is on a feature branch, not OpenMS mainline. Adopting it now would
  make FasTag depend on two unmerged branches at once.

Revisit when the format stabilises and `MzPeakFile` lands in mainline. The cheap
first step, whenever that is: round-trip a benchmark mzML through
`MzPeakFile::store` and back, and confirm precursor m/z and charge survive.

## Gapped tags were overranked — FIXED, with a caveat

`-gap_penalty` (default 100) multiplies a gapped tag's E-value for **ranking
only**; the cutoff still sees the uncorrected value, so no tag is ever removed.
Rank-1 accuracy improved in all 12 profile × tag-length cells of
`bench/benchmark.cpp`, and on real ddaPASEF the tag set is byte-identical while
the rank-1 gapped share falls 96.4% → 76.7%.

**The caveat**: 100 is calibrated on synthetic spectra, not derived. The
multiplicity argument (a gap edge picks from 190 pairs where an ordinary edge
tests 19 residues) justifies only ~10; measurement says the true over-scoring is
roughly an order of magnitude larger, and rank-1 accuracy is still creeping up at
500. The curve has no clean optimum, which is the signature of "just rank all
contiguous tags above all gapped ones" — near enough true, since a contiguous
rank-1 tag is right ~55% of the time against ~15% for a gapped one.

Worth revisiting **only with real ground truth**. Regenerating the Sage PSMs (see
`TEST-DATA.md`) is the prerequisite; until then this constant rests on a
generator this project wrote to test itself, which is a real limitation.

Related, and unfixed: within the gapped family the E-value carries almost no
information about correctness — r1-gapped accuracy sits at 15-20% regardless of
penalty. Pushing the family down helps; ordering *inside* it does not seem to.

## Gap-crossing extension — MEASURED, rejected

Task "rework gaps as merge semantics" is closed on the evidence. Extension
crossing a hole was implemented fully (reverse gap edges, shared gap budget, gap
tried only where no ordinary edge continues) and lost in all six cells tested:
rank-1 down up to 13.1 points at extension 4, with **zero** recall gain. See the
comment on `extendPath` for the numbers. The stale `wip/gap-merge` branch is
superseded and would revert a large amount of later work if merged; it is kept
only as provenance.

Subsumption (collapsing a contiguous tag into a gapped superset) was **not**
implemented. Its motivation was reducing redundancy between the two families,
which `-gap_penalty` now addresses more directly by separating them in the
ranking. Reconsider only if redundancy is shown to cost something measurable.

## Peak-cap defaults — sweep DONE, default deliberately unchanged

The sweep this entry asked for has been run across all four acquisition classes
(`bench/benchmark.cpp`, TL=4, gaps on).

**The flat 100-peak cap is strongly dataset-dependent — the defect suspected here
is real, and larger than expected.** rank-1 / total recall:

| cap | astral | ddapasef | diatracer | iontrap |
|---|---|---|---|---|
| 50 | 12.9 / 35.1 | 37.6 / 78.2 | 26.6 / 64.7 | 14.5 / 62.2 |
| **100 (default)** | **27.5 / 69.9** | **45.7 / 87.4** | **44.3 / 86.6** | **22.7 / 76.3** |
| 200 | 46.6 / 90.8 | 46.5 / 89.3 | 55.7 / 94.7 | 25.3 / 80.8 |
| 400 | 63.3 / 97.2 | 46.7 / 89.3 | 60.3 / 96.0 | 25.3 / 80.9 |
| 800 | 67.3 / 97.8 | 46.7 / 89.3 | 60.2 / 96.2 | 25.3 / 80.9 |

Sparse classes saturate by 200 because their spectra have no more peaks to give.
Dense ones do not saturate until 400-800: **the default starves Astral data by a
factor of two in rank-1 accuracy.** Where a run lands on this table is decided
entirely by the instrument, which is the dependence worth removing.

`-peaks_per_window` removes it, as designed. The quota binds only where there is
density to spend it, so sparse classes barely move while dense ones gain a lot:

| setting | astral | ddapasef | diatracer | iontrap |
|---|---|---|---|---|
| cap 100 | 27.5 / 69.9 | 45.7 / 87.4 | 44.3 / 86.6 | 22.7 / 76.3 |
| ppw 10, cap 400 | 44.3 / 84.6 | 45.8 / 88.8 | 50.9 / 89.6 | 23.8 / 80.6 |
| ppw 15, cap 800 | 54.8 / 92.6 | 46.6 / 89.2 | 56.6 / 94.8 | 24.6 / 80.8 |

It is not free. Real-file wall time, 8 threads:

| | astral_lf (102k spectra) | ddaPASEF (371k spectra) |
|---|---|---|
| default | 10.3 s | 35.2 s |
| ppw 10 / cap 400 | 20.3 s | 44.0 s |
| ppw 15 / cap 800 | 27.7 s | 49.6 s |

### Decision: default unchanged, pending real ground truth

Deliberate, and the conservative call. Every number above comes from spectra this
project generated to test itself; none is a measurement on real data with real
identifications. Changing a released tool's output and doubling its runtime is
not something synthetic evidence should decide, however consistent it is.

**Prerequisite to revisiting**: regenerate the Sage ground truth. That needs Sage
0.14.7, a human FASTA with reversed decoys, the S23 diaTracer file, and the
precursor-injection fix in `vault/03-Design/SAGE-Tag-Confirmation.md` (Sage
otherwise panics with `missing MS1 precursor for frame=1`, because diaTracer
writes precursor m/z only on the isolation window). None of it is currently on
this machine — see `TEST-DATA.md`. Rerun the two tables above against PSMs and
the decision makes itself.

**Recommendation when that happens**: `-peaks_per_window 10 -max_peaks 400`. Most
of the available gain, at the knee of the time curve, and its worst case across
the four classes still beats today's default. Users who know their data is dense
can pass it today.
