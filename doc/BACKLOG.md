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

## Gapped tags are overranked

Rank-1 recall is 81.5% for contiguous tags against 33.3% with `-gaps 1`, and
enabling gaps makes the single best tag *worse* in absolute terms (2,836 to
1,676 spectra). The mechanism is understood: a gap's mass is the minimum over
~190 candidate pairs, and the m/z-fidelity null prices it as though it were one
free measurement, so gapped tags are systematically flattered.

A proper fix needs a decoy-calibrated gap penalty, or scoring the two families
separately before they compete for the per-spectrum quota. Both are research, not
patches. Until then the feature is off by default and the README says to consume
the top-N list rather than the top-1 when gaps are on.

## Peak-cap defaults rest on one acquisition type

Raising `max_peaks` from 100 to 500 gained 56% more spectra with a correctly
placed tag — but only measured on diaTracer pseudo-MS2, whose peak count is
clamped flat at 500 so the cap always binds. Real ddaPASEF runs 95-273 peaks per
spectrum, where a 100-peak cap will rarely bind at all.

`-peaks_per_window` exists as the density-adaptive alternative and is at least as
good where it could be measured, but the defaults were deliberately left alone
rather than tuned on data whose defining property is being unusually dense.

Redo the cap and windowing sweeps across the four benchmark datasets before
changing any default.
