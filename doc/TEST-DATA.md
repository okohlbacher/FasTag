# Test data

Two tiers. `ctest` needs none of the second — it runs on synthetic spectra built
in-process, so the suite is self-contained and offline.

## Tier 1 — ctest fixtures (in-repo, synthetic)

Built by `synth()` in `test/fastag_test.cpp` from a known peptide via
`AASequence`, so residue masses come from the same OpenMS the tool uses. If the
two ever disagreed the tests would catch it.

| what | why |
|---|---|
| ideal b/y ladders, 5 peptides | top-ranked tag must read the peptide |
| + 60 noise peaks | precision must stay above 50% |
| y-only, one interior ion deleted | gap must span the hole, in the right order |
| repeated identical runs | determinism |
| empty / 2-peak / no-precursor / charge 0 | must not crash or fabricate |

**There is deliberately no committed mzML fixture.** The local file used during
development (`local_indexed.mzML`) carries an `indexedmzML` wrapper *without* an
`<indexList>`, so `openFile` returns false and the run silently takes the
full-load path. A fixture like that makes a streaming test pass while never
exercising streaming — it happened here. Any committed mzML fixture must be
verified to contain a real `<indexList>`, e.g. via
`FileConverter -write_scan_index true`.

## Tier 2 — benchmark corpus (external, not in-repo)

Four acquisitions. **One run per instrument**, so this is coverage of the input
space, not a sampled benchmark. Every instrument below was read out of the
file's own metadata, not from a catalogue — see the Eclipse entry for why.

| # | file | instrument | analyser | acquisition | spectra | peaks/MS2 (p25/med/p75) | tolerance |
|---|---|---|---|---|---|---|---|
| 1 | S23_sage.mzML | Bruker timsTOF | TOF | DIA to diaTracer pseudo-MS2 | 632,677 | 500/500/500 (clamped) | 20 ppm |
| 2 | HeLa ddaPASEF | Bruker timsTOF Pro | TOF | real ddaPASEF | 371,187 | 95/163/273 | 20 ppm |
| 3 | Astral label-free | Orbitrap Astral | FTMS | real DDA | 102,236 | 956/1435/1682 | 20 ppm |
| 4 | Eclipse DDA-TMT | Orbitrap **Eclipse** | **ion trap** | real DDA, TMT | 41,407 | 183/273/371 | **0.3 Da** |

### Provenance

1. **S23** — local, diaTracer output. Note this is a *constructed* spectrum type:
   diaTracer aggregates correlated signal across the DIA window and retention
   time. Its 500-peak count is a hard output cap, not a property of the ions,
   which is why it is the wrong dataset to tune a peak budget on.
2. **HeLa ddaPASEF** — PRIDE `PXD054073`,
   `HeLa_T2_DDA_E_A_nLC_timspro_R02.d.7z` (3.17 GB). Confirmed by its own
   metadata: `MethodName: DDA PASEF-standard_1.1sec_cycletime.m`. Convert with
   `brfp convert -i <.d> -f 2 -o out.mzML` (writes indexed mzML directly).
3. **Astral label-free** — PRIDE `PXD076528`,
   `7673_YD12_P116861_S00_R1.raw` (2.4 GB). Human, DDA, label-free, tryptic.
   Convert with `ThermoRawFileParser -i in.raw -b out.mzML -f 2`.
4. **Eclipse DDA-TMT** — MassIVE `MSV000096674` (no PXD accession),
   `ec04479_qy_4cell_SanJose_A1.mzpeak`.
   **The mzPeak catalogue labels this "Orbitrap Astral". The file says
   `MS:1003029 name="Orbitrap Eclipse"`.** The deposit contains Astral, Eclipse
   and Ascend data and the index labels by deposit. Trusting it, and then running
   ion-trap MS2 at 20 ppm, produced 3,007 tags where 0.3 Da produces 824,959 — a
   factor of 274 that looked exactly like bad data. Read the instrument out of
   the file.

### What each one is for

- **1 vs 2** separates *DIA-derived* from *real DDA* on the same vendor. Answers
  whether pseudo-spectra flatter the tagger. (They do not: 2 is more laddered.)
- **2 vs 3** separates *vendor and analyser* at similar acquisition. Also the
  only pair whose peak density genuinely varies, so the only one that can test
  whether a peak budget generalises.
- **4** is the low-resolution cornerstone. Every other file is high-res; without
  one ion-trap acquisition the tolerance-mismatch failure mode is invisible.
- **1** additionally has 14,867 Sage PSMs at 1% FDR as ground truth, which is
  what every recall number in the README is measured against.

### Reproducing a measurement

```bash
FasTag -in <file> -out tags.tsv -tag_length 6 -threads N \
       -fragment_tolerance <see table> -fragment_tolerance_unit <ppm|Da>
```

Ground-truth recall counts spectra gaining a correctly placed tag: a tag is
correct if its sequence occurs in the Sage-identified peptide at a position whose
prefix mass matches `nterm_mass` (y-derived), or whose reverse matches with
`cterm_mass + H2O` (b-derived), within 0.05 Da.

### Known gap

No test, in either tier, exercises `-out_spectra`. It is the one output path
with no coverage.
