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

## Gapped tags were overranked — FIXED, and now validated on real data

`-gap_penalty` (default 100) multiplies a gapped tag's E-value for **ranking
only**; the cutoff still sees the uncorrected value, so no tag is ever removed.
Rank-1 accuracy improved in all 12 profile × tag-length cells of
`bench/benchmark.cpp`, and on real ddaPASEF the tag set is byte-identical while
the rank-1 gapped share falls 96.4% → 76.7%.

**Real-ground-truth validation, done.** The synthetic calibration this entry used
to caveat has been checked against a real search: [PXD000001](
https://www.ebi.ac.uk/pride/archive/projects/PXD000001), the canonical
ProteomeXchange demo dataset (Erwinia carotovora, TMT 6-plex, LTQ Orbitrap
Velos HCD) — chosen because it ships a Mascot `.dat` search of its own mzML,
giving real per-spectrum PSMs rather than a re-run search of unknown quality.
Ground truth: top hit per query from the `.dat`'s `peptides`/`decoy_peptides`
sections, target-decoy FDR ≤ 1% (score ≥ 16.85, 2,254 of 6,103 queries).
Tag correctness re-derives `nterm_mass`/`cterm_mass` against the identified
peptide's real fragment masses, which for this dataset means accounting for the
search's fixed mods (TMT6plex on the N-terminus and every K, Methylthio on every
C) — omitting them silently zeroes out every real match, which is exactly what
happened on the first pass.

`-tag_length 4 -gaps 1`, threads varied only for speed:

| `-gap_penalty` | rank-1 | top-5 | total recall | gapped share of rank-1 | gapped rank-1 correct |
|---|---|---|---|---|---|
| 1 | 69.0% | 90.8% | 98.05% | 34.7% | 25.4% |
| 10 | 83.0% | 95.2% | 98.05% | 10.3% | 27.6% |
| 30 | 85.1% | 95.7% | 98.05% | 6.0% | 29.9% |
| **100 (default)** | **86.4%** | **96.1%** | **98.05%** | **3.2%** | **27.8%** |
| 300 | 86.8% | 96.2% | 98.05% | 2.2% | 24.0% |
| 1000 | 86.9% | 96.4% | 98.05% | 2.0% | 22.2% |
| 10000 | 87.1% | 96.4% | 98.05% | 1.6% | 18.9% |

Total recall is flat, exactly as designed — the penalty reorders and never
removes a tag. **Real data confirms the shape the synthetic curve predicted**:
rank-1 keeps climbing with no clean optimum, and 100 already sits at its knee —
10x to 1000 buys +0.5 rank-1 points, another 10x to 10000 buys +0.2 more. The
multiplicity argument that originally suggested ~10 undersold it; measurement
(synthetic and now real) agrees the true over-scoring is roughly an order of
magnitude larger, which is where 100 came from. **No change to the default**:
real evidence supports the existing calibration rather than overturning it.

**Scope of this validation, stated plainly**: one dataset, one (older,
TMT-labelled) instrument, 2,254 PSMs. The gapped-family rank-1 counts at high
penalty are small (37-72 spectra), so the 19-30% correctness range within that
family is noisy, not a precise curve. Confirms direction and rough magnitude;
does not replace validation on additional real datasets if the question is
reopened.

**A second, non-TMT dataset was attempted and found too sparse to add a real
data point.** [PXD059878](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/10/PXD059878/)
(`thermo-ltq-xl-iontrap`, a PC4-acetylation/phosphorylation site-mapping study,
label-free, real linear-ion-trap CID) ships a Proteome Discoverer `.msf` — a
SQLite database, directly queryable — with 53 high-confidence (PD "High",
Rank 1) PSMs. Real ground truth, correctly scored (including per-PSM Phospho/
Oxidation position deltas, not the fixed-mod rule PXD000001 needed), but only
2 of the 50 matching spectra recovered a correct tag at any `-gap_penalty`,
let alone at rank 1 — not enough signal for a curve. Two compounding reasons,
both real properties of this dataset rather than a bug: it is a targeted PTM
study, so most identified peptides are short and carry a variable
modification FasTag cannot see (no modification support — a tag spanning a
phosphorylated residue simply will not be found); and it is old, low-resolution
linear-ion-trap CID, which produces markedly sparser fragment ladders than the
HCD/Orbitrap data used above even before that. Left as a documented dead end
rather than reported as a second confirmation it cannot honestly be.

Related, and now also confirmed on real data rather than only synthetic: within
the gapped family the E-value carries almost no information about correctness —
rank-1-gapped accuracy sits at roughly 19-30% regardless of penalty, in the same
range synthetic data found (15-20%). Pushing the family down helps; ordering
*inside* it does not seem to.

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

## Peak-cap defaults — CHANGED to `-peaks_per_window 10 -max_peaks 400`, validated on real data

The default is now `-peaks_per_window 10 -max_peaks 400` (was `0` / `100`). The
synthetic sweep below recommended it; a real search then confirmed the flat-cap
increase this entry originally suspected of helping does **not** help real data,
while the per-window quota does.

### Real ground truth (PXD000001, added later)

The synthetic sweep could not be trusted alone (see the original decision at the
bottom). Rerun against real PSMs — [PXD000001](https://www.ebi.ac.uk/pride/archive/projects/PXD000001),
the same Mascot-searched Erwinia TMT dataset that settled `-gap_penalty`,
2,254 PSMs at 1% target-decoy FDR, `TL=4 -gaps 1`, scored against the identified
peptides with the search's fixed mods (TMT6plex N-term/K, Methylthio C):

| setting | wall | rank-1 | top-5 | total |
|---|---|---|---|---|
| cap 50 | 3.5 s | 81.9 | 92.8 | 95.8 |
| cap 100 (old default) | 10.4 s | 86.4 | 96.1 | 98.0 |
| cap 200 | 36.6 s | 86.3 | 96.5 | 98.2 |
| cap 400 | 64.0 s | 86.0 | 96.4 | 98.3 |
| cap 800 | 97.7 s | 86.0 | 96.4 | 98.3 |
| **ppw 10 / cap 400 (new default)** | **14.1 s** | **87.0** | 96.0 | 98.0 |
| ppw 15 / cap 800 | 42.0 s | 87.0 | 96.3 | 98.1 |

Two things the real data settles that the synthetic could not:

1. **A flat cap increase is counterproductive here.** 100 → 800 leaves rank-1
   flat-to-worse (86.4 → 86.0) for +0.3 pp total recall at ~10x the runtime. This
   dataset (LTQ Orbitrap Velos, older TMT) simply is not dense enough to be
   starved the way the synthetic Astral profile is — so a naive "raise max_peaks"
   would have slowed every run for nothing.
2. **The per-window quota is the setting that helps.** `-peaks_per_window 10`
   lifts rank-1 +0.6 pp (86.4 → 87.0) at only 1.35x the time, because it keeps
   the strongest peaks *per m/z window* rather than the strongest overall.

The gain on this one (non-dense) dataset is small; the large gains remain
synthetic (dense Astral, below). But the direction now agrees across synthetic
and real — ppw helps, flat cap does not — and the cost is modest, so the change
is made rather than deferred. **Honest caveat**: one real dataset, and not a
dense one; the strong dense-instrument gain is still synthetic-only, since no
dense real ground truth is on this machine. Revert with
`-peaks_per_window 0 -max_peaks 100` for the old speed.

### The synthetic sweep (what recommended it)

Run across all four acquisition classes (`bench/benchmark.cpp`, TL=4, gaps on).

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

### Original decision (superseded): default unchanged, pending real ground truth

Kept for the record. The call at the time was conservative — every number in the
synthetic tables comes from spectra this project generated to test itself, and
changing a released tool's output and runtime is not something synthetic evidence
should decide alone. That prerequisite (real ground truth) has since been met via
PXD000001 above, which is why the default was changed. The synthetic
recommendation it named — `-peaks_per_window 10 -max_peaks 400` — is exactly what
the real data confirmed, so it is now the default rather than a suggestion.

The Sage ground truth this entry originally waited for is still gone (see
`TEST-DATA.md`); PXD000001 stood in for it, as it did for `-gap_penalty`.

## Research: how MS taggers are used in the literature (2021+) — candidate features

Full survey with citations and adversarial verdicts in
[RESEARCH-tagger-uses-2021plus.md](RESEARCH-tagger-uses-2021plus.md). Distilled
here as a feature backlog. Two standing principles gate everything: tagging is a
sensitive **prefilter** (recall-biased; specificity comes from later stages), and
synthetic evidence cannot flip a shipped default.

**The real competitor is the fragment-ion index** (MSFragger/Sage/Comet-FI), not
other taggers. FasTag should not claim to beat it on ordinary open search. Its
defensible ground is (i) database-free operation, (ii) multi-PTM / mutation
reconciliation where a precursor-offset open search struggles, and (iii) being a
calibrated *tag/feature generator* for the rescoring and assembly ecosystems.

Candidate features, with verdicts (ADOPT / MAYBE / SKIP):

| # | Feature | Driving use | Verdict |
|---|---|---|---|
| F1 | Flank-mass reconciliation, delta = mod/mutation | open/blind PTM, variants | ADOPT — this is TagRecon stage A+; validated by PIPI2 (2024), Open-pFind |
| F2 | FM-/substring index for tag→DB lookup | scale of F1 | ADOPT — reuse a succinct-index lib; enabling substrate |
| F3 | Gap tags (mass-gap residues) | noisy/ion-starved spectra | ADOPT (already have `-gaps`); cap at 1 |
| F4 | **Calibrated tag confidence / per-residue FDR** | rescoring, taxonomy, variant QC | PARTLY DONE — E-value validated + per-residue `min_conf`/`mean_conf` shipped; true decoy q-value is research-grade (see below) |
| F5 | Mass-shift localization (shifted-fragment method) | open PTM; the mod-localization roadmap | ADOPT for localization; interoperate with PTM-Shepherd, don't reimplement |
| F6 | Multi-length tags as a recall/specificity knob | DB search, HLA, taxonomy | ADOPT — but default recall-y (short); long tags opt-in |
| F7 | Tag-agreement features for MS2Rescore/Oktoberfest | rescoring (+10-30% IDs) | ADOPT — cheap, high-ROI, repositions FasTag as a feature source |
| F8 | **Native tag→taxonomic/species detector** | metaproteomics | ADOPT — most differentiated; no mainstream *tagger* ships this (cf. MegaPX, NovoLign). In progress |
| F9 | Single-substitution (mutation) tags | proteogenomics | MAYBE→ADOPT — nearly free once F1 exists; emit as flagged *unvalidated* variants only |
| F10 | Low-latency single-spectrum tagging | real-time / instrument-control search | MAYBE — build the streaming path + benchmark; vendor integration out of scope |
| F11 | Glyco spectrum flag (oxonium detector) | glycoproteomics | SKIP full glyco; MAYBE the flag only |
| F12 | Top-N tags per spectrum for chimeric/DIA | chimeric DDA, DIA | MAYBE — top-N yes (already have `-max_tags`); full DIA deconvolution SKIP |
| F13 | **Standard self-describing output** (ProForma / mzTab / mzIdentML / USI) | all downstreams | ADOPT — unglamorous, force-multiplying; low cost |
| F14 | Overlap-friendly long reads for antibody assembly | mAb/repertoire de novo | MAYBE — nearly free given F4+F13 (feed Stitch/ALPS); don't build an assembler |
| F15 | Cross-link diagnostic-ion tagging | XL-MS | SKIP — crowded specialist space |

**Top 5 to prioritize** (from the report): **F4** (calibrated tag confidence —
foundational for everything trustworthy downstream), **F1+F2** (reconciliation on
a substring index — the core extension, already begun as TagRecon stage A),
**F8** (native taxonomic detector — most differentiated; in progress), **F7+F13**
(rescoring-feature output in a standard format — cheap, force-multiplying), **F5**
(mass-shift localization — on the roadmap, complements F1).

Several F1/F5/F8 items are already under way (TagRecon stage A; the modification
support in v0.13.0; the species detector on `feature/species-detector`). F4, F7
and F13 are the highest-leverage *not-yet-started* items.

## F4 tag confidence — E-value validated, per-residue shipped, decoy FDR deferred

Two of the three parts of F4 landed; the third is a documented research problem.

**Shipped.** The `evalue` was validated as a strong per-tag confidence on real
identifications (PXD000001, 1% FDR PSMs, TMT/Methylthio fixed mods): correct
tags separate from incorrect by ~400x in the median E-value (0.006 vs 2.306),
and the best-E-value tag per spectrum is correct **88.5%** of the time. Two new
output columns, `min_conf` and `mean_conf` (both in [0,1], higher better), score
each residue's own support as m/z-fit x endpoint-intensity, so `min_conf`
localises the weakest residue — correct tags carry ~0.45 median `min_conf`
against ~0.17 for incorrect, and ranking by `mean_conf` alone recovers 86.5%
rank-1, nearly matching the E-value. Cheap (computed during scoring) and
always on.

**Deferred, and why.** A calibrated per-tag *q-value / FDR* was attempted via a
target-decoy scheme (tag decoy spectra, read the false-positive rate off the
decoy E-value distribution) and **does not calibrate with a single-spectrum
decoy**. Both a uniform-random-m/z decoy and a gap-shuffle decoy produce an
essentially empty null — q collapses to ~0 for every tag — because FasTag's real
false tags are not random-peak coincidences but reads of *chimeric co-isolated
peptides* and real noise structure (the 66,516 "incorrect" tags on PXD000001
have median E-value 2.3, i.e. they are real-ladder reads that simply do not
match the one identified peptide). Destroying one spectrum's ladder cannot
reproduce those. This matches the literature: de novo / tag FDR is its own
research topic (NovoBoard, Winnow). A real q-value needs a chimeric-aware or
entrapment null, or a reconcile-to-peptide-first target-decoy — tracked, not
shipped, so no miscalibrated FDR reaches users. The `TagFDR` machinery (an
empirical target-decoy q-curve) was written and unit-tested but is not wired in
until a null that calibrates exists.
