# FasTag GUI — phased implementation plan

A cross-platform Electron desktop app wrapping the FasTag CLI: pick input files,
set parameters, save/load presets, run batches, browse tag results. The CLI is
the source of truth; the GUI shells out to it.

Synthesized from an adversarial review (codex), the `electron-app` packaging
skill, and hands-on P0 findings. Each phase is a shippable increment,
walking-skeleton first.

## Decisions made up front (the hard ones)

| Question | Decision | Why |
|---|---|---|
| **Binary integration** | **Bundle** the CLI + its libs + `share/OpenMS` via `extraResources`. Not detect-installed, not download-on-first-run. | Reproducibility and "it just works". A detected install drifts in version; download-on-first-run is package-manager hell and a signing hole. Advanced users can still override with `FASTAG_BIN`. |
| **Nested-binary signing** | Sign every Mach-O/PE in the bundle in an `afterPack` hook, **before** the outer signing pass; notarize in `afterSign`. | An unsigned bundled binary is the #1 macOS notarization rejection. Signing order is load-bearing. |
| **Large results (100k–1M rows)** | **DuckDB in an Electron utility process**, querying the TSV directly; renderer virtualizes display only. Not DuckDB-WASM, not load-all-into-DOM. | Native DuckDB reads/sorts/filters a million-row TSV out of process without freezing the UI. WASM is slower and memory-bound. P0/P1 ship a bounded streaming preview as the honest stand-in. |
| **Progress** | **Done.** CLI has a `-progress` flag emitting `FASTAG_PROGRESS done=<n> total=<n>` on stderr; the GUI parses it into a determinate bar (mzML) or indeterminate bar with a spectrum count (streaming mzPeak, total unknown). | A structured line is a real contract, not scraped prose. Off by default so the CLI stays quiet for non-GUI use. |
| **Parameter form** | Extract the CLI's `-write_ini` into a **versioned manifest at build time**; a curated UI overlay keys off canonical option names. Not a runtime-generated form. | Single source of truth, CI fails on drift. P0 proved why: a hand-picked default (`0.02` with no unit) silently meant 0.02 ppm → **zero tags**. The manifest kills that whole class of bug. |
| **Batch concurrency** | Default **concurrency = 1** (the CLI is already multithreaded); advanced: bounded by a global CPU budget. | Running N multithreaded jobs at once oversubscribes cores and is slower. Sequential is the correct default. |
| **Security** | `sandbox:true`, `contextIsolation:true`, `nodeIntegration:false`, strict CSP, typed allowlisted IPC, `shell:false` + argv array (never a shell string), external links via `shell.openExternal`. | Standard Electron hardening; already in place from P0. |
| **Repo layout** | Same repo, **`gui/` subdir**. | Atomic evolution with the CLI: a param change and its GUI overlay land in one commit. Separate release workflow, same repo. |
| **Framework** | **electron-vite + React + TS**, packaged with **electron-builder**. Not Forge (its Vite path is experimental), not plain HTML. | Cheapest maintainable path to a real form + results table; electron-builder's updater is written against its own output. |

## Phases

### P0 — architecture canary ✅ (done)
- electron-vite + React + TS scaffold in `gui/`, hardened window, typed IPC.
- **Binary canary**: resolve the bundled/dev binary, run `--help`, verify it
  executes and links its libs on *this* machine before offering to run anything.
- Findings that shaped the plan: the ppm-unit default bug; the mzPeak
  capability gap between local builds; the conda-build `libomp` loader path.

### P1 — secure walking skeleton ✅ (done)
Pick one file → set core params (tag_length, fragment_tolerance + **unit**,
max_tags) → run → streamed stderr log → **bounded** result preview (first N rows).
Cancel kills the process tree. Defaults mirror the CLI exactly.

### P2 — parameter contract + presets ✅ (presets + last-used done)
- Build-time `-write_ini` → `params.manifest.json` (name, type, default,
  restrictions, description, advanced flag). CI fails if the tool and manifest
  disagree.
- Curated UI overlay: core params up front, advanced collapsible, tooltips from
  the manifest descriptions, client-side validation from `restrictions`.
- Named presets + last-used, persisted atomically in `userData` with a
  `schemaVersion` + `cliVersion`; migrate on load.

### P3 — observable jobs
- ✅ `-progress` in the CLI; GUI shows a real progress bar (determinate for
  indexed input, indeterminate + count for streaming mzPeak).
- Run IDs + a state machine (queued→running→done/failed/cancelled).
- Write output to `*.partial`, rename on success (no half-written TSVs).

### P4 — batch ✅ (sequential queue done)
- Queue N files, per-file status, sequential by default, cancel one/all,
  output-collision handling, failure isolation (one file failing doesn't stop the queue).

### P5 — million-row results browser
- DuckDB utility process over the TSV; paged/sorted/filtered queries via IPC;
  renderer virtualization. Cache the on-disk DB keyed by source identity;
  invalidate on file change.

### P6 — distribution
- electron-builder: dmg **+ zip** (updater needs the zip) & both arches on macOS,
  nsis on Windows, AppImage on Linux.
- GH Actions matrix: build → sign nested binaries (afterPack) → sign+notarize+staple
  (macOS) / Authenticode via SignPath (Windows) → `--publish always` per leg.
- electron-updater via GH Releases; test an N→N+1 update before enabling.
- **Release-build prerequisite**: one FasTag built from current `main` against
  the mzPeak-capable OpenMS (no local build currently has both mzPeak *and* the
  latest features).

## Biggest risk
Reliable cross-platform delivery of the native CLI **+ its runtime libs + OpenMS
data + signing + cancellation** as one signed, notarized, self-updating bundle.
Everything else is ordinary app work; this is the part that bites.

## Smallest honest MVP
P1–P3 on all three platforms: a bundled, version-verified FasTag that runs one
file with the core parameters, shows progress and logs, cancels cleanly, and
displays a bounded result preview. Batch and the million-row browser are strict
add-ons on top.
