# FasTag GUI

Cross-platform desktop app (Electron) wrapping the FasTag CLI. Pick input files,
set parameters, run, and browse tag results. The CLI stays the source of truth;
the GUI shells out to it.

See [PLAN.md](PLAN.md) for the phased implementation plan and the architectural
decisions (binary bundling, million-row results, signing, auto-update).

Status: **P0/P1 — walking skeleton.** One file → core params → run → streamed
log → bounded result preview, with a hardened renderer and typed IPC.

## Develop

```bash
cd gui
npm install
./dev.sh          # electron-vite dev + the env the local binary needs
```

`dev.sh` launches the app and puts the local build's runtime libs on the loader
path. The app finds FasTag through `resources/fastag/` (git-ignored symlinks):

```
resources/fastag/bin/FasTag       -> a local FasTag build
resources/fastag/share/OpenMS     -> the matching OpenMS data dir
```

Point these at whichever build you want to drive:

- `build-mzpeak` — reads `.mzpeak` **and** `.mzML`, but an older FasTag revision
  (no subsample / calibrated confidence / species columns).
- `build-main` — latest features, `.mzML` only (its OpenMS lacks mzPeak).

No local build has *both* yet; the release will bundle one FasTag built from
current `main` against the mzPeak-capable OpenMS (see PLAN.md P6).

Override the resolved binary/data explicitly with `FASTAG_BIN` and
`OPENMS_DATA_PATH` — these win over the bundled/symlinked paths.

## Build

```bash
npm run build       # compile main + preload + renderer into out/
npm run typecheck    # tsc, both node and web projects
```

Packaging (electron-builder), signing, notarization, and auto-update land in
PLAN.md phase P6.

## Layout

```
src/main/       Node side: window, IPC, spawn FasTag, TSV preview
src/preload/    the single allowlisted context bridge (window.fastag)
src/renderer/   React UI (form, log, results)
resources/      bundled FasTag (dev: symlinks; release: real files)
```
