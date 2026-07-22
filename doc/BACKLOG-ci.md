# CI history: Linux, macOS and Windows

All three platforms are green: Linux x64/arm64 and macOS x64/arm64 gate `main`
in `ci.yml`; Windows x64 gates `main` in its own `windows.yml`, kept separate
purely because building OpenMS from source costs ~100 minutes against seconds
for the conda path everywhere else, not because it is in doubt. This records
what each one needed, so none of it has to be rediscovered.

**Why a slow-to-schedule or previously-red platform stays out of a shared
matrix until proven**: a job that has never passed teaches the habit of
ignoring red. Better to add a target once it means something than to carry
several where some are permanently amber.

## What Linux needed, and why it matters to the others

Six issues, in the order CI found them. Numbers 1-4 are platform-independent and
are **already fixed** in the workflow, so macOS and Windows inherit them:

| # | issue | fix | scope |
|---|---|---|---|
| 1 | Eigen 3.4 rejects the CMake version *range* `3.4.0...<6` that `OpenMSConfig.cmake` requests | pass `cmake/eigen3-range-shim` | all |
| 2 | `OpenMS/config.h` includes `boost/current_function.hpp`, but bioconda's `openms` does not depend on Boost | add `libboost-headers` | all |
| 3 | `FileTypes::MZPEAK` referenced outside its `#ifdef` — **a FasTag bug** | gate it; detect the enum, not just the header | all |
| 4 | `find_dependency(Qt6)` against an incomplete Qt6 | add `qt6-main` | all |
| 5 | conda's ARM Qt6 is marked cross-compiled and demands a host Qt | `-DQT_HOST_PATH=$CONDA_PREFIX` | ARM |
| 6 | contrib shells out to `cl.exe` and needs the MSVC developer env | `ilammy/msvc-dev-cmd` | Windows |

Not one was visible on any machine that had built FasTag before — laptop and HPC
both had Eigen worked around by hand and Boost/Qt supplied by contrib. These are
what a new user hits on a clean machine.

## macOS — BOTH GREEN

Resolved. `linux-x64`, `linux-arm64`, `macos-arm64` and `macos-x64` all passed in
run 29891819334, tests included. The earlier diagnosis in this file was wrong in
a way worth recording.

**macos-arm64 passes.** Not a prediction — it completed successfully in run
29857387912, in the very matrix that was then removed. bioconda's `osx-arm64`
OpenMS 3.5.0 works with the conda path unchanged; `QT_HOST_PATH` (issue 5) was
all it needed.

**macos-x64 was never a scheduling delay. `macos-13` is retired.** This file
previously said the runner was "scarce" and that turning off `cancel-in-progress`
"should let it schedule". It did not, because no runner was ever going to pick it
up — and turning off cancellation converted an impossible job into an **outage**:
that run sat `queued` holding the `refs/heads/main` concurrency group, and every
push to main was blocked behind it for about a day. It took cancelling the run by
hand to clear. The v0.9.0 tag build passed throughout only because tags get their
own concurrency group, which made main look merely slow rather than stuck.

`macos-15-intel`, GitHub's replacement Intel label, **schedules and passes**. The
entire "scarce runner" theory was wrong: one label change fixed what six runs and
a concurrency rewrite could not. The lesson is not about macOS — it is that
*"queued" was read as a capacity symptom for a day without anyone checking
whether the label still existed.*

Two changes make a repeat impossible to hide:

- `cancel-in-progress: true` — a wedge is now superseded by the next push instead
  of blocking it. The starvation this was meant to prevent had a different cause.
- `timeout-minutes: 40` — bounds a job that hangs. Note it does **not** bound a
  job waiting for a runner; GitHub starts that clock only once the job is picked
  up. The only defence against an unschedulable label is not to name one.

There is a local caveat CI does not exercise: the macOS **curl framework**
problem. CMake searches frameworks first, so a third-party
`/Library/Frameworks/libcurl.framework` shadows the SDK libcurl that
`libOpenMS.dylib` links, and the binary dies at startup. `CMakeLists.txt`
already pins curl to the SDK; a clean runner has no such framework, so a
regression there would pass CI and fail on a developer's machine.

## Windows — GREEN, in `windows.yml`

Resolved. `windows-x64-openms` and `windows-x64` both passed in run
29932193776 — configure, build, ctest (7/7) and the smoke test all green, no
`continue-on-error` left on either job.

**Two jobs, not one — discovered by hitting it, not designed in.**
`actions/cache`'s save is a post-hook gated on the *job* having succeeded, and
that gate is not overridden by `continue-on-error` at the job level (which
only changes how the workflow reports the job's outcome, not the runtime
`success()` a post-hook's own condition evaluates). The first attempt cached
contrib and the OpenMS install in the same job as FasTag's own
configure/build/test; every time FasTag's side failed, *both* caches silently
failed to save, so the next push repaid the full ~100-minute build to test a
one-line CMake flag. Split into `windows-x64-openms` (contrib + OpenMS
install, nothing FasTag-specific, so its own success — and therefore its cache
saves — no longer depends on FasTag configuring cleanly) and `windows-x64`
(needs the first job, restores both caches, redoes only the cheap per-runner
installs: choco, Qt, conda). This is the same two-job pattern OpenMS's own CI
uses for contrib, now understood to be load-bearing rather than incidental.

**Real timings, not estimates**: contrib 37 min, OpenMS itself 63 min, both on
a cold cache. Warm-cache reruns (only the `windows-x64` job, choco/Qt/conda
reinstalled fresh each time) take a few minutes.

**What FasTag's own configure needed, once OpenMS itself built** — all found
by direct diagnosis on a warm cache, not guessed:

| # | issue | fix |
|---|---|---|
| 1 | `OpenMSConfig.cmake` re-runs `find_dependency(XercesC)` etc. at every downstream consumer's configure time, using ordinary `find_package` — not via `OPENMS_CONTRIB_LIBS`, which is a hint OpenMS's own CMakeLists reads only at its own configure time | add `contrib-build` to `CMAKE_PREFIX_PATH` for FasTag's configure too |
| 2 | OpenMS's Windows install splits one logical CMake package directory across two destinations, neither complete alone: `<prefix>/CMake/` has only `OpenMSConfig(Version).cmake`; `<prefix>/bin/cmake/OpenMS/` has `OpenMSTargets(.cmake\|-release.cmake)` and `Modules/` (`FindLIBSVM.cmake` etc.) | copy `OpenMSConfig(Version).cmake` **into** the Targets/Modules directory, not the reverse — `OpenMSTargets.cmake` computes `${_IMPORT_PREFIX}` by walking a fixed parent-directory count baked in at its *original* install depth, so moving it instead silently corrupts every DLL path it exports |
| 3 | Test/benchmark binaries linking OpenMS failed at run time, `0xc0000135` (`STATUS_DLL_NOT_FOUND`) — Windows has no RPATH/RUNPATH; the loader only searches the `.exe`'s own directory and `PATH` | add OpenMS's, Qt's, contrib's and conda's `bin`/`Library/bin` directories to `PATH` before `ctest` |

Two more issues cost real iterations before being understood correctly:
`xargs dirname` silently corrupts a Windows path (backslash is xargs's own
quote-escape character), and `$GITHUB_WORKSPACE`/`QT_ROOT_DIR` are
backslash-separated in a way that CMake's own generated `try_compile` scratch
files reject outright (`Invalid character escape '\a'`) even though the
filesystem tolerates the mixed separators fine — both fixed by never piping a
Windows path through `xargs`, and by normalizing every workspace/Qt path to
forward slashes once, up front.

## Release binaries: built, not signed

Every tag build attaches its own binary to the GitHub Release automatically
(`gh release upload ... --clobber`, one per platform, racing harmlessly
across two workflows and five jobs — whichever finishes first creates the
release itself if it doesn't exist yet). Filenames carry no version
(`FasTag-<platform>.tar.gz`/`.zip`) so the README's
`/releases/latest/download/<name>` links stay valid across every future
release without editing.

**Signing scaffolding exists in both workflows (`ci.yml` for macOS,
`windows.yml` for Windows) but is inactive** — neither
`MACOS_CERTIFICATE_BASE64` nor `SIGNPATH_API_TOKEN` exists as a repo secret,
so every tag build takes the "not configured" branch and ships an unsigned
binary with a `::warning::`, exactly as before this scaffolding was added.
Adapted from a working reference (BALL/BALLView's `release.yml`), not
invented — but that project ships a Qt `.app` bundle and an NSIS installer;
FasTag ships one bare CLI executable, so the macOS path has no
`.app`/`.dmg`/entitlements/Sparkle machinery, and — importantly — no
`stapler staple` step, because stapling only works on bundles, installers
and disk images, never a lone Mach-O binary. A notarized standalone binary is
instead verified by Gatekeeper via an online ticket lookup on first run.

**To enable macOS signing**, add 7 repo secrets: `MACOS_CERTIFICATE_BASE64`
(a Developer ID Application cert+key, `.p12` export, base64-encoded),
`MACOS_CERTIFICATE_PASSWORD` (that `.p12`'s password), `MACOS_KEYCHAIN_PASSWORD`
(any password for the throwaway CI keychain), `MACOS_SIGNING_IDENTITY`
(`Developer ID Application: Name (TEAMID)`, read via
`security find-identity -p codesigning -v`), `MACOS_APPLE_ID`, `MACOS_TEAM_ID`,
`MACOS_NOTARY_PASSWORD` (an app-specific password for `notarytool`). Requires
an active Apple Developer Program membership. A Developer ID Application
certificate is scoped to the **Team**, not to one app — if the same Apple
Developer account already has one (for another project), it signs FasTag too;
no need for a second certificate.

**To enable Windows signing**, apply to
[SignPath Foundation](https://signpath.org/)'s free open-source signing
program for FasTag specifically (enrollment is per-project — a certificate
issued for a different project's slug does not cover this one), then add 2
repo secrets: `SIGNPATH_API_TOKEN` and `SIGNPATH_ORG_ID`. The private key
never leaves SignPath's HSM; there is no `.pfx` file to manage or lose.
Also requires, in SignPath's own console (not a repo file), an
artifact-configuration matching what `windows.yml` requests (slug `initial`
by default — edit the workflow if enrollment assigns a different one):
```xml
<artifact-configuration xmlns="http://signpath.io/artifact-configuration/v1">
  <zip-file><pe-file path="*.exe"><authenticode-sign /></pe-file></zip-file>
</artifact-configuration>
```
The release signing policy should keep manual approval on — `wait-for-completion:
true` blocks the job until a human clicks Approve in SignPath's UI, a
deliberate anti-supply-chain control, not a bug to route around.

**Two gotchas already avoided, not rediscovered**: `secrets.X != ''` used
directly in a step or job `if:` makes the **whole workflow** invalid — GitHub
emits a 0-job `startup_failure` on every push, tag included. Both signing
paths instead map the secret into a step `env:` (legal), and every downstream
step gates on a plain step *output* instead. And `actions/upload-artifact`
always zips its payload, so an artifact-configuration expecting a bare `.exe`
fails with "file does not correspond to the specified file type" — the
`<zip-file>` wrapper above is required, not optional.

## Anything else worth knowing

- `cancel-in-progress: true` in both workflows, for the reason given above: a
  wedge is superseded by the next push rather than blocking it.
- Every job verifies `OpenMSConfig.cmake` exists before configuring, because a
  package with tools but no development files otherwise fails deep inside CMake
  with a misleading message. bioconda's `openms` **does** ship it — that question
  is settled.
- The conda OpenMS ships `MzPeakFile.h` but its `FileTypes` enum has no `MZPEAK`,
  so CI builds **without** mzPeak support. The mzPeak read path is therefore
  exercised only locally, and nothing in CI covers it.
