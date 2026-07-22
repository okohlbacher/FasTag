# CI backlog: macOS and Windows

Linux x64 and arm64 are green and gate `main`. macOS and Windows are deferred,
but were taken far enough that most of the work is diagnosis rather than
discovery. This records what is known so it does not have to be rediscovered.

**Why they are not in the matrix**: a job that has never passed teaches the habit
of ignoring red. Better to have two targets that mean something than six where
four are permanently amber.

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

The local macOS **curl framework** caveat still stands and CI cannot catch it: a
third-party `/Library/Frameworks/libcurl.framework` shadows the SDK libcurl that
`libOpenMS.dylib` links, and the binary dies at startup. `CMakeLists.txt` pins
curl to the SDK; a clean runner has no such framework, so a regression there would
pass CI and fail on a developer's machine.

There is a local caveat the CI does not exercise: the macOS **curl framework**
problem. CMake searches frameworks first, so a third-party
`/Library/Frameworks/libcurl.framework` shadows the SDK libcurl that
`libOpenMS.dylib` links, and the binary dies at startup. `CMakeLists.txt` already
pins curl to the SDK; a clean runner has no such framework, so CI would not catch
a regression there.

## Windows

**Status: further away, and the obstacle is OpenMS, not FasTag.**

There is no OpenMS conda package for Windows and the official installer carries
**no SDK** — only end-user executables. So Windows must build OpenMS *and*
contrib from source: 30-60 minutes even cached, against a couple of minutes for
the conda path.

The last failure was issue 6 above (MSVC environment), fixed but never
re-verified. The next unknowns, in the order they will appear:

1. Does contrib build cleanly on a hosted runner within the time limit?
2. Does OpenMS build with `WITH_GUI=OFF`?
3. Does FasTag compile under MSVC? `__uint128_t` was ported to `Kmer128` for
   exactly this, and that was the only known blocker — but it has never been
   compiled by MSVC, so others may surface. OpenMP under MSVC is the next most
   likely: it supports only OpenMP 2.0, and the parallel loop uses `SignedSize`
   induction variables partly for that reason.

**windows-arm64 is likely not achievable.** OpenMS has no ARM Windows support I
am aware of, and contrib may not configure at all. Attempt it only after
windows-x64 is green, and be prepared to declare it unsupported rather than leave
it permanently red.

A cached OpenMS build is essential — key on the OpenMS tag so it is paid once per
revision. The workflow that was removed did this; recover it from git history at
`4952e2d` rather than rewriting it.

## Anything else worth knowing

- `cancel-in-progress: false` is deliberate. See above.
- Every job verifies `OpenMSConfig.cmake` exists before configuring, because a
  package with tools but no development files otherwise fails deep inside CMake
  with a misleading message. bioconda's `openms` **does** ship it — that question
  is settled.
- The conda OpenMS ships `MzPeakFile.h` but its `FileTypes` enum has no `MZPEAK`,
  so CI builds **without** mzPeak support. The mzPeak read path is therefore
  exercised only locally, and nothing in CI covers it.
