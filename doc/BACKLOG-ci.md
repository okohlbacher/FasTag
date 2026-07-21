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

## macOS

**Status: unverified, and closer than it looks.** bioconda ships OpenMS 3.5.0 for
both `osx-64` and `osx-arm64`, so the conda path applies unchanged. macos-arm64
failed on issue 5, which is now fixed — it plausibly passes as-is.

The real obstacle is **macos-x64 was never once built**. `macos-13` is the last
Intel runner and is scarce; across six consecutive runs it sat queued while each
new push cancelled the run out from under it. `cancel-in-progress` is now off,
which should let it schedule.

Note it read as *"queued"* rather than *"unverified"* — an untested platform
looking merely slow in an otherwise-green matrix is the more dangerous of the
two states, and is the reason it is called out here rather than left implicit.

To re-enable, restore to the matrix:

```yaml
- { name: macos-x64,   runner: macos-13 }
- { name: macos-arm64, runner: macos-14 }
```

Expect macos-arm64 to pass. Watch macos-13 for scheduling rather than build
failure — and if Intel runners are retired, drop that target rather than pretend
it is covered.

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
