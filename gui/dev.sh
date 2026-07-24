#!/usr/bin/env bash
# Dev launcher for the FasTag GUI.
#
# In a real release the bundled FasTag is self-contained (dylibbundler fixes
# every @rpath and it links a single OpenMP runtime), so NONE of the env below
# is needed. But the *local dev* binaries are raw build-tree outputs whose
# runtime libs live in the OpenMS build tree + Homebrew, so we put those on the
# loader path here.
#
# resources/fastag/{bin/FasTag,share/OpenMS} are symlinks (see README); the GUI
# resolves them as the "bundled" binary. Point them at whichever build you want:
#   build-rel     -> mzPeak + latest features + progress (current main vs patched
#                    OpenMS) — what the release ships. Needs the env below.
#   build-main    -> latest features + progress, mzML only, no extra env needed.
set -euo pipefail
cd "$(dirname "$0")"

# OpenMS build-tree libs + Homebrew (Arrow 25, libsvm, xerces, libomp).
export DYLD_FALLBACK_LIBRARY_PATH="${DYLD_FALLBACK_LIBRARY_PATH:-}:/Users/kohlbach/Claude/OpenMS-mzpeak-build/build/lib:/opt/homebrew/lib:/opt/homebrew/opt/libomp/lib:/opt/homebrew/opt/libsvm/lib:/opt/homebrew/opt/xerces-c/lib"

# build-rel links Homebrew libomp while OpenMS/Arrow bring their own — two copies
# in one process. Harmless here; a single-toolchain release build won't hit it.
# ponytail: dev-only workaround, documented; the release build avoids it entirely.
export KMP_DUPLICATE_LIB_OK=TRUE

exec npm run dev
