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

# OpenMS build-tree libs + Homebrew (Arrow, libsvm, xerces, libomp).
#
# OPENMS_LIB_DIR points at the lib/ of the OpenMS build or install tree the
# symlinked binary was linked against. Derived from the binary's own OpenMS data
# symlink when not set, so this works on any checkout rather than one machine.
if [ -z "${OPENMS_LIB_DIR:-}" ] && [ -L resources/fastag/share/OpenMS ]; then
  # .../<tree>/share/OpenMS -> <tree>. An INSTALL tree keeps its libs in
  # <tree>/lib, a BUILD tree in <tree>/build/lib; take whichever exists.
  _t="$(cd "$(dirname "$(readlink resources/fastag/share/OpenMS)")/.." 2>/dev/null && pwd)"
  for _c in "$_t/lib" "$_t/build/lib"; do
    [ -d "$_c" ] && OPENMS_LIB_DIR="$_c"
  done
fi
BREW_PREFIX="$(brew --prefix 2>/dev/null || echo /opt/homebrew)"
export DYLD_FALLBACK_LIBRARY_PATH="${DYLD_FALLBACK_LIBRARY_PATH:-}:${OPENMS_LIB_DIR:-}:$BREW_PREFIX/lib:$BREW_PREFIX/opt/libomp/lib:$BREW_PREFIX/opt/libsvm/lib:$BREW_PREFIX/opt/xerces-c/lib"

# build-rel links Homebrew libomp while OpenMS/Arrow bring their own — two copies
# in one process. Harmless here; a single-toolchain release build won't hit it.
# ponytail: dev-only workaround, documented; the release build avoids it entirely.
export KMP_DUPLICATE_LIB_OK=TRUE

exec npm run dev
