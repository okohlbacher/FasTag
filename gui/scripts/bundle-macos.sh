#!/usr/bin/env bash
# Assemble the self-contained macOS resource tree the Tauri app bundles:
#
#   resources/fastag/
#     bin/FASTag              the CLI, install names rewritten to @executable_path/../lib
#     lib/*.dylib             its full dylib closure, with exactly ONE libomp
#     share/OpenMS/           OpenMS data (OPENMS_DATA_PATH)
#     share/FASTag/taxonomy/  the k-mer index + NCBI dumps (FASTAG_TAXONOMY_DIR)
#
# Bundling libomp exactly once is what fixes `OMP: Error #15 (libomp already
# initialized)` — verified: the bundled binary runs species detection in a clean
# env with no KMP_DUPLICATE_LIB_OK. Tauri's resource bundler skips symlinked
# dirs, so everything here must be real files.
#
# Inputs (env, with local-dev defaults):
#   FASTAG_BIN     the built CLI                     [build-rel/FASTag]
#   OPENMS_LIB     dylibbundler search path (OpenMS) [OpenMS-mzpeak-build/build/lib]
#   DEPS_LIB       dylibbundler search path (deps)   [fastag-mm/envs/omsbuild/lib]  (CI: $CONDA_PREFIX/lib)
#   OPENMS_SHARE   OpenMS data dir                   [OpenMS-mzpeak-build/share/OpenMS]
#   TAXONOMY_DIR   share/FASTag/taxonomy             [share/FASTag/taxonomy]
#   OUT            output resource root              [gui/src-tauri/resources/fastag]
set -euo pipefail

repo="$(cd "$(dirname "$0")/../.." && pwd)"
FASTAG_BIN="${FASTAG_BIN:-$repo/build-rel/FASTag}"
OPENMS_LIB="${OPENMS_LIB:-$repo/../OpenMS-mzpeak-build/build/lib}"
DEPS_LIB="${DEPS_LIB:-$repo/../fastag-mm/envs/omsbuild/lib}"
OPENMS_SHARE="${OPENMS_SHARE:-$repo/../OpenMS-mzpeak-build/share/OpenMS}"
TAXONOMY_DIR="${TAXONOMY_DIR:-$repo/share/FASTag/taxonomy}"
OUT="${OUT:-$repo/gui/src-tauri/resources/fastag}"

command -v dylibbundler >/dev/null || { echo "need dylibbundler (brew install dylibbundler)"; exit 1; }
[ -x "$FASTAG_BIN" ] || { echo "no CLI at $FASTAG_BIN"; exit 1; }

echo ">> staging binary + dylib closure"
rm -rf "$OUT"; mkdir -p "$OUT/bin" "$OUT/lib" "$OUT/share/FASTag"
cp "$FASTAG_BIN" "$OUT/bin/FASTag"

# -od overwrite, -b bundle, -p install-name prefix, -s extra search paths for @rpath deps.
dylibbundler -od -b \
  -x "$OUT/bin/FASTag" \
  -d "$OUT/lib" \
  -p '@executable_path/../lib/' \
  -s "$OPENMS_LIB" \
  -s "$DEPS_LIB" \
  </dev/null

# dylibbundler rewrites each pre-existing LC_RPATH to its -p prefix individually,
# so a Mach-O that had several rpaths ends up with duplicate '@executable_path/../lib/'
# entries — and a duplicate LC_RPATH is FATAL on arm64 (dyld refuses to launch).
# Collapse duplicates on every Mach-O, then ad-hoc re-sign (install_name_tool
# invalidates the signature; real Developer ID signing downstream overwrites it).
dedup_rpath() {
  local f="$1" dup
  for dup in $(otool -l "$f" | awk '/LC_RPATH/{g=1;next} g&&/path /{print $2; g=0}' | sort | uniq -d); do
    while [ "$(otool -l "$f" | awk '/LC_RPATH/{g=1;next} g&&/path /{print $2; g=0}' | grep -cx "$dup")" -gt 1 ]; do
      install_name_tool -delete_rpath "$dup" "$f"
    done
  done
  codesign --force --sign - "$f" 2>/dev/null || true
}
echo ">> dedup LC_RPATH on $(ls "$OUT"/lib/*.dylib | wc -l | tr -d ' ') libs + the binary"
dedup_rpath "$OUT/bin/FASTag"
for so in "$OUT"/lib/*.dylib; do [ -L "$so" ] && continue; dedup_rpath "$so"; done

echo ">> copying OpenMS data + taxonomy"
cp -R "$OPENMS_SHARE" "$OUT/share/OpenMS"
cp -R "$TAXONOMY_DIR" "$OUT/share/FASTag/taxonomy"

libomp=$(ls "$OUT"/lib | grep -c '^libomp\.dylib$' || true)
echo ">> done: $(du -sh "$OUT" | cut -f1) at $OUT (libomp copies: $libomp — must be 1)"
[ "$libomp" = "1" ] || { echo "!! expected exactly one libomp"; exit 1; }
