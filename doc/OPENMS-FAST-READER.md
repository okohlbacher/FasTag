# Building an OpenMS that lets FasTag use the fast reader

FasTag runs correctly against stock OpenMS. It just spends most of its wall time
waiting for the reader, and stops scaling past about 16 threads no matter how
many cores it has. On a 5.3 GB mzML with 632,677 spectra that is the difference
between **60 s and 9 s**.

This is optional. Skip it and everything still works — FasTag detects the
situation at configure time and again at runtime, and takes the slow path.

## Why

`OnDiscMSExperiment::openFile()` calls `loadMetaData_()`, which streams the whole
file to collect per-spectrum metadata *before any tagging starts*: measured at
~55 s of a 60 s run, single-threaded. That prologue is ~95% of the runtime, so
adding cores does nothing.

`openFile(file, skipMetaData=true)` removes it, but on stock OpenMS
`MzMLSpectrumDecoder::domParseSpectrum()` fills only the binary arrays and the
native ID. MS level is 0 and there is no precursor, so every spectrum fails the
MS2 test and the run produces nothing.

The patch fixes that at the source: `domParseSpectrum` already builds a DOM of
the entire `<spectrum>` element and reads three things from it, discarding MS
level, retention time and precursors. It now keeps them. No extra I/O, no second
parse — 101 lines of production code.

## Applying it

The patch is in this repo's handoff bundle, or as OpenMS branch
`feature/ondisc-metadata-from-block`.

```bash
cd <your-openms-checkout>
git apply /path/to/0001-ondisc-metadata-from-block.patch
```

Then build. **On macOS, pin curl first** — see below, or the build will link and
then die at startup.

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DWITH_GUI=OFF \
  -DCURL_LIBRARY_RELEASE="$(xcrun --show-sdk-path)/usr/lib/libcurl.tbd" \
  -DCURL_INCLUDE_DIR="$(xcrun --show-sdk-path)/usr/include"
cmake --build build --target OpenMS -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
```

Point FasTag at it:

```bash
cmake -S . -B build -DOpenMS_DIR=<your-openms-checkout>/build
```

Configure prints which path you get:

```
-- FasTag: fast reader available -- this OpenMS reports per-spectrum metadata
```

## The macOS curl trap

Worth its own section because the failure is late and misleading: the build
compiles and links cleanly, then **every binary dies at startup** with

```
dyld: Library not loaded: @rpath/libcurl.framework/Versions/Release-7.61.1/libcurl
  Referenced from: .../libOpenMS.dylib
```

CMake searches frameworks before libraries on macOS, so a third-party
`/Library/Frameworks/libcurl.framework` shadows the SDK's curl. libOpenMS ends up
with an `@rpath` reference to a framework that is on no rpath. Pinning
`CURL_LIBRARY_RELEASE` and `CURL_INCLUDE_DIR` to the SDK avoids it.

Note the SDK ships a `.tbd` stub, not a `.dylib`. Pointing at
`/usr/lib/libcurl.4.dylib` instead leaves `CURL::libcurl` undefined and fails
differently. FasTag's own `CMakeLists.txt` already does this for its build; the
OpenMS build has no such guard, which is why it must be passed explicitly.

Check with:

```bash
otool -L build/lib/libOpenMS.dylib | grep -i curl
# want:  /usr/lib/libcurl.4.dylib
# not:   @rpath/libcurl.framework/...
```

## Verifying the patch is actually live

Three checks, in increasing order of how much they prove. Do at least the second
— a patched *source tree* and a patched *linked library* are not the same thing,
and confusing them has already produced one wrong benchmark here.

```bash
# 1. is it in the library you will link?
nm -gU build/lib/libOpenMS.dylib | grep -c handleSpectrumMetadata_   # macOS
nm -DC build/lib/libOpenMS.so    | grep -c handleSpectrumMetadata_   # Linux
# expect 1

# 2. does FasTag take the fast path? Run on an INDEXED mzML and look for the
#    absence of this warning:
#      "This OpenMS does not report spectrum metadata from the per-spectrum
#       reader, so the fast path is unavailable"

# 3. does the output still match? Run the same data indexed and non-indexed;
#    the non-indexed input forces the full-load fallback.
cmp fast.tsv slow.tsv && echo identical
```

A file needs an index for any of this to matter. An `indexedmzML` wrapper is not
enough — it must contain an actual `<indexList>`, or `openFile` returns false and
FasTag loads the whole run regardless. This bites test fixtures in particular:
a before/after comparison can pass while never exercising the streaming path at
all.

```bash
FileConverter -in in.mzML -out indexed.mzML -write_scan_index true
```

## What it does not change

Only per-spectrum fields move. Run-level `ExperimentalSettings` — instrument
configuration, sample, and so on — still needs the full parse, so `getMetaData()`
returns null under `skipMetaData`. That is why FasTag's `-out_spectra` path keeps
the slow load: it needs those settings to write a valid mzML.

## Status

Not upstream. Do not open a pull request against OpenMS without the repository
owner's explicit permission.
