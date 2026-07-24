cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(FasTag VERSION 0.17.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE AND NOT DEFINED CURL_LIBRARY)
  # Resolve curl to the macOS SDK, as OpenMS itself does.
  #
  # OpenMSConfig.cmake pulls in CURL transitively, and CMake searches frameworks
  # first on macOS, so a third-party /Library/Frameworks/libcurl.framework
  # shadows the SDK libcurl that libOpenMS.dylib actually links. The result is a
  # binary carrying an @rpath reference to a framework on no rpath, which dies at
  # startup. Note the SDK ships a .tbd stub, not a .dylib; pointing at
  # /usr/lib/libcurl.4.dylib instead leaves CURL::libcurl undefined.
  execute_process(COMMAND xcrun --show-sdk-path
                  OUTPUT_VARIABLE _sdk OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(_sdk AND EXISTS "${_sdk}/usr/lib/libcurl.tbd")
    set(CURL_LIBRARY "${_sdk}/usr/lib/libcurl.tbd" CACHE FILEPATH "curl (macOS SDK)")
    set(CURL_INCLUDE_DIR "${_sdk}/usr/include" CACHE PATH "curl headers (macOS SDK)")
  endif()
endif()

# Version is checked after the fact, not in find_package.
#
# OpenMS ships a SameMinorVersion config-version file, so find_package(OpenMS 3.5)
# rejects 3.6 outright -- a MINIMUM version cannot be expressed through
# find_package at all. Ask for any OpenMS and enforce the floor ourselves.
find_package(OpenMS REQUIRED)
if(OpenMS_VERSION VERSION_LESS 3.5)
  message(FATAL_ERROR "FasTag needs OpenMS >= 3.5 (found ${OpenMS_VERSION}).")
endif()
message(STATUS "FasTag: using OpenMS ${OpenMS_VERSION} (${OpenMS_DIR})")

# Report at configure time which reader path this OpenMS will allow.
#
# Stock OpenMS parses every spectrum's metadata serially in openFile() before
# tagging starts -- ~55 s of a 60 s run on a 5 GB file, and the reason wall time
# stops improving past 16 threads. FasTag skips that pass when the OpenMS it
# links can report MS level and precursor from the spectrum block itself, which
# needs a patch to MzMLSpectrumDecoder (see doc/OPENMS-FAST-READER.md).
#
# FasTag detects this at RUNTIME too and falls back with a warning, so an
# unpatched OpenMS is correct, merely slow. This check exists so that is known
# when the build is configured rather than discovered from a log much later.
set(_fastag_fast_reader OFF)
get_target_property(_openms_inc OpenMS INTERFACE_INCLUDE_DIRECTORIES)
if(_openms_inc)
  foreach(_d IN LISTS _openms_inc)
    string(GENEX_STRIP "${_d}" _d)
    if(_d AND EXISTS "${_d}/OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h")
      file(READ "${_d}/OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h" _hdr)
      string(FIND "${_hdr}" "handleSpectrumMetadata_" _found)
      if(NOT _found EQUAL -1)
        set(_fastag_fast_reader ON)
      endif()
      break()
    endif()
  endforeach()
endif()

if(_fastag_fast_reader)
  message(STATUS "FasTag: fast reader available -- this OpenMS reports "
                 "per-spectrum metadata")
else()
  message(STATUS "FasTag: fast reader NOT available -- this OpenMS parses all "
                 "metadata up front. Correct, but expect roughly a minute of "
                 "extra single-threaded startup per large file, and no scaling "
                 "past ~16 threads. See doc/OPENMS-FAST-READER.md.")
endif()

# mzPeak input, if this OpenMS provides it.
#
# MzPeakFile is not in OpenMS mainline yet (github.com/OpenMS/mzpeak, and the
# OpenMS-side handler on a feature branch), so FasTag must build either way. When
# the header is present, -in also accepts .mzpeak; when it is not, the tool is
# identical minus that format.
#
# Detected by header presence rather than version, because there is no version to
# test against while the API is pre-1.0.
# Requires BOTH the MzPeakFile header and a MZPEAK member in the FileTypes enum.
#
# Those are not the same thing. bioconda's OpenMS ships MzPeakFile.h while its
# FileTypes enum has no MZPEAK, so detecting only the header enabled code that
# then failed to compile on every CI platform. Detect what the code actually
# uses.
set(_fastag_mzpeak OFF)
if(_openms_inc)
  foreach(_d IN LISTS _openms_inc)
    string(GENEX_STRIP "${_d}" _d)
    if(_d AND EXISTS "${_d}/OpenMS/FORMAT/MzPeakFile.h"
          AND EXISTS "${_d}/OpenMS/FORMAT/FileTypes.h")
      file(READ "${_d}/OpenMS/FORMAT/FileTypes.h" _ft)
      string(FIND "${_ft}" "MZPEAK" _has_enum)
      # store() as well as the enum, because -out_spectra now writes mzPeak too.
      # Same principle as the enum check above: detect every API the code calls,
      # not a proxy for them.
      file(READ "${_d}/OpenMS/FORMAT/MzPeakFile.h" _mp)
      string(FIND "${_mp}" "void store(" _has_store)
      if(NOT _has_enum EQUAL -1 AND NOT _has_store EQUAL -1)
        set(_fastag_mzpeak ON)
      endif()
      break()
    endif()
  endforeach()
endif()

if(_fastag_mzpeak)
  message(STATUS "FasTag: mzPeak input available (-in accepts .mzpeak)")
else()
  message(STATUS "FasTag: mzPeak input NOT available -- this OpenMS lacks "
                 "MzPeakFile or a MZPEAK entry in FileTypes. mzML input is "
                 "unaffected.")
endif()

# OpenMP, declared rather than inherited.
#
# The parallel loop is the tool's central performance claim, and it was working
# only because OpenMS happened to propagate the flags transitively. A clean build
# against an OpenMS that did not would compile silently -- the #pragmas become
# no-ops and the omp_* fallback stubs in FasTag.cpp return 1 -- producing a
# correct but SERIAL binary that still claims to scale. Ask explicitly, and say
# so when it is missing.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "FasTag: OpenMP ${OpenMP_CXX_VERSION} -- parallel over spectra")
else()
  message(WARNING "FasTag: no OpenMP. The build is correct but SERIAL; -threads "
                  "will have no effect and README.md's performance figures do "
                  "not apply.")
endif()

add_library(fastag_core STATIC src/FasTagger.cpp src/FastaFilter.cpp src/TagRecon.cpp src/TaxIndex.cpp src/TaxStats.cpp)
target_include_directories(fastag_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
target_link_libraries(fastag_core PUBLIC OpenMS)

add_executable(FasTag src/FasTag.cpp)
target_link_libraries(FasTag PRIVATE fastag_core)
if(OpenMP_CXX_FOUND)
  target_link_libraries(FasTag PRIVATE OpenMP::OpenMP_CXX)
endif()
if(_fastag_mzpeak)
  target_compile_definitions(FasTag PRIVATE FASTAG_HAVE_MZPEAK)
endif()
install(TARGETS FasTag RUNTIME DESTINATION bin)

# The taxonomy -species resolves by default, at <bin>/../share/FasTag/taxonomy.
#
# Only the pruned NCBI dumps (~30 KB) are installed here. The k-mer index is
# ~1 GB and platform-independent, so it ships once as the release asset
# FasTag-taxonomy-k7.tar.gz rather than inside every per-platform tarball -- see
# the README beside these files. FasTag names the directory it searched when the
# index is absent, so that is a clear error, never a silent empty report.
# reference-set.tsv is the manifest; it is installed too so a user can rebuild.
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/data/taxonomy/nodes.dmp"
              "${CMAKE_CURRENT_SOURCE_DIR}/data/taxonomy/names.dmp"
              "${CMAKE_CURRENT_SOURCE_DIR}/data/taxonomy/README.md"
              "${CMAKE_CURRENT_SOURCE_DIR}/data/taxonomy/reference-set.tsv"
        DESTINATION share/FasTag/taxonomy)

add_executable(buildtaxdb src/buildtaxdb.cpp)
target_link_libraries(buildtaxdb PRIVATE fastag_core)
install(TARGETS buildtaxdb RUNTIME DESTINATION bin)

add_executable(tagfeatures src/tagfeatures.cpp)
install(TARGETS tagfeatures RUNTIME DESTINATION bin)

option(FASTAG_BUILD_BENCH "Build the synthetic benchmark" ON)
if(FASTAG_BUILD_BENCH)
  add_executable(benchmark bench/benchmark.cpp)
  target_link_libraries(benchmark PRIVATE fastag_core)
endif()

option(FASTAG_BUILD_TESTS "Build tests" ON)
if(FASTAG_BUILD_TESTS)
  enable_testing()
  foreach(t ranksum_test fastag_test fasta_filter_test tagrecon_test taxindex_test taxstats_test sampler_test proforma_test)
    add_executable(${t} test/${t}.cpp)
    target_link_libraries(${t} PRIVATE fastag_core)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()

  add_executable(tagfeatures_test test/tagfeatures_test.cpp src/tagfeatures.cpp)
  target_include_directories(tagfeatures_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_compile_definitions(tagfeatures_test PRIVATE TAGFEATURES_NO_MAIN)
  add_test(NAME tagfeatures_test COMMAND tagfeatures_test)

  # One short run per acquisition class, gaps and deisotoping ON.
  #
  # These are smoke tests, not quality gates -- the benchmark exits non-zero only
  # if a profile produces no tags at all. That is deliberately weak, and it is
  # still the check that would have caught `-deisotope` throwing an uncaught
  # exception at the 0.3 Da tolerance the iontrap profile uses. Every profile
  # other than iontrap passed that day.
  if(FASTAG_BUILD_BENCH)
    foreach(prof astral ddapasef diatracer iontrap)
      add_test(NAME bench_${prof}
               COMMAND benchmark synth ${prof} 150 4 1 2 100)
    endforeach()
  endif()
endif()
