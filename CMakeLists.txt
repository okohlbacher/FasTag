cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(FasTag VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE AND NOT DEFINED CURL_LIBRARY)
  # Resolve curl to the macOS SDK, as OpenMS itself does.
  #
  # OpenMSConfig.cmake pulls in CURL transitively. CMake searches frameworks
  # first on macOS, so a third-party /Library/Frameworks/libcurl.framework
  # shadows the SDK libcurl that libOpenMS.dylib actually links. The resulting
  # binary carries an @rpath reference to a framework on no rpath and dies at
  # startup with "Library not loaded: @rpath/libcurl.framework/...".
  #
  # Note the SDK ships a .tbd text stub rather than a .dylib; pointing at
  # /usr/lib/libcurl.4.dylib instead leaves CURL::libcurl undefined.
  execute_process(COMMAND xcrun --show-sdk-path
                  OUTPUT_VARIABLE _fastag_sdk
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(_fastag_sdk AND EXISTS "${_fastag_sdk}/usr/lib/libcurl.tbd")
    set(CURL_LIBRARY "${_fastag_sdk}/usr/lib/libcurl.tbd" CACHE FILEPATH "curl (macOS SDK)")
    set(CURL_INCLUDE_DIR "${_fastag_sdk}/usr/include" CACHE PATH "curl headers (macOS SDK)")
  endif()
endif()

option(FASTAG_BUILD_TOOL  "Build the FasTag TOPP tool (requires OpenMS)" ON)
option(FASTAG_BUILD_TESTS "Build the algorithm tests (no OpenMS needed)"  ON)

# The algorithm core is header-only and deliberately free of OpenMS (and of every
# other) dependency, so it can be tested without a full OpenMS build.
add_library(fastag_core INTERFACE)
target_include_directories(fastag_core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(fastag_core INTERFACE cxx_std_17)

if(FASTAG_BUILD_TOOL)
  # OpenMS 3.5 or newer.
  # Point CMake at your build with -DOpenMS_DIR=/path/to/OpenMS-build
  find_package(OpenMS REQUIRED)
  add_executable(FasTag src/FasTag.cpp)
  target_link_libraries(FasTag PRIVATE fastag_core OpenMS)
  install(TARGETS FasTag RUNTIME DESTINATION bin)
endif()

if(FASTAG_BUILD_TESTS)
  enable_testing()
  find_package(Threads REQUIRED)
  foreach(t ranksum_test fastag_test fasta_filter_test)
    add_executable(${t} test/${t}.cpp)
    target_link_libraries(${t} PRIVATE fastag_core Threads::Threads)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
endif()

install(DIRECTORY include/FasTag DESTINATION include)
