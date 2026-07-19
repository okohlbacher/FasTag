cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(FasTag VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
  # OpenMS 3.6 (currently the develop branch; 3.5 lacks APIs used here).
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
