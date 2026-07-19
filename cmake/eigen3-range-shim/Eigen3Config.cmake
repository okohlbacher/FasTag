# Shim for Eigen 3.4's config-version file, which rejects CMake version RANGES.
#
# OpenMSConfig.cmake asks for `Eigen3 3.4.0...<6`. Eigen 3.4's shipped
# Eigen3ConfigVersion.cmake enforces a same-major-version policy: it marks any
# range whose maximum exceeds the next major version (4) as incompatible, so a
# perfectly good Eigen 3.4.0 is refused and configuration fails before OpenMS's
# own fallback can run.
#
# This shim simply declares the imported target against the real headers.
# Use it with:  -DEigen3_DIR=<this directory> -DEIGEN3_INCLUDE_DIR=<eigen headers>
if(NOT DEFINED EIGEN3_INCLUDE_DIR OR NOT EXISTS "${EIGEN3_INCLUDE_DIR}/signature_of_eigen3_matrix_library")
  find_path(EIGEN3_INCLUDE_DIR signature_of_eigen3_matrix_library
            PATH_SUFFIXES eigen3 include/eigen3)
endif()

if(NOT EIGEN3_INCLUDE_DIR)
  set(Eigen3_FOUND FALSE)
  return()
endif()

if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()

set(EIGEN3_INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")
set(Eigen3_FOUND TRUE)
