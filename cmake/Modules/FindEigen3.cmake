# Find Eigen, the expression template matrix library.
#
# This module will define the following variables:
#  EIGEN3_FOUND - true if Eigen was found on the system
#  EIGEN3_INCLUDE_DIRS - Location of the Eigen includes
#  EIGEN3_VERSION - Full version string
#
# This module will also create the Eigen3::Eigen target.

# Find the include directory
find_path(EIGEN3_HEADER LLT.h
        PATHS ${CMAKE_SOURCE_DIR}/thirdparty/eigen/Eigen/src/Cholesky/
        )
set(EIGEN3_INCLUDE_DIR "${EIGEN3_HEADER}/../../../")

# Extract the version
if(EIGEN3_INCLUDE_DIR)
    file(READ "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)
    string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen3_world_version_match "${_eigen3_version_header}")
    set(EIGEN_VERSION_MAJOR "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen3_major_version_match "${_eigen3_version_header}")
    set(EIGEN_VERSION_MINOR "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen3_minor_version_match "${_eigen3_version_header}")
    set(EIGEN_VERSION_PATCH "${CMAKE_MATCH_1}")

    # Create the full version
    set(EIGEN3_VERSION "${EIGEN_VERSION_MAJOR}.${EIGEN_VERSION_MINOR}.${EIGEN_VERSION_PATCH}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3
        REQUIRED_VARS EIGEN3_INCLUDE_DIR
        VERSION_VAR EIGEN3_VERSION
        FOUND_VAR EIGEN3_FOUND
        HANDLE_COMPONENTS
        )

include(FeatureSummary)
set_package_properties(EIGEN3 PROPERTIES
        URL "http://eigen.tuxfamily.org/"
        DESCRIPTION "A C++ template library for linear algebra: matrices,
                 vectors, numerical solvers, and related algorithms"
        )

if (EIGEN3_FOUND AND NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
endif()