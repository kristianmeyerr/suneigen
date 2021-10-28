# Find Eigen, the expression template matrix library.
#
# This module will define the following variables:
#  EIGEN_FOUND - true if Eigen was found on the system
#  EIGEN_INCLUDE_DIRS - Location of the Eigen includes
#  EIGEN_VERSION_MAJOR - Major version
#  EIGEN_VERSION_MINOR - Minor version
#  EIGEN_VERSION_PATCH - Patch level
#  EIGEN_VERSION - Full version string
#
# This module will also create the Eigen3::Eigen target.

# Set the include dir
set(EIGEN_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/thirdparty/eigen3")

# Extract the version
file(READ "${EIGEN_INCLUDE_DIRS}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)
string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen3_world_version_match "${_eigen3_version_header}")
set(EIGEN_VERSION_MAJOR "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen3_major_version_match "${_eigen3_version_header}")
set(EIGEN_VERSION_MINOR "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen3_minor_version_match "${_eigen3_version_header}")
set(EIGEN_VERSION_PATCH "${CMAKE_MATCH_1}")

# Create the full version
set(EIGEN_VERSION "${EIGEN_VERSION_MAJOR}.${EIGEN_VERSION_MINOR}.${EIGEN_VERSION_PATCH}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3
    REQUIRED_VARS
        EIGEN_INCLUDE_DIRS
        EIGEN_VERSION
)

include(FeatureSummary)
set_package_properties(EIGEN3 PROPERTIES
    URL
        "http://eigen.tuxfamily.org/"
    DESCRIPTION
        "A C++ template library for linear algebra: matrices, vectors,
                 numerical solvers, and related algorithms"
        )

if (EIGEN3_FOUND AND NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN_INCLUDE_DIRS})
endif()
