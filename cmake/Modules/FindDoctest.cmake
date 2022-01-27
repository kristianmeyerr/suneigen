# Find Doctest, The fastest feature-rich C++11/14/17/20 single-header testing framework
#
# This module will define the following variables:
#  DOCTEST_FOUND - true if autodiff was found
#  DOCTEST_INCLUDE_DIR - Location of the autodiff includes
#  DOCTEST_VERSION - Full version string
#
# This module will also create the doctest::doctest target.

# Find the include directory
find_path(DOCTEST_HEADER doctest.h
    PATHS
        ${CMAKE_SOURCE_DIR}/thirdparty/doctest/doctest/
)
set(DOCTEST_INCLUDE_DIR "${DOCTEST_HEADER}/../")

# Extract the version
if (DOCTEST_INCLUDE_DIR)
    file(READ "${DOCTEST_INCLUDE_DIR}/scripts/version.txt" _DOCTEST_VERSION_FILE)
    string(REGEX MATCH "([0-9]+).([0-9]+).([0-9]+)" DOCTEST_VERSION "${_DOCTEST_VERSION_FILE}")
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Doctest
    REQUIRED_VARS DOCTEST_INCLUDE_DIR
    VERSION_VAR DOCTEST_VERSION
    HANDLE_COMPONENTS
)

include(FeatureSummary)
set_package_properties(Doctest PROPERTIES
    URL "https://github.com/doctest/doctest.git"
    DESCRIPTION "The fastest feature-rich C++11/14/17/20 single-header testing framework"
)

if (Doctest_FOUND AND NOT TARGET doctest::doctest)
    add_library(doctest::doctest INTERFACE IMPORTED)
    target_include_directories(doctest::doctest INTERFACE ${DOCTEST_INCLUDE_DIR})
endif()

mark_as_advanced(DOCTEST_INCLUDE_DIR)