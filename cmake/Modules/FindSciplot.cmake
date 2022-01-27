# Find sciplot, A modern C++ scientific plotting library powered by gnuplot
#
# This module will define the following variables:
#  Sciplot_FOUND - true if Sciplot was found
#  SCIPLOT_INCLUDE_DIR - Location of the sciplot includes
#  SCIPLOT_VERSION - Full version string
#
# This module will also create the sciplot::sciplot target.

# Find include directories
find_path(SCIPLOT_SCR_DIR Utils.hpp
        PATHS
        ${CMAKE_SOURCE_DIR}/thirdparty/sciplot/sciplot/
        )
set(SCIPLOT_INCLUDE_DIR "${SCIPLOT_SCR_DIR}/..")

# Extract the version
if (SCIPLOT_INCLUDE_DIR)
    file(READ "${SCIPLOT_INCLUDE_DIR}/CMakeLists.txt" _SCIPLOT_CMAKE_FILE)
    string(REGEX MATCH "([0-9]+).([0-9]+).([0-9]+)" SCIPLOT_VERSION "${_SCIPLOT_CMAKE_FILE}")
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sciplot
        REQUIRED_VARS SCIPLOT_INCLUDE_DIR
        VERSION_VAR SCIPLOT_VERSION
        HANDLE_COMPONENTS
        )

include(FeatureSummary)
set_package_properties(Sciplot PROPERTIES
        URL
        https://github.com/sciplot/sciplot.git
        DESCRIPTION
        "A modern C++ scientific plotting library powered by gnuplot"
        )

if (Sciplot_FOUND AND NOT TARGET sciplot::sciplot)
    add_library(sciplot::sciplot INTERFACE IMPORTED)
    target_include_directories(sciplot::sciplot INTERFACE ${SCIPLOT_INCLUDE_DIR})
endif()