# Find SUNDIALS, the SUite of Nonlinear and DIfferential/ALgebraic equation Solvers.
#
# This module will define the following variables:
#  SUNDIALS_FOUND - true if SUNDIALS was found on the system
#  SUNDIALS_INCLUDE_DIRS - Location of the SUNDIALS includes
#  SUNDIALS_LIBRARIES - Required libraries for all requested components
#  SUNDIALS_VERSION_MAJOR - Major version
#  SUNDIALS_VERSION_MINOR - Minro version
#  SUNDIALS_VERSION_PATCH - Patch level
#  SUNDIALS_VERSION - Full version string
#
# This module exports the target SUNDIALS::<component>.

set(SUNDIALS_INSTALL_DIR "${CMAKE_SOURCE_DIR}/thirdparty/sundials/install/")
set(SUNDIALS_INCLUDE_DIR "${SUNDIALS_INSTALL_DIR}/include/")

# extract version
file(READ "${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h" _SUNDIALS_VERSION_FILE)
string(REGEX REPLACE ".*#define SUNDIALS_VERSION \"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\1" SUNDIALS_VERSION_MAJOR "${_SUNDIALS_VERSION_FILE}")
string(REGEX REPLACE ".*#define SUNDIALS_VERSION \"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\2" SUNDIALS_VERSION_MINOR "${_SUNDIALS_VERSION_FILE}")
string(REGEX REPLACE ".*#define SUNDIALS_VERSION \"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\3" SUNDIALS_VERSION_PATCH "${_SUNDIALS_VERSION_FILE}")

# Create the full version
set(SUNDIALS_VERSION "${SUNDIALS_VERSION_MAJOR}.${SUNDIALS_VERSION_MINOR}.${SUNDIALS_VERSION_PATCH}")

# Library settings
set(SUNDIALS_WANT_COMPONENTS
    sundials_cvodes
    sundials_nvecserial
)
set(SUNDIALS_PREFER_STATIC_LIBRARIES ON)

# find the SUNDIALS libraries
foreach(LIB ${SUNDIALS_WANT_COMPONENTS})
    if (UNIX AND SUNDIALS_PREFER_STATIC_LIBRARIES)
        set(THIS_LIBRARY_SEARCH lib${LIB}.a ${LIB})
    else()
        set(THIS_LIBRARY_SEARCH ${LIB})
    endif()

    find_library(SUNDIALS_${LIB}_LIBRARY
        NAMES ${THIS_LIBRARY_SEARCH}
        PATHS
            ${SUNDIALS_INSTALL_DIR}
        PATH_SUFFIXES
            lib
        )

    set(SUNDIALS_${LIB}_FOUND FALSE)
    if (SUNDIALS_${LIB}_LIBRARY)
        list(APPEND SUNDIALS_LIBRARIES ${SUNDIALS_${LIB}_LIBRARY})
        set(SUNDIALS_${LIB}_FOUND TRUE)
    endif()
    mark_as_advanced(SUNDIALS_${LIB}_LIBRARY)
endforeach()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS
    REQUIRED_VARS SUNDIALS_LIBRARIES SUNDIALS_INCLUDE_DIR
    VERSION_VAR SUNDIALS_VERSION
    HANDLE_COMPONENTS
)

include(FeatureSummary)
set_package_properties(SUNDIALS PROPERTIES
        URL "https://computation.llnl.gov/projects/sundials"
        DESCRIPTION "SUNDIALS Suite of nonlinear and differential/algebraic equation solvers"
        )

if (SUNDIALS_FOUND)
    foreach(LIB ${SUNDIALS_WANT_COMPONENTS})
        if (SUNDIALS_${LIB}_LIBRARY)
            add_library(SUNDIALS::${LIB} UNKNOWN IMPORTED)
            set_target_properties(SUNDIALS::${LIB} PROPERTIES
                    INTERFACE_INCLUDE_DIRECTORIES "${SUNDIALS_INCLUDE_DIR}"
                    IMPORTED_LOCATION ${SUNDIALS_${LIB}_LIBRARY}
                    )
        endif()
    endforeach()
endif()

mark_as_advanced(
    SUNDIALS_LIBRARIES
    SUNDIALS_INCLUDE_DIR
    SUNDIALS_INCLUDE_DIRS
)
