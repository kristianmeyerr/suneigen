cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(suneigen VERSION 0.0.1 LANGUAGES CXX)

# --------------------------------------------------------------------------------
#                         Options
# --------------------------------------------------------------------------------

# Set options via commandline (e.g. -DENABLE_LTO=[ON|OFF])
option(ENABLE_WARNINGS_SETTINGS "Allow target_set_warnings to add flags and defines.
                                 Set this to OFF if you want to provide your own warning parameters." ON)
option(ENABLE_LTO "Enable link time optimization" ON)

option(SUNEIGEN_BUILD_PYTHON "Enable the compilation of the python module." ON)
option(SUNEIGEN_BUILD_TESTS "Enable the compilation of the test files." ON)
option(SUNEIGEN_BUILD_APP "Enable the compilation of the app files." ON)
option(SUNEIGEN_USE_DOXYGEN "Add a doxygen target to generate the documentation." ON)

if(SUNEIGEN_BUILD_PYTHON)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# --------------------------------------------------------------------------------
#                         CMake includes
# --------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# Configure the project for testing with CTest
include(CTest)

# Include custom modules
include(cmake/ConfigSafeGuards.cmake)
include(cmake/LTO.cmake)
include(cmake/Warnings.cmake)

find_package(Doctest 2 REQUIRED)

find_package(Eigen3 3.4 REQUIRED)
find_package(SUNDIALS 5.8 REQUIRED COMPONENTS sundials_cvodes sundials_nvecserial)
set(GSL_LITE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/thirdparty/gsl/")

# Check for LTO support.
find_lto(CXX)

if (${SUNEIGEN_USE_DOXYGEN})
    find_package(Doxygen) # Triggers -Wpoison-system-directories
    if(DOXYGEN_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(doc
                ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile &> doxygen.log
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "${BoldMagenta}Generating API documentation with Doxygen (open ./html/index.html to view).${ColourReset}" VERBATIM
                )
    endif()
endif ()

# --------------------------------------------------------------------------------
#                            Build the library
# --------------------------------------------------------------------------------

set(LIBRARY_NAME suneigen)

set(LIB_SOURCES
        src/sunlinsol_superlu.cpp
        src/solver.cpp
        src/solver_cvodes.cpp
        src/misc.cpp
        src/vector.cpp
        src/exception.cpp
        src/model.cpp
        src/model_ode.cpp
        src/suneigen.cpp
        src/return_data.cpp
        src/forward_problem.cpp
        src/sundials_matrix_wrapper.cpp
        src/sundials_linsol_wrapper.cpp
        src/model_state.cpp
        src/abstract_model.cpp)

# Compile all sources into a library.
add_library(${LIBRARY_NAME} OBJECT ${LIB_SOURCES})

# Add an alias suneigen::suneigen to the target library suneigen
add_library(${LIBRARY_NAME}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

# Lib needs its header files, and users of the library must also see these (PUBLIC).
target_include_directories(${LIBRARY_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${GSL_LITE_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/thirdparty/sundials/src  # This is needed for cvodes/cvodes_impl.h
)

# Use the same visiliblity flag as in pybind11
target_compile_options(${LIBRARY_NAME}
    PRIVATE
        -fvisibility=hidden
)

# Install suneigen interface library
install(TARGETS ${LIBRARY_NAME} DESTINATION bin)

# Link libraries used.
target_link_libraries(${LIBRARY_NAME}
    PUBLIC
        Eigen3::Eigen
        SUNDIALS::sundials_cvodes
        SUNDIALS::sundials_nvecserial
        ${CMAKE_DL_LIBS}  # for dladdr function
)

# Set the compile options you want (change as needed).
target_set_warnings(${LIBRARY_NAME} ENABLE ALL AS_ERROR ALL DISABLE Annoying)
# target_compile_options(${LIBRARY_NAME} ... )  # For setting manually.

set_target_properties(${LIBRARY_NAME}
    PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

# --------------------------------------------------------------------------------
#                            Build Library tests
# --------------------------------------------------------------------------------

if(SUNEIGEN_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# --------------------------------------------------------------------------------
#                            Build python module
# --------------------------------------------------------------------------------

if(SUNEIGEN_BUILD_PYTHON)

    add_subdirectory(thirdparty/pybind11)

    set(PYTHON_MODULE_NAME suneigen4py)

    set(PYTHON_SOURCES
            python/main.cpp
            python/py_model_ode.cpp
            python/pycvodes.cpp
    )

    # Compile suneigen source files into a Python module
    pybind11_add_module(${PYTHON_MODULE_NAME} ${PYTHON_SOURCES})

    # Link autodiff against autodiff C++ library
    target_link_libraries(${PYTHON_MODULE_NAME}
        PRIVATE
            ${LIBRARY_NAME}
            Eigen3::Eigen
    )

    # SETUP_VERSION_INFO is defined by setup.py and passed into the C++ code as a
    # define (VERSION_INFO) here.
    target_compile_definitions(${PYTHON_MODULE_NAME}
        PRIVATE
            VERSION_INFO=${SETUP_VERSION_INFO}
    )

    # https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
    target_compile_options(${PYTHON_MODULE_NAME}
        PRIVATE
            -fvisibility=hidden
    )

    # Enable link-time-optimization if available for non-debug configurations
    target_enable_lto(${PYTHON_MODULE_NAME} optimized)

    set_target_properties(${PYTHON_MODULE_NAME}
        PROPERTIES
            CXX_STANDARD 20
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
    )
endif()

# --------------------------------------------------------------------------------
#                            Build the c++ application
# --------------------------------------------------------------------------------

if(SUNEIGEN_BUILD_APP)

    add_subdirectory(app/)

endif()
