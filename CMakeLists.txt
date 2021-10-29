cmake_minimum_required(VERSION 3.15)

project(SUNEIGEN VERSION 0.0.1 LANGUAGES CXX)

# --------------------------------------------------------------------------------
#                         Options
# --------------------------------------------------------------------------------

# Set options via commandline (e.g. -DENABLE_LTO=[ON|OFF])
option(ENABLE_WARNINGS_SETTINGS "Allow target_set_warnings to add flags and defines.
                                 Set this to OFF if you want to provide your own warning parameters." ON)
option(ENABLE_LTO "Enable link time optimization" ON)
option(ENABLE_DOCTESTS "Include tests in the library. Setting this to OFF will remove all doctest related code.
                        Tests in tests/*.cpp will still be enabled." ON)
option(SUNEIGEN_BUILD_PYTHON "Enable the compilation of the python module." ON)
option(SUNEIGEN_BUILD_TESTS "Enable the compilation of the test files." ON)
option(SUNEIGEN_BUILD_APP "Enable the compilation of the app files." ON)

# --------------------------------------------------------------------------------
#                         CMake modules
# --------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Configure the project for testing with CTest
include(CTest)

# Include custom modules
include(ConfigSafeGuards)
include(LTO)
include(Misc)
include(Warnings)

# Check for LTO support.
find_lto(CXX)

# --------------------------------------------------------------------------------
#                         Thirdparty dependencies
# --------------------------------------------------------------------------------
set(t t)
include(Doctest)
find_package(Eigen3 3.4 REQUIRED)
find_package(SUNDIALS 5.8 REQUIRED COMPONENTS sundials_cvodes sundials_nvecserial)

# --------------------------------------------------------------------------------
#                            Build the library
# --------------------------------------------------------------------------------

set(LIBRARY_NAME suneigen)

set(LIB_SOURCES
        src/example.cpp
)

# Compile all sources into a library.
add_library(${LIBRARY_NAME} OBJECT ${LIB_SOURCES})

# Add an alias suneigen::suneigen to the target library suneigen
add_library(${LIBRARY_NAME}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

# Lib needs its header files, and users of the library must also see these (PUBLIC).
target_include_directories(${LIBRARY_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Install suneigen interface library
install(TARGETS ${LIBRARY_NAME} DESTINATION bin)

# Link libraries used.
target_link_libraries(${LIBRARY_NAME}
    PUBLIC
        doctest
        Eigen3::Eigen
        SUNDIALS::sundials_cvodes
        SUNDIALS::sundials_nvecserial
)

# Set the compile options you want (change as needed).
target_set_warnings(${LIBRARY_NAME} ENABLE ALL AS_ERROR ALL DISABLE Annoying)
# target_compile_options(${LIBRARY_NAME} ... )  # For setting manually.

set_target_properties(${LIBRARY_NAME}
    PROPERTIES
        CXX_STANDARD 17
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

    set(PYTHON_SOURCES
            python/main.cpp
    )

    # Compile suneigen source files into a Python module
    pybind11_add_module(suneigen4py ${PYTHON_SOURCES})

    # SETUP_VERSION_INFO is defined by setup.py and passed into the C++ code as a
    # define (VERSION_INFO) here.
    target_compile_definitions(suneigen4py PRIVATE VERSION_INFO=${SETUP_VERSION_INFO})

    set_target_properties(suneigen4py
        PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
    )
endif()

# --------------------------------------------------------------------------------
#                            Build the c++ application
# --------------------------------------------------------------------------------

if(SUNEIGEN_BUILD_APP)

    set(APP_SOURCES
            app/main.cpp
    )

    add_executable(app ${APP_SOURCES})

    target_link_libraries(app
        PRIVATE
            ${LIBRARY_NAME}
            Eigen3::Eigen
    )

    target_set_warnings(app ENABLE ALL AS_ERROR ALL DISABLE Annoying)

    # Enable link-time-optimization if available for non-debug configurations
    target_enable_lto(app optimized)

    set_target_properties(app
        PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
        )

endif()
