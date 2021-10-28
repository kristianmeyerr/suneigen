cmake_minimum_required(VERSION 3.15)

# Set project name, version and c++ language
# Version numbers are available by including "config.h" in
# the source. See config.h.in for some more details.
project(suneigen)

# --------------------------------------------------------------------------------
#                         Options
# --------------------------------------------------------------------------------
#[[
# Set options via commandline (e.g. -DENABLE_LTO=[ON|OFF])
option(ENABLE_WARNINGS_SETTINGS "Allow target_set_warnings to add flags and defines.
                                 Set this to OFF if you want to provide your own warning parameters." ON)
option(ENABLE_LTO "Enable link time optimization" ON)
option(ENABLE_DOCTESTS "Include tests in the library. Setting this to OFF will remove all doctest related code.
                        Tests in tests/*.cpp will still be enabled." ON)
option(SUNEIGEN_BUILD_PYTHON "Enable the compilation of the python module." ON)
option(SUNEIGEN_BUILD_TESTS "Enable the compilation of the test files." ON)
option(SUNEIGEN_BUILD_APP "Enable the compilation of the app files." ON)
]]
# --------------------------------------------------------------------------------
#                         CMake modules
# --------------------------------------------------------------------------------
#[[
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
include(ConfigSafeGuards)
include(Colors)
include(CTest)
include(Doctest)
include(LTO)
include(Misc)
include(Warnings)

# Check for LTO support.
find_lto(CXX)

# --------------------------------------------------------------------------------
#                         Thirdparty dependencies
# --------------------------------------------------------------------------------

# Find eigen library
find_package(Eigen3 REQUIRED)
]]
# --------------------------------------------------------------------------------
#                            Build the library
# --------------------------------------------------------------------------------
#[[
set(LIBRARY_NAME suneigen)

set(LIB_SOURCES
        src/example.cpp
)

# Compile all sources into a library.
add_library(${LIBRARY_NAME} OBJECT ${LIB_SOURCES})

# Add an alias autodiff::autodiff to the target library autodiff
add_library(${LIBRARY_NAME}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

# Lib needs its header files, and users of the library must also see these (PUBLIC).
target_include_directories(${LIBRARY_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Install suneigen interface library
install(TARGETS ${LIBRARY_NAME} DESTINATION bin)

# There's also (probably) doctests within the library, so we need to see this as well.
target_link_libraries(${LIBRARY_NAME}
    PUBLIC
        doctest
        Eigen3::Eigen
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
#                            Build tests
# --------------------------------------------------------------------------------

if(SUNEIGEN_BUILD_TESTS)
    add_subdirectory(tests)
endif()
]]
# --------------------------------------------------------------------------------
#                            Build python module
# --------------------------------------------------------------------------------

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


# --------------------------------------------------------------------------------
#                            Build the c++ application
# --------------------------------------------------------------------------------
#[[
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

# ---------------------------------------------------
#   Summarize the build
# ---------------------------------------------------

message("")
message("------------------------------- Summary -------------------------------")
message("C++ compiler name: ${CMAKE_CXX_COMPILER_ID} at ${CMAKE_CXX_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
message("Source dir: ${CMAKE_SOURCE_DIR}")
message("Binary dir: ${CMAKE_BINARY_DIR}")
message("Install dir: ${CMAKE_INSTALL_PREFIX}")
message("C Flags: ${CMAKE_C_FLAGS}")
message("C++ Flags: ${CMAKE_CXX_FLAGS}")
message("----------------------------- Dependencies ----------------------------")

message("Found Eigen3: ${Eigen3_FOUND}")
if (Eigen3_FOUND)
    message("  Version: ${EIGEN_VERSION}")
    message("  Include dir: ${EIGEN_INCLUDE_DIRS}")
endif()
]]