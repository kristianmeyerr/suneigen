
set(MODELFILES
        models/model_robertson.cpp)

# List all files containing tests.
set(TESTFILES
        main.cpp
        test_sunlinsol_superlu.cpp
        test_model_robertson.cpp
)

set(TEST_MAIN unit_tests)   # Default name for test executable
set(TEST_RUNNER_PARAMS "")  # Any arguemnts to feed the test runner

# Make unit_tests target
add_executable(${TEST_MAIN} ${TESTFILES} ${MODELFILES})

target_link_libraries(${TEST_MAIN}
    PRIVATE
        ${LIBRARY_NAME}
        doctest::doctest
        Eigen3::Eigen
)

set_target_properties(${TEST_MAIN}
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

target_set_warnings(${TEST_MAIN} ENABLE ALL AS_ERROR ALL DISABLE Annoying)

# Use the same visiliblity flag as in pybind11
target_compile_options(${TEST_MAIN}
    PRIVATE
        -fvisibility=hidden
)

set_target_properties(${TEST_MAIN} PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

add_test(
        # Use some per-module/project prefix so that it is easier to run only tests for this module
        NAME ${LIBRARY_NAME}.${TEST_MAIN}
        COMMAND ${TEST_MAIN} ${TEST_RUNNER_PARAMS}
)

# Adds a 'coverage' target.
include(../cmake/CodeCoverage.cmake)
