
# List all files containing tests.
set(TESTFILES
        main.cpp
        test_sunlinsol_superlu.cpp
        cvsRoberts_dns.cpp
        model_robertson.cpp
)

# Make unit_tests target
add_executable(unit_tests ${TESTFILES})

target_link_libraries(unit_tests
    PRIVATE
        ${LIBRARY_NAME}
        doctest
        Eigen3::Eigen
)

set_target_properties(unit_tests
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)
target_set_warnings(unit_tests
        ENABLE ALL
        AS_ERROR ALL
        DISABLE Annoying
)

set_target_properties(unit_tests PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

set(TEST_RUNNER_PARAMS "")
add_test(
        # Use some per-module/project prefix so that it is easier to run only tests for this module
        NAME ${LIBRARY_NAME}.unit_tests
        COMMAND unit_tests ${TEST_RUNNER_PARAMS}
)

# Adds a 'coverage' target.
include(CodeCoverage)
