cmake_minimum_required(VERSION 3.15)

# List all files containing tests.
set(TESTFILES
        main.cpp
        dummy.cpp
)

set(TEST_MAIN unit_tests)
set(TEST_RUNNER_PARAMS "")

# --------------------------------------------------------------------------------
#                         Make Tests.
# --------------------------------------------------------------------------------
add_executable(${TEST_MAIN} ${TESTFILES})
target_link_libraries(${TEST_MAIN} PRIVATE ${LIBRARY_NAME} doctest)
set_target_properties(${TEST_MAIN} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_set_warnings(${TEST_MAIN} ENABLE ALL AS_ERROR ALL DISABLE Annoying) # Set warnings (if needed).

set_target_properties(${TEST_MAIN} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

add_test(
        # Use some per-module/project prefix so that it is easier to run only tests for this module
        NAME ${LIBRARY_NAME}.${TEST_MAIN}
        COMMAND ${TEST_MAIN} ${TEST_RUNNER_PARAMS}
)



# Adds a 'coverage' target.
include(CodeCoverage)