set(APP_NAME suneigen-app)

set(APP_SOURCES
    main.cpp
    model_robertson.cpp
)

add_executable(${APP_NAME} ${APP_SOURCES})

target_link_libraries(${APP_NAME}
    PRIVATE
        ${LIBRARY_NAME}
        Eigen3::Eigen
)

target_include_directories(${APP_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/app
)

target_set_warnings(${APP_NAME} ENABLE ALL AS_ERROR ALL DISABLE Annoying)

# Enable link-time-optimization if available for non-debug configurations
target_enable_lto(${APP_NAME} optimized)

# Use the same visiliblity flag as in pybind11
target_compile_options(${APP_NAME}
    PRIVATE
        -fvisibility=hidden
)

set_target_properties(${APP_NAME}
    PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)
