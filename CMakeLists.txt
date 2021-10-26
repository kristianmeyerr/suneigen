cmake_minimum_required(VERSION 3.15)

project(suneigen)

add_subdirectory(thirdparty/pybind11)
pybind11_add_module(suneigen4py python/main.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(suneigen4py PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
