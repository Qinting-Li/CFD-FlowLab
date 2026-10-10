cmake_minimum_required(VERSION 3.10)
project(entropy_stable_euler_flux LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

add_executable(entropy_stable_euler_tests
    tests/test_entropy_stable_euler.cpp
)

target_include_directories(entropy_stable_euler_tests PRIVATE include)

add_test(
    NAME entropy_stable_euler_tests
    COMMAND entropy_stable_euler_tests
)
