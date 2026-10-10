cmake_minimum_required(VERSION 3.10)
project(EntropyStableEulerTests CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_executable(EntropyStableEulerTests EntropyStableEulerTests.cpp)
target_include_directories(EntropyStableEulerTests PRIVATE ../include)
add_test(NAME EntropyStableEulerTests COMMAND EntropyStableEulerTests)
