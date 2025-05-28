cmake_minimum_required(VERSION 3.10)
project(RecursiveOptimizers)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(/usr/include/eigen3)

add_subdirectory(src)
add_subdirectory(tests)