cmake_minimum_required(VERSION 2.8)
# set(CMAKE_VERBOSE_MAKEFILE 1)

# basic
project(cuKolmogorov2D)
include_directories(.)
add_definitions(-std=c++11 -Wall)

# libraries
find_package(CUDA REQUIRED)
add_subdirectory(cujak)

cuda_add_library(cuKolmogorov2D STATIC
  Converter.cu
)
target_link_libraries(cuKolmogorov2D cujak)

add_subdirectory(test)
