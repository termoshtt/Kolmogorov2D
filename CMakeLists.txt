cmake_minimum_required(VERSION 2.8)
# set(CMAKE_VERBOSE_MAKEFILE 1)

# basic
project(cuKolmogorov2D)
include_directories(.)
add_definitions(-std=c++11 -Wall)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++11)

# libraries
find_package(CUDA REQUIRED)
add_subdirectory(cujak)

add_custom_command(
  OUTPUT Kolmogorov2D.pb.cc Kolmogorov2D.pb.h
  DEPENDS Kolmogorov2D.proto
  COMMAND "protoc" "Kolmogorov2D.proto" "--cpp_out=." "--python_out=."
)

cuda_add_library(cuKolmogorov2D STATIC
  Converter.cu
  Kolmogorov2D.pb.cc
)
target_link_libraries(cuKolmogorov2D cujak)

add_subdirectory(test)

add_custom_target(
    document
    COMMAND "doxygen" "doc/Doxyfile"
)
