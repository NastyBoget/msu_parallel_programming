cmake_minimum_required(VERSION 3.6)
project(task3)
include_directories(/usr/include/local)
set(CMAKE_CXX_COMPILER "/usr/local/bin/mpic++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(task3.out main_mpi.cpp include/solver_mpi.h include/grid.h include/index.h include/functions.h include/splitter.h)
