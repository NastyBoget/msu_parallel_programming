cmake_minimum_required(VERSION 3.6)
project(task3)
include_directories(/usr/include/local)
set(CMAKE_CXX_COMPILER "/usr/local/bin/mpic++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")

add_executable(task3.out main.cpp include/grid.h include/solver.h include/functions.h include/saver.h include/index.h)
