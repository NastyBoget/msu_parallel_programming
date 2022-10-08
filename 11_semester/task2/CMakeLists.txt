cmake_minimum_required(VERSION 3.6)
project(task2)
include_directories(/usr/include/local)
set(CMAKE_CXX_COMPILER "/usr/local/bin/mpic++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(task2.out task2.cpp)
