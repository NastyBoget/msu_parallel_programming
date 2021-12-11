cmake_minimum_required(VERSION 3.6)
project(task3)
# cmake file for clion
# orterun path: /usr/local/bin/orterun
include_directories(/usr/include/local)

set(CMAKE_CXX_COMPILER "/usr/local/bin/mpic++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(task3.out task3.cpp data_structures.h utils.cpp batcher_sort.cpp)
