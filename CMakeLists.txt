cmake_minimum_required(VERSION 3.10)
project(dupcnt)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

add_executable(dupcnt
        dupcnt_core.cpp
        dupcnt_core.h
        kseq.h
        main.cpp )
target_link_libraries(dupcnt -lz)