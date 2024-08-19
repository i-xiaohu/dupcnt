cmake_minimum_required(VERSION 3.10)
project(dupcnt)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

set(cstl
        cstl/kvec.h
        cstl/khash.h
        cstl/ksort.h
        cstl/kbtree.h
        cstl/kstring.c cstl/kstring.h
        cstl/kthread.c cstl/kthread.h)

set(bwalib
        bwalib/kseq.h
        bwalib/kopen.c bwalib/kopen.h
        bwalib/utils.c bwalib/utils.h)

set(FM_index2
        FM_index2/bseq1.h
        FM_index2/FMI_search.cpp FM_index2/FMI_search.h
        FM_index2/read_index_ele.cpp FM_index2/read_index_ele.h
        FM_index2/bntseq.cpp FM_index2/bntseq.h
        FM_index2/sais.h)

add_executable(mem2idx ${FM_index2} bwalib/utils.c bwalib/utils.h FM_index2/index_main.cpp)
target_link_libraries(mem2idx -lz -lrt)

add_executable(dupcnt
        ${FM_index2}
        bwalib/utils.c  bwalib/utils.h
        input.cpp input.h
        ${cstl}
        dupcnt_core.cpp
        dupcnt_core.h
        main.cpp )
target_link_libraries(dupcnt -lz -lrt -lpthread)