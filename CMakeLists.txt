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

set(FM_index
        FM_index/bntseq.c FM_index/bntseq.h
        FM_index/bwt.c FM_index/bwt.h
        FM_index/bwt_gen.c
        FM_index/is.c
        FM_index/QSufSort.c FM_index/QSufSort.h
        FM_index/rle.c FM_index/rle.h
        FM_index/rope.c FM_index/rope.h)

set(bwalib
        bwalib/kseq.h
        bwalib/bwashm.c
        bwalib/bwa.c bwalib/bwa.h
        bwalib/kopen.c bwalib/kopen.h
        bwalib/ksw.c bwalib/ksw.h
        bwalib/utils.c bwalib/utils.h)

add_executable(bwaidx ${FM_index} ${bwalib} FM_index/index_main.c)
target_link_libraries(bwaidx -lz -lrt)

add_executable(dupcnt
        ${FM_index}
        ${bwalib}
        ${cstl}
        dupcnt_core.cpp
        dupcnt_core.h
        main.cpp )
target_link_libraries(dupcnt -lz -lrt -lpthread)