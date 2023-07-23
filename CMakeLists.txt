cmake_minimum_required(VERSION 2.8)
project(tiebrush)

# set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

set(CMAKE_CXX_STANDARD 11)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive -DNOCURL=1")

option(TIEBRUSH_BUILD_LIBBIGWIG "Turn on/off building step of libBigWig (used for bioconda)." ON)
if (NOT TIEBRUSH_BUILD_LIBBIGWIG)
    add_definitions(-DTIEBRUSH_BUILD_LIBBIGWIG=0)
    set(LIBBIGWIG_MAKE_CMD "")
else ()
    add_definitions(-DTIEBRUSH_BUILD_LIBBIGWIG=1)
    set(LIBBIGWIG_MAKE_CMD "make")
endif ()

option(TIEBRUSH_NATIVE_BUILD "Architecture-specific optimizations, i.e., g++ -march=native." OFF)
if (TIEBRUSH_NATIVE_BUILD)
    add_definitions(-DTIEBRUSH_NATIVE_BUILD=1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -ipo -no-prec-div -fp-model fast=2")
    endif ()
endif ()

option(TIEBRUSH_STATIC_BUILD "Static build." OFF)
if (TIEBRUSH_STATIC_BUILD)
    add_definitions(-DTIEBRUSH_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        set (CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        message(WARNING "WARNING: Builds on Mac are never fully static.")
    else ()
        set (CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
    endif ()
    # on linux cmake adds -rdynamic automatically which clang can't handle in static builds
    # if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    #     set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    # endif ()
endif ()

include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")

include(ExternalProject)
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --without-libdeflate --disable-libcurl
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )

add_executable(tiebrush src/tiebrush.cpp
        src/commons.h
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

add_dependencies(tiebrush htslib)
if (APPLE AND TIEBRUSH_STATIC_BUILD)
    target_link_libraries(tiebrush ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
else ()
    target_link_libraries(tiebrush ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
endif ()

ExternalProject_Add(libBigWig
        GIT_REPOSITORY "https://github.com/alevar/libBigWig"
        UPDATE_COMMAND ""
        DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/include/"
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/include/libBigWig"
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND "${LIBBIGWIG_MAKE_CMD}"
        INSTALL_COMMAND ""
        )

include_directories(${CMAKE_SOURCE_DIR}/include/libBigWig)
link_directories(${CMAKE_SOURCE_DIR}/include/libBigWig)

add_executable(tiecov src/tiecov.cpp
        src/commons.h
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

add_dependencies(tiecov htslib libBigWig)

if (APPLE AND TIEBRUSH_STATIC_BUILD)
    target_link_libraries(tiecov ${CMAKE_SOURCE_DIR}/include/libBigWig/libBigWig.a ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
else ()
    target_link_libraries(tiecov ${CMAKE_SOURCE_DIR}/include/libBigWig/libBigWig.a ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
endif ()

install(TARGETS tiebrush DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS tiecov DESTINATION ${CMAKE_INSTALL_PREFIX})
install(PROGRAMS tiewrap.py DESTINATION ${CMAKE_INSTALL_PREFIX})
