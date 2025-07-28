cmake_minimum_required(VERSION 3.10)
project(tiebrush VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set debug build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Options
option(TIEBRUSH_BUILD_LIBBIGWIG "Turn on/off building step of libBigWig (used for bioconda)." ON)
option(TIEBRUSH_NATIVE_BUILD "Architecture-specific optimizations, i.e., g++ -march=native." OFF)
option(TIEBRUSH_STATIC_BUILD "Static build." OFF)

# Compiler definitions
if (NOT TIEBRUSH_BUILD_LIBBIGWIG)
    add_compile_definitions(TIEBRUSH_BUILD_LIBBIGWIG=0)
    set(LIBBIGWIG_MAKE_CMD "")
else ()
    add_compile_definitions(TIEBRUSH_BUILD_LIBBIGWIG=1)
    set(LIBBIGWIG_MAKE_CMD "make")
endif ()

if (TIEBRUSH_NATIVE_BUILD)
    add_compile_definitions(TIEBRUSH_NATIVE_BUILD=1)
    add_compile_options(-march=native)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        add_compile_options(-xHOST -ipo -no-prec-div -fp-model fast=2)
    endif ()
endif ()

if (TIEBRUSH_STATIC_BUILD)
    add_compile_definitions(TIEBRUSH_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        message(WARNING "WARNING: Builds on Mac are never fully static.")
    else ()
        set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
    endif ()
endif ()

# Common compiler flags
add_compile_options(-fpermissive -DNOCURL=1 -g -O0)

# Find required packages
find_package(Threads REQUIRED)

# Try to find pkg-config, but make it optional
find_package(PkgConfig QUIET)

# Find compression libraries with fallbacks
if(PkgConfig_FOUND)
    pkg_check_modules(ZLIB QUIET zlib)
    pkg_check_modules(LZMA QUIET liblzma)
    pkg_check_modules(BZ2 QUIET libbz2)
endif()

# Fallback: find libraries directly if pkg-config didn't work
if(NOT ZLIB_FOUND)
    find_library(ZLIB_LIBRARIES NAMES z)
    if(ZLIB_LIBRARIES)
        set(ZLIB_FOUND TRUE)
    endif()
endif()

if(NOT LZMA_FOUND)
    find_library(LZMA_LIBRARIES NAMES lzma)
    if(LZMA_LIBRARIES)
        set(LZMA_FOUND TRUE)
    endif()
endif()

if(NOT BZ2_FOUND)
    find_library(BZ2_LIBRARIES NAMES bz2)
    if(BZ2_LIBRARIES)
        set(BZ2_FOUND TRUE)
    endif()
endif()

# Set default libraries if not found (these are usually available on macOS)
if(NOT ZLIB_FOUND)
    set(ZLIB_LIBRARIES z)
    set(ZLIB_FOUND TRUE)
endif()

if(NOT LZMA_FOUND)
    set(LZMA_LIBRARIES lzma)
    set(LZMA_FOUND TRUE)
endif()

if(NOT BZ2_FOUND)
    set(BZ2_LIBRARIES bz2)
    set(BZ2_FOUND TRUE)
endif()

# External projects
include(ExternalProject)

# htslib external project
ExternalProject_Add(htslib
    GIT_REPOSITORY "https://github.com/samtools/htslib.git"
    GIT_TAG "1.18"
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --without-libdeflate --disable-libcurl
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a
)

# libBigWig external project
ExternalProject_Add(libBigWig
    GIT_REPOSITORY "https://github.com/alevar/libBigWig"
    UPDATE_COMMAND ""
    DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/include/"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/include/libBigWig"
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "${LIBBIGWIG_MAKE_CMD}"
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${CMAKE_SOURCE_DIR}/include/libBigWig/libBigWig.a
)

# Create imported targets for external libraries
add_library(hts STATIC IMPORTED)
set_target_properties(hts PROPERTIES
    IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/include/htslib
)

add_library(bigwig STATIC IMPORTED)
set_target_properties(bigwig PROPERTIES
    IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/include/libBigWig/libBigWig.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/include/libBigWig
)

# tiebrush executable
add_executable(tiebrush
    src/tiebrush.cpp
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
    ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp
)

target_include_directories(tiebrush PRIVATE
    ${CMAKE_SOURCE_DIR}/include/
    ${CMAKE_SOURCE_DIR}/include/htslib
)

target_link_libraries(tiebrush
    hts
    Threads::Threads
    ${ZLIB_LIBRARIES}
    ${LZMA_LIBRARIES}
    ${BZ2_LIBRARIES}
)

add_dependencies(tiebrush htslib)

# tiecov executable
add_executable(tiecov
    src/tiecov.cpp
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
    ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp
)

target_include_directories(tiecov PRIVATE
    ${CMAKE_SOURCE_DIR}/include/
    ${CMAKE_SOURCE_DIR}/include/htslib
    ${CMAKE_SOURCE_DIR}/include/libBigWig
)

target_link_libraries(tiecov
    bigwig
    hts
    Threads::Threads
    ${ZLIB_LIBRARIES}
    ${LZMA_LIBRARIES}
    ${BZ2_LIBRARIES}
)

add_dependencies(tiecov htslib libBigWig)

# Installation
install(TARGETS tiebrush tiecov
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}
)
install(PROGRAMS tiewrap.py DESTINATION ${CMAKE_INSTALL_PREFIX})
