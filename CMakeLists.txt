cmake_minimum_required(VERSION 2.8)
project(tiebrush)

# set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

set(CMAKE_CXX_STANDARD 11)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -fpermissive -DNOCURL=1")

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

set(INC_DIR "${CMAKE_SOURCE_DIR}/include")
set(HTS_DIR "${INC_DIR}/htslib")

include_directories("${INC_DIR}")
link_directories("${HTS_DIR}")

include(ExternalProject)

ExternalProject_Add(libdeflate
        SOURCE_DIR ${INC_DIR}/libdeflate
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND $(MAKE) libdeflate.a COMMAND cp libdeflate.h libdeflate.a ${HTS_DIR}
        INSTALL_COMMAND ""
        )

ExternalProject_Add(htslib
        SOURCE_DIR ${INC_DIR}/htslib
        BUILD_IN_SOURCE 1
        DEPENDS libdeflate
        CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/config_htslib.sh
        BUILD_COMMAND $(MAKE) lib-static
        INSTALL_COMMAND ""
        )

add_executable(tiebrush src/tiebrush.cpp
        src/commons.h
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${INC_DIR}/gclib/GStr.h
        ${INC_DIR}/gclib/GStr.cpp
        ${INC_DIR}/gclib/GArgs.h
        ${INC_DIR}/gclib/GArgs.cpp
        ${INC_DIR}/gclib/GBase.h
        ${INC_DIR}/gclib/GBase.cpp)

add_dependencies(tiebrush libdeflate htslib)

if (APPLE AND TIEBRUSH_STATIC_BUILD)
    target_link_libraries(tiebrush ${HTS_DIR}/libhts.a ${HTS_DIR}/libdeflate.a z lzma bz2 pthread)
else ()
    target_link_libraries(tiebrush ${HTS_DIR}/libhts.a ${HTS_DIR}/libdeflate.a z lzma bz2 pthread)
endif ()

ExternalProject_Add(libBigWig
        GIT_REPOSITORY "https://github.com/alevar/libBigWig"
        UPDATE_COMMAND ""
        DOWNLOAD_DIR "${INC_DIR}/include/"
        SOURCE_DIR "${INC_DIR}/libBigWig"
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND "${LIBBIGWIG_MAKE_CMD}"
        INSTALL_COMMAND ""
        )

include_directories(${INC_DIR}/libBigWig)
link_directories(${INC_DIR}/libBigWig)

add_executable(tiecov src/tiecov.cpp
        src/commons.h
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${INC_DIR}/gclib/GStr.h
        ${INC_DIR}/gclib/GStr.cpp
        ${INC_DIR}/gclib/GArgs.h
        ${INC_DIR}/gclib/GArgs.cpp
        ${INC_DIR}/gclib/GBase.h
        ${INC_DIR}/gclib/GBase.cpp)

add_dependencies(tiecov libdeflate htslib libBigWig)

if (APPLE AND TIEBRUSH_STATIC_BUILD)
    target_link_libraries(tiecov ${INC_DIR}/libBigWig/libBigWig.a ${HTS_DIR}/libhts.a ${HTS_DIR}/libdeflate.a z lzma bz2 pthread)
else ()
    target_link_libraries(tiecov ${INC_DIR}/libBigWig/libBigWig.a ${HTS_DIR}/libhts.a ${HTS_DIR}/libdeflate.a z lzma bz2 pthread)
endif ()

install(TARGETS tiebrush DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS tiecov DESTINATION ${CMAKE_INSTALL_PREFIX})
install(PROGRAMS tiewrap.py DESTINATION ${CMAKE_INSTALL_PREFIX})
