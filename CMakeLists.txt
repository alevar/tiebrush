cmake_minimum_required(VERSION 2.8)
project(tiebrush)

# set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

set(CMAKE_CXX_STANDARD 11)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DNOCURL=1")

include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")

include(ExternalProject)
ExternalProject_Add(htslib
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND autoheader COMMAND autoconf COMMAND ./configure --without-libdeflate --disable-libcurl
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND ""
)

# very jankie solution to fix issues with linking libbigwig
ExternalProject_Add(libBigWig
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/libBigWig/
    BUILD_ALWAYS 1
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)
execute_process(
        COMMAND make -C ${CMAKE_SOURCE_DIR}/include/libBigWig
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/include/libBigWig
)

find_library(LBW_LIB BigWig "${CMAKE_SOURCE_DIR}/include/libBigWig")

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
target_link_libraries(tiebrush libhts.a -lz -llzma -lbz2 -pthread)

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

add_dependencies(tiecov htslib)
target_link_libraries(tiecov ${LBW_LIB} libhts.a -lz -llzma -lbz2 -pthread) # -lBigWig -lcrypto

install(TARGETS tiebrush tiecov DESTINATION bin)
