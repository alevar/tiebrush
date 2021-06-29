#!/bin/sh
git submodule update --init --recursive
mkdir -p build 
cd build && /bin/rm -rf *
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j2
