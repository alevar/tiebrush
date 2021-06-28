#!/bin/sh
autoheader
autoconf
export CPPFLAGS="$CPPFLAGS -I." 
export LDFLAGS="$LDFLAGS -L."
./configure --with-libdeflate --disable-libcurl
