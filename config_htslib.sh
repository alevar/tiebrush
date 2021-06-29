#!/bin/sh
#autoheader
#autoconf
export CPPFLAGS="$CPPFLAGS -I." 
export LDFLAGS="$LDFLAGS -L."
#./configure --with-libdeflate --disable-libcurl
/bin/rm -f config.h
sed -i -e 's/-lcurl/-ldeflate/' Makefile
sed -i -e 's/HAVE_LIBCURL 1/HAVE_LIBDEFLATE 1/' Makefile
