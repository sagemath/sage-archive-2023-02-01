AUTOMAKE_OPTIONS = gnu

lib_LTLIBRARIES = libcsage.la

libcsage_la_SOURCES = src/interrupt.c src/stdsage.c src/ntl_wrap.cc src/mpn_pylong.c src/mpz_pylong.c src/ZZ_pylong.cc

include_HEADERS = src/interrupt.h src/stdsage.h src/ccobject.h src/ntl_wrap.h src/mpn_pylong.h src/mpz_pylong.h src/ZZ_pylong.h

libcsage_la_LIBADD = -lntl -lgmp
libcsage_la_LDFLAGS = -version-info 0:0:0
