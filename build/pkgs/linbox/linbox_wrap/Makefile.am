AUTOMAKE_OPTIONS = gnu

lib_LTLIBRARIES = liblinboxwrap.la

liblinboxwrap_la_SOURCES = src/linbox_wrap.cpp

include_HEADERS = src/linbox_wrap.h

liblinboxwrap_la_LIBADD = -llinbox $(BLAS_LIBS)
liblinboxwrap_la_LDFLAGS = -version-info 0:0:0
