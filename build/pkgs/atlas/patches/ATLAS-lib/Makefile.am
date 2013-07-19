ACLOCAL_AMFLAGS = -I m4

# Upstream doesn't version libraries, so this is kind of arbitrary
SO_VERSION=3:0

# You shouldn't have to customize anything from here on

ATLAS_SERIAL_LIBS=libatlas.la libcblas.la libf77blas.la liblapack.la

ATLAS_PARALLEL_LIBS=libptcblas.la libptf77blas.la libptlapack.la

all: all_parallel

all_serial: $(ATLAS_SERIAL_LIBS)

all_parallel: all_serial $(ATLAS_PARALLEL_LIBS)

clean:
	$(RM) -r -f .libs *.la *-obj

install: install_parallel

install_serial: $(ATLAS_SERIAL_LIBS)
	$(MKDIR_P) @libdir@
	@LIBTOOL@ --mode=install $(INSTALL) -c libatlas.la @libdir@/libatlas.la
	@LIBTOOL@ --mode=install $(INSTALL) -c libcblas.la @libdir@/libcblas.la
	@LIBTOOL@ --mode=install $(INSTALL) -c libf77blas.la @libdir@/libf77blas.la
	@LIBTOOL@ --mode=install $(INSTALL) -c liblapack.la @libdir@/liblapack.la
	@LIBTOOL@ --finish @libdir@

install_parallel: install_serial $(ATLAS_PARALLEL_LIBS)
	$(MKDIR_P) @libdir@
	@LIBTOOL@ --mode=install $(INSTALL) -c libptcblas.la @libdir@/libptcblas.la
	@LIBTOOL@ --mode=install $(INSTALL) -c libptf77blas.la @libdir@/libptf77blas.la
	@LIBTOOL@ --mode=install $(INSTALL) -c libptlapack.la @libdir@/libptlapack.la
	@LIBTOOL@ --finish @libdir@


# In presence of multiple definitions of thread related functions, use the ones from *_mut.o
# Works with all nm output formats (BSD/POSIX/System V)
libatlas.la: libatlas.a
	-$(RM) -r -f libatlas-obj
	$(MKDIR) libatlas-obj
	cd libatlas-obj && $(AR) x ../$<
	if [ `$(NM) -g $< | $(EGREP) -w 'ATL_(Set|Reset|Free|Dec)AtomicCount' | $(GREP) -w T | wc -l` -gt 4 ]; then \
	  $(RM) `ls -1 libatlas-obj/ATL_{Set,Reset,Free,Dec}AtomicCount_*.o | $(GREP) -v '_mut.o$$'`; \
	fi
	@LIBTOOL@ --tag=CC --mode=link $(CC) @LIBTOOL_TYPE@ \
	    -o $@ libatlas-obj/*.o \
	    @PTHREAD_LIB@ -lm \
	    -rpath @libdir@ -version-info $(SO_VERSION)


libcblas.la: libcblas.a libatlas.la
	-$(RM) -r -f libcblas-obj
	$(MKDIR) libcblas-obj
	cd libcblas-obj && $(AR) x ../$<
	@LIBTOOL@ --tag=CC --mode=link $(CC) @LIBTOOL_TYPE@ \
	    -o $@ libcblas-obj/*.o \
	    -latlas \
	    -rpath @libdir@ -version-info $(SO_VERSION)


libptcblas.la: libptcblas.a libatlas.la
	-$(RM) -r -f libptcblas-obj
	$(MKDIR) libptcblas-obj
	cd libptcblas-obj && $(AR) x ../$<
	@LIBTOOL@ --tag=CC --mode=link $(CC) @LIBTOOL_TYPE@ \
	    -o $@ libptcblas-obj/*.o \
	    @PTHREAD_LIB@ -latlas \
	    -rpath @libdir@ -version-info $(SO_VERSION)


libf77blas.la: libf77blas.a libatlas.la
	-$(RM) -r -f libf77blas-obj
	$(MKDIR) libf77blas-obj
	cd libf77blas-obj && $(AR) x ../$<
	@LIBTOOL@ --tag=F77 --mode=link $(F77) @LIBTOOL_TYPE@ \
	    -o $@ libf77blas-obj/*.o \
	    -latlas \
	    -rpath @libdir@ -version-info $(SO_VERSION)

libptf77blas.la: libptf77blas.a libatlas.la
	-$(RM) -r -f libptf77blas-obj
	$(MKDIR) libptf77blas-obj
	cd libptf77blas-obj && $(AR) x ../$<
	@LIBTOOL@ --tag=F77 --mode=link $(F77) @LIBTOOL_TYPE@ \
	    -o $@ libptf77blas-obj/*.o \
	    @PTHREAD_LIB@ -latlas \
	    -rpath @libdir@ -version-info $(SO_VERSION)

libf77refblas.la: libf77refblas.a
	-$(RM) -r -f libf77refblas-obj
	$(MKDIR) libf77refblas-obj
	cd libf77refblas-obj && $(AR) x ../$<
	@LIBTOOL@ --tag=F77 --mode=link $(F77) @LIBTOOL_TYPE@ \
	    -o $@ libf77refblas-obj/*.o \
	    -latlas \
	    -rpath @libdir@ -version-info $(SO_VERSION)

liblapack.la: liblapack.a libatlas.la libcblas.la libf77blas.la
	-$(RM) -r -f liblapack-obj
	$(MKDIR) liblapack-obj
	cd liblapack-obj && $(AR) x ../$<
	@LIBTOOL@ --tag=F77 --mode=link $(F77) @LIBTOOL_TYPE@ \
	    -o $@ liblapack-obj/*.o \
	    -latlas -lcblas -lf77blas -lm \
	    -rpath @libdir@ -version-info $(SO_VERSION)

libptlapack.la: libptlapack.a libatlas.la libptcblas.la libptf77blas.la
	-$(RM) -r -f libptlapack-obj
	$(MKDIR) libptlapack-obj
	cd libptlapack-obj && $(AR) x ../$<
	@LIBTOOL@ --tag=F77 --mode=link $(F77) @LIBTOOL_TYPE@ \
	    -o $@ libptlapack-obj/*.o \
	    @PTHREAD_LIB@ -latlas -lptcblas -lptf77blas -lm \
	    -rpath @libdir@ -version-info $(SO_VERSION)


.PHONY: all all_serial all_parallel build install install_serial install_parallel clean
