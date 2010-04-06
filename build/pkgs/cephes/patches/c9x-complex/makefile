# Makefile for C9X complex variables library
# Edit endian-ness and computer type properly in mconf.h before compiling.
# For additional test programs see also celefunt.mak.

# Caution: the tests contain some references to the cephes math library
# or references to standard float or long double functions that might
# not be available.
# This should be enough to run at least the double precision functions:
LIBS = stubs.o -lm
# These are intended to get the cephes libraries:
#LIBS = -lml -lmd -lmf
#LIBS = ../libml.a ../libmd.a ../libmf.a

CC = gcc
#CC = /b/gnu/gcc-build/gcc/xgcc -B/b/gnu/gcc-build/gcc/ -I/b/gnu/gcc-build/gcc
# For 128-bit wide long double
#PRECFLAGS = -DLD128BITS
PRECFLAGS =
CFLAGS = $(PRECFLAGS) -I. -g -O2 -Wall
#CFLAGS = $(PRECFLAGS) -I. -g -O1 -Wall
#CFLAGS = $(PRECFLAGS) -I. -g -O0 -Wall

# float precision
SFILES = cmplxf.o clogf.o cgammaf.o
# double precision
DFILES = cmplx.o clog.o cgamma.o
# long double precision
#LFILES =
LFILES = cmplxl.o clogl.o cgammal.o
LIBMCFILES = $(SFILES) $(DFILES) $(LFILES)

# mtstf and mtstl test the float and long double functions.
all: libmc.a mtst cvect dccalc mtstf whitebxf whitebx

mtst: mtst.o $(DFILES) stubs.o
	$(CC) -v -o mtst mtst.o $(DFILES) $(LIBS)

mtstf: mtstf.o $(SFILES) stubs.o
	$(CC) -o mtstf mtstf.o $(SFILES) $(LIBS)

mtstl: mtstl.o $(LFILES) stubs.o
	$(CC) -o mtstl mtstl.o $(LFILES) $(LIBS)

cvect: cvect.o clog.o cmplx.o stubs.o
	$(CC) -o cvect cvect.o clog.o cmplx.o $(LIBS)

#cvect-nans: cvect-nans.o clog-nans.o cmplx-nans.o
#	$(CC) -o cvect-nans cvect-nans.o clog-nans.o cmplx-nans.o $(LIBS)

whitebxf: whitebxf.o libmc.a
	$(CC) $(CFLAGS) -o whitebxf whitebxf.o libmc.a $(LIBS)

whitebxf.o: whitebxf.c complex.h

whitebx: whitebx.o  libmc.a
	$(CC) $(CFLAGS) -o whitebx whitebx.o libmc.a $(LIBS)

whitebx.o: whitebx.c complex.h


libmc.a: $(LIBMCFILES)
	rm -f libmc.a
	ar rv libmc.a $(LIBMCFILES)
	ranlib libmc.a

# Complex variable calculator program.
# If you don't have the readline library, leave out the references to it
# (and curses) or get the library from ftp.gnu.org.
# If compiling on a 16-bit computer, you may need to set LARGEMEM 0
# in dccalc.h.
dccalc: dccalc.o
#	$(CC) -o dccalc dccalc.o libmc.a -lreadline -lncurses $(LIBS)
	$(CC) -o dccalc dccalc.o libmc.a $(LIBS)

dccalc.o: dccalc.c dccalc.h
#	$(CC) $(CFLAGS) -DUSE_READLINE=1 -c dccalc.c
	$(CC) $(CFLAGS) -c dccalc.c

clog.o: clog.c complex.h protos.h
cmplx.o: cmplx.c complex.h protos.h
clogf.o: clogf.c complex.h protos.h
cmplxf.o: cmplxf.c complex.h protos.h
clogl.o: clogl.c complex.h protos.h
cmplxl.o: cmplxl.c complex.h protos.h
#clog-nans.o: clog-nans.c complex.h protos.h
#cmplx-nans.o: cmplx-nans.c complex.h protos.h

check: mtst cvect mtstf whitebxf whitebx
	-./mtstf
	-./cvect
	-./mtst
	-./whitebxf
	-./whitebx
	echo For additional tests, make -f celefunt.mak

clean:
	rm -f *.o
	rm -f mtst
	rm -f cvect
	rm -f mtstf
	rm -f mtstl
	rm -f cvect-nans
	rm -f dccalc
	rm -f libmc.a
	rm -f tcabs
	rm -f tcabsf
	rm -f tcexp
	rm -f tcexpf
	rm -f tclog
	rm -f tclogf
	rm -f tcpow
	rm -f tcpowf
	rm -f tcsin
	rm -f tcsinf
	rm -f tcsqrt
	rm -f tcsqrtf
	rm -f whitebxf
	rm -f whitebx

#zip:
#	rm -f c9x-complex.zip
#	rm -f c9x-complex.shar
#	rm -f c9x-complex.shar.gz
#	zip c9x-complex.zip *.c *.h *.pat *.doc *.mak makefile read.me
#	shar *.c *.h *.pat *.doc *.mak makefile read.me > c9x-complex.shar
#	gzip c9x-complex.shar
