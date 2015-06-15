#new makefile for lrslib-050    2012.9.27

#contains multithread version of lrs called plrs, wrapper written by Gary Roumanis
#if Boost libraries are not available, comment out plrs compiles     http://www.boost.org/

#make all               uses gmp and long arithmetic 
#make allmp             uses native mp and long arithmetic

#Select a path below to give location of boost atomic library
# versions of gcc since at least 4.6.4 include boost atomic library
# sony, mune
BOOSTINC = /usr/include
BOOSTLIB = /usr/lib
#mai64 mai12
#BOOSTINC = /usr/include/boost_1_57_0
#BOOSTLIB = /usr/include/boost_1_57_0/stage/lib

# cgm-server.mcgill.ca
#BOOSTINC = /home/cgm/avis/C/ve/boost/cgm/boost_1_57_0/
#BOOSTLIB = /home/cgm/avis/C/ve/boost/cgm/boost_1_57_0/stage/lib
# obsolete version for cgm
#BOOSTINC = /usr/include/boost151/boost/include/
#BOOSTLIB = /usr/include/boost151/boost/lib/

#Select a path below to give location of gmp library

#cygwin
INCLUDEDIR = /usr/include
LIBDIR     = /usr/lib

#linux at mcgill with gmp version 3
#INCLUDEDIR = /usr/local/include
#LIBDIR     = /usr/local/lib

CFLAGS=-O3
# These flags should *not* include the arithmetic selecting define.
CPPFLAGS=-DLRS_QUIET -DTIMES -DSIGNALS
# set to something more useful if your system has a ranlib command
RANLIB ?= /bin/true

# default set of executables to build
BINARIES=2nash lrs lrs1 nash redund redund1 setnash setnash2
# default set of libraries to build
LIB=liblrsgmp.a
LIBRARIES=$(LIB) $(SHLIB) $(SHLINK)

# where to install binaries, libraries, include files
prefix := /usr/local

# Shared library
SONAME ?=liblrsgmp.so.0
SOMINOR ?=.0.0
SHLIB ?=$(SONAME)$(SOMINOR)
SHLINK ?=liblrsgmp.so

LRSGMPLIB=$(LIB)

# rule to build gmp arithmetic using object files.
%-GMP.o: %.c lrsgmp.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -DGMP -o $@ $<

# rule to build gmp arithmetic and relocatable object files.
%-GMP-SHR.o: %.c lrsgmp.h
	$(CC) -c -fPIC $(CFLAGS) $(CPPFLAGS) -DGMP -o $@ $<

# Cancel built in rules
%: %.c
%: %.cpp

# How to build a gmp arithmetic using tool
%: %-GMP.o $(LRSGMPLIB)
	$(CC) $< -L. -llrsgmp -L${LIBDIR}  -lgmp -o $@

all: $(BINARIES)

all-static: $(LIB)
	$(MAKE) LRSGMPLIB=$(LIB) all

all-shared: $(SHLINK)
	$(MAKE) LRSGMPLIB=$(SHLINK)

$(BINARIES): $(LRSGMPLIB)

liblrsgmp.a: lrslib-GMP.o lrsgmp-GMP.o
	ar r $@ $^
	$(RANLIB) $@

$(SHLIB): lrslib-GMP-SHR.o lrsgmp-GMP-SHR.o
	$(CC) -shared -Wl,-soname=$(SONAME) $(SHLIBFLAGS) -o $@ $^ -lgmp

$(SHLINK): $(SHLIB)
	ln -sf $< $@

lrs1:	lrs.c lrslib.c lrslib.c lrslong.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DLONG -o lrs1 lrs.c lrslib.c lrslong.c

redund1: redund.c lrslib.c lrslib.c lrslong.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DLONG -o redund1 redund.c lrslib.c lrslong.c

setnash: setupnash.c lrslib.c lrsmp.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o setnash setupnash.c lrslib.c lrsmp.c

setnash2: setupnash2.c lrslib.c lrsmp.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o setnash2 setupnash2.c lrslib.c lrsmp.c

install: install-static

install-shared: all-shared install-common
	 mkdir -p $(DESTDIR)${prefix}/lib
	 install -t $(DESTDIR)${prefix}/lib $(SHLIB)
	 cd $(DESTDIR)${prefix}/lib && ln -sf $(SHLIB) $(SHLINK)
	 cd $(DESTDIR)${prefix}/lib && ln -sf $(SHLIB) $(SONAME)

install-static: all-static install-common
	 mkdir -p $(DESTDIR)${prefix}/lib
	 install -t $(DESTDIR)${prefix}/lib $(LIB)

install-common:
	 mkdir -p $(DESTDIR)${prefix}/bin
	 install -t $(DESTDIR)${prefix}/bin $(BINARIES)
	 mkdir -p $(DESTDIR)${prefix}/include
	 install -t $(DESTDIR)${prefix}/include lrslib.h lrsgmp.h

fourier:	fourier.c lrslib.h lrslib.c lrsgmp.h lrsgmp.c
	gcc -O3 -DTIMES -DSIGNALS  -DGMP -I${INCLUDEDIR} fourier.c lrslib.c lrsgmp.c -L${LIBDIR}  -lgmp -o fourier


plrs:	  	plrs.cpp plrs.hpp lrslib.c lrslib.h lrsmp.c lrsmp.h
		g++ -DGMP  -Wall -Wno-write-strings -Wno-sign-compare -I${BOOSTINC}  -Wl,-rpath=${BOOSTLIB} -O3 -DPLRS -DGMP -o plrs plrs.cpp lrslib.c lrsgmp.c -L${BOOSTLIB} -lboost_thread -lboost_system -lgmp

plrsmp:	  	plrs.cpp plrs.hpp lrslib.c lrslib.h lrsmp.c lrsmp.h
		g++ -Wall -Wno-write-strings -Wno-sign-compare -Wno-unused-variable -I${BOOSTINC}  -L${BOOSTLIB} -Wl,-rpath=${BOOSTLIB} -O3 -DPLRS -DLONG -o plrs1 plrs.cpp lrslib.c lrslong.c -lboost_thread -lboost_system
		g++ -Wall -Wno-write-strings -Wno-sign-compare -Wno-unused-variable -I${BOOSTINC}  -L${BOOSTLIB} -Wl,-rpath=${BOOSTLIB} -O3 -DPLRS -o plrsmp plrs.cpp lrslib.c lrsmp.c  -lboost_thread -lboost_system

allmp:		lrs.c lrslib.c lrslib.h lrsmp.c lrsmp.h
		gcc -Wall -O3 -DTIMES -DSIGNALS -o lrs lrs.c lrslib.c lrsmp.c
		gcc -Wall -O3 -DTIMES -DSIGNALS -DLONG -o lrs1 lrs.c lrslib.c lrslong.c
		gcc -O3 -DTIMES -DSIGNALS -o redund  redund.c lrslib.c lrsmp.c
		gcc -O3 -DTIMES -DSIGNALS -DLONG -o redund1  redund.c lrslib.c lrslong.c
		gcc -O3 -DLRS_QUIET  -DTIMES -DSIGNALS -o nash nash.c lrslib.c lrsmp.c
		gcc -O3 -o setnash setupnash.c lrslib.c lrsmp.c
		gcc -O3 -o setnash2 setupnash2.c lrslib.c lrsmp.c
		gcc -O3 -o 2nash 2nash.c


clean:		
		rm -f $(BINARIES) $(LIBRARIES) plrs *.o *.exe
