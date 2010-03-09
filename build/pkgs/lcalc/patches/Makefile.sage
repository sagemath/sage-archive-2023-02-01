#
# Makefile for the L-function calculator version 0.2
#

#todo: set up a configuration file that will detect
#the operating system, whether gmp, or pari are installed and
#the location of all the relevant files,
#the c compiler, that will generate options which are specific to the
#compilers, optimization options depending on the chip, etc


# Comment out the following line to remove the use of pari's
# elliptic curve routines. Doing so disables the -e option.
# g++ with -DINCLUDE_PARI sends a #define INCLUDE_PARI to the preprocessor.

PARI_DEFINE = -DINCLUDE_PARI
#PREPROCESSOR_DEFINE = -DUSE_LONG_DOUBLE

#OPENMP_FLAG = -fopenmp

#PREPROCESSOR_DEFINE = -DUSE_MPFR
ifeq ($(PREPROCESSOR_DEFINE),-DUSE_MPFR)
  #GMP_FLAGS= -L/usr/local/lib -lmpfrcpp -lgmp -lmpfr
  GMP_FLAGS= gmpfrxx.o mpfr_mul_d.o -lmpfr -lgmp -lgmpxx -lm
else
  GMP_FLAGS=
endif


OS_NAME := $(shell uname)

CC = g++
#cc = /home/mrubinst/local/bin/gcc
#CC = /home/mrubinst/local/bin/g++
#LD = /home/mrubinst/local/bin/g++

#CC = /Users/michaelrubinstein/math/L/packages/gcc4.3/usr/local/bin/g++
#EXTRA= -pg
#EXTRA = -ftree-vectorize -ftree-vectorizer-verbose=5 -funroll-loops



#MACHINE_SPECIFIC_FLAGS = -Wno-long-double
MACHINE_SPECIFIC_FLAGS = -ffast-math -fPIC

#G4 = FALSE
#ifeq ($(G4),TRUE)
#MACHINE_SPECIFIC_FLAGS = -fast -mcpu=G4 -mtune=G4
#endif

#G5 = TRUE
ifeq ($(G5),TRUE)
   #MACHINE_SPECIFIC_FLAGS = -fast -mcpu=G5 -mtune=G5
   #MACHINE_SPECIFIC_FLAGS = -ffast-math -maltivec -mpowerpc -mpowerpc64 -ftree-vectorize -ftree-vectorizer-verbose=5 -funroll-loops

   #MACHINE_SPECIFIC_FLAGS = -ffast-math -mpowerpc -mpowerpc64 -ftree-vectorize -ftree-vectorizer-verbose=5 -funroll-loops
   MACHINE_SPECIFIC_FLAGS = -ffast-math -mpowerpc -mpowerpc64 -m64
   #MACHINE_SPECIFIC_FLAGS = -mpowerpc -mpowerpc64 -m64
endif

CCFLAGS =  $(CXXFLAGS) -O3  $(OPENMP_FLAG)  $(PREPROCESSOR_DEFINE) $(MACHINE_SPECIFIC_FLAGS) $(EXTRA)

#warning- O2 doesn't help with -DUSE_LONG_DOUBLE on mac, and actually seems to hurt, making runtime longer
#by a factor of 1.5


ifeq ($(PARI_DEFINE),-DINCLUDE_PARI)
    #location of pari.h.
    LOCATION_PARI_H = $(SAGE_LOCAL)/include/pari #usual location

    #location of libpari.a or of libpari.so
    #depending on whether static or dynamic libraries are being used.
    #On mac os x it's the former, on linux I think usually the latter.
    LOCATION_PARI_LIBRARY = $(SAGE_LOCAL)/lib #usual location
else
    #supplied as a dummy so as to avoid more ifeq's below
    LOCATION_PARI_H = .
    LOCATION_PARI_LIBRARY = .
endif



#INCLUDEFILES= -I../include -I../../packages/gcc4.3/usr/local/include
INCLUDEFILES= -I../include

#For Mac os x we omit shared library options

ifeq ($(OS_NAME),Darwin)
    DYN_OPTION=dynamiclib
else
    DYN_OPTION=shared
endif

LDFLAGS = -L$(LOCATION_PARI_LIBRARY) -lpari



#NOTES:
#for caedmon: the shared pari library is in a funny location: /usr/local/pari/pari-2.1.5/lib
#At compile time we need to specify that location with:
#    -L/usr/local/pari/pari-2.1.5/lib -lpari
#At runtime, the computer attempts to load the pari shared library, and if it isn't in a standard
#location, we can do two things.
#1) One page suggested:
#    -Xlinker -rpath -Xlinker /usr/local/pari/pari-2.1.5/lib
#2) The other option, not recommended, is to type at the unix prompt:
#    LD_LIBRARY_PATH=/usr/local/pari/pari-2.1.5/lib:$LD_LIBRARY_PATH
#    export LD_LIBRARY_PATH
#If this is not done correctly, at runntime one gets the error
#    ./lcalc: error while loading shared libraries: libpari.so.1: cannot open shared
#    object file: No such file or directory
#One can list, after compiling, dynamic dependencies with the command: ldd lcalc and it will
#become clear which libraries the computer can find.


INSTALL_DIR= /usr/local

#object files for the libLfunction library
OBJ_L = Lglobals.o Lgamma.o Lriemannsiegel.o Lriemannsiegel_blfi.o Ldokchitser.o

#object files for the command line program
OBJ2=$(OBJ_L) Lcommandline_globals.o Lcommandline_misc.o Lcommandline_numbertheory.o Lcommandline_values_zeros.o
OBJ3=$(OBJ2) Lcommandline_elliptic.o Lcommandline_twist.o Lcommandline.o cmdline.o
OBJECTS = $(OBJ3)

all:
#	make print_vars
	make libLfunction.so
	make lcalc
	make examples
#	make find_L
#	make test

print_vars:
	@echo OS_NAME = $(OS_NAME)

lcalc: $(OBJECTS)
	$(CC) $(CCFLAGS) $(INCLUDEFILES) $(OBJECTS) $(LDFLAGS) -o lcalc $(GMP_FLAGS)

examples:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) example_programs/example.cc libLfunction.so -o example_programs/example $(GMP_FLAGS)


proc:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) example_programs/proc.cc libLfunction.so -o example_programs/proc $(GMP_FLAGS)

test:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) example_programs/test.cc libLfunction.so -o example_programs/test $(GMP_FLAGS)

find_L:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) find_L_functions/find_L_functions.cc libLfunction.so -o find_L_functions/find_L $(GMP_FLAGS)

.cc.o:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -c $<
.c.o:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -c $<


Lglobals.o: ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h

#Lmisc.o: ../include/Lmisc.h ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h

Lgamma.o: ../include/Lgamma.h ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h

Lriemannsiegel.o: ../include/Lriemannsiegel.h ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h
Lriemannsiegel.o: ../include/Lint_complex.h ../include/Lgamma.h
Lriemannsiegel.o: ../include/Lmisc.h

Ldokchitser.o: ../include/Ldokchitser.h ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h

#all I need here is dependence on the definition of the L-function class
#and the Complex typedef
Lcommandline_globals.o: ../include/Lcommandline_globals.h ../include/L.h
Lcommandline_globals.o: ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h
#Lcommandline_globals.o: ../include/Lmisc.h ../include/Lgamma.h
#Lcommandline_globals.o: ../include/Lriemannsiegel.h
#Lcommandline_globals.o: ../include/Ldirichlet_series.h ../include/Lprint.h
#Lcommandline_globals.o: ../include/Lnumberzeros.h ../include/Lgram.h
#Lcommandline_globals.o: ../include/Lvalue.h ../include/Lfind_zeros.h

Lcommandline_misc.o: ../include/Lcommandline_misc.h ../include/L.h
Lcommandline_misc.o: ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h
Lcommandline_misc.o: ../include/Lmisc.h ../include/Lgamma.h
Lcommandline_misc.o: ../include/Lriemannsiegel.h
Lcommandline_misc.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lcommandline_misc.o: ../include/Lnumberzeros.h ../include/Lgram.h
Lcommandline_misc.o: ../include/Lvalue.h ../include/Lfind_zeros.h
Lcommandline_misc.o: ../include/Lcommandline_numbertheory.h
Lcommandline_misc.o: ../include/Lcommandline_globals.h

Lcommandline_numbertheory.o: ../include/Lcommandline_numbertheory.h
Lcommandline_numbertheory.o: ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h

Lcommandline_values_zeros.o: ../include/Lcommandline_values_zeros.h
Lcommandline_values_zeros.o: ../include/L.h ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h
Lcommandline_values_zeros.o: ../include/Lint_complex.h ../include/Lmisc.h
Lcommandline_values_zeros.o: ../include/Lgamma.h ../include/Lriemannsiegel.h
Lcommandline_values_zeros.o: ../include/Ldirichlet_series.h
Lcommandline_values_zeros.o: ../include/Lprint.h ../include/Lnumberzeros.h
Lcommandline_values_zeros.o: ../include/Lgram.h ../include/Lvalue.h
Lcommandline_values_zeros.o: ../include/Lfind_zeros.h
Lcommandline_values_zeros.o: ../include/Lcommandline_numbertheory.h
Lcommandline_values_zeros.o: ../include/Lcommandline_globals.h

Lcommandline_elliptic.o: ../include/Lcommandline_elliptic.h ../include/L.h
Lcommandline_elliptic.o: ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h
#Lcommandline_elliptic.o: ../include/Lmisc.h ../include/Lgamma.h
Lcommandline_elliptic.o: ../include/Lriemannsiegel.h
Lcommandline_elliptic.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lcommandline_elliptic.o: ../include/Lnumberzeros.h ../include/Lgram.h
Lcommandline_elliptic.o: ../include/Lvalue.h ../include/Lfind_zeros.h
Lcommandline_elliptic.o: ../include/Lcommandline_numbertheory.h
Lcommandline_elliptic.o: ../include/Lcommandline_globals.h
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(LOCATION_PARI_H) $(PARI_DEFINE) -c Lcommandline_elliptic.cc

Lcommandline_twist.o: ../include/Lcommandline_twist.h ../include/L.h
Lcommandline_twist.o: ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h
#Lcommandline_twist.o: ../include/Lmisc.h ../include/Lgamma.h
Lcommandline_twist.o: ../include/Lriemannsiegel.h
Lcommandline_twist.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lcommandline_twist.o: ../include/Lnumberzeros.h ../include/Lgram.h
Lcommandline_twist.o: ../include/Lvalue.h ../include/Lfind_zeros.h
Lcommandline_twist.o: ../include/Lcommandline_numbertheory.h
Lcommandline_twist.o: ../include/Lcommandline_globals.h
Lcommandline_twist.o: ../include/Lcommandline_elliptic.h
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(LOCATION_PARI_H) $(PARI_DEFINE) -c Lcommandline_twist.cc

cmdline.o: ../include/cmdline.h ../include/getopt.h
#$(CC) $(CCFLAGS) $(INCLUDEFILES) -DHAVE_LONG_LONG -c cmdline.c


Lcommandline.o: ../include/Lcommandline.h ../include/L.h
Lcommandline.o: ../include/Lglobals.h ../include/Lcommon.h ../include/Lcomplex.h ../include/Lnumeric.h ../include/Lint_complex.h
Lcommandline.o: ../include/Lmisc.h ../include/Lgamma.h
Lcommandline.o: ../include/Lriemannsiegel.h ../include/Ldirichlet_series.h
Lcommandline.o: ../include/Lprint.h ../include/Lnumberzeros.h
Lcommandline.o: ../include/Lgram.h ../include/Lvalue.h
Lcommandline.o: ../include/Lfind_zeros.h
Lcommandline.o: ../include/Lcommandline_numbertheory.h
Lcommandline.o: ../include/Lcommandline_globals.h
Lcommandline.o: ../include/Lcommandline_misc.h
Lcommandline.o: ../include/Lcommandline_elliptic.h
Lcommandline.o: ../include/Lcommandline_twist.h
Lcommandline.o: ../include/Lcommandline_values_zeros.h
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(LOCATION_PARI_H) $(PARI_DEFINE) -c Lcommandline.cc


libLfunction.so: $(OBJ_L)
	g++ -$(DYN_OPTION)  -o libLfunction.so $(OBJ_L)

clean:
	rm -f *.o lcalc libLfunction.so example_programs/example

install:
	cp -f lcalc $(INSTALL_DIR)/bin/.
	cp -f libLfunction.so $(INSTALL_DIR)/lib/.
	cp -rf ../include $(INSTALL_DIR)/include/Lfunction


SRCS = Lcommandline.cc Lcommandline_elliptic.cc Lcommandline_globals.cc Lcommandline_misc.cc Lcommandline_numbertheory.cc Lcommandline_twist.cc Lcommandline_values_zeros.cc Lgamma.cc Lglobals.cc Lmisc.cc Lriemannsiegel.cc Lriemannsiegel_blfi.cc cmdline.c
depend:
	makedepend -f depends -- $(CCFLAGS) -Y../include -- $(SRCS)
