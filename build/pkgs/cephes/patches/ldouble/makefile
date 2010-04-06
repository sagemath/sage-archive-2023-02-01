# Unix make file for long double precision function library
# Be sure to check the size of a pointer, LARGEMEM, in lcalc.h;
# the size of long double, XPD, in mconf.h; and the type of
# computer and endianness of the numbers in mconf.h.
CC= gcc
#CC = /a/gnu/linux/gcc/xgcc -B/a/gnu/linux/gcc/ -I/a/gnu/linux/gcc/include
#CFLAGS= -O0 -fno-builtin
CFLAGS= -O2 -Wall -Wuninitialized -fno-builtin
LIBS=
OBJS= acoshl.o asinhl.o asinl.o atanhl.o atanl.o bdtrl.o btdtrl.o cbrtl.o \
chdtrl.o coshl.o ellpel.o ellpkl.o elliel.o ellikl.o ellpjl.o \
exp10l.o exp2l.o expl.o fdtrl.o gammal.o gdtrl.o igamil.o igaml.o \
incbetl.o incbil.o isnanl.o j0l.o j1l.o jnl.o ldrand.o log10l.o log2l.o \
logl.o nbdtrl.o ndtril.o ndtrl.o pdtrl.o powl.o powil.o sinhl.o sinl.o \
sqrtl.o stdtrl.o tanhl.o tanl.o unityl.o ynl.o \
floorl.o polevll.o unityl.o mtherr.o
# cmplxl.o clogl.o

all: mtstl lcalc fltestl nantst testvect monotl libml.a

mtstl: libml.a mtstl.o $(OBJS)
	$(CC) $(CFLAGS) -o mtstl mtstl.o libml.a $(LIBS)
#	aout2exe mtstl

mtstl.o: mtstl.c

lparanoi: libml.a lparanoi.o setprec.o ieee.o econst.o $(OBJS)
	$(CC) $(CFLAGS) -o lparanoi lparanoi.o setprec.o ieee.o econst.o libml.a $(LIBS)
#	aout2exe lparanoi

lparanoi.o: lparanoi.c
	$(CC) $(CFLAGS) -Wno-implicit -c lparanoi.c

econst.o: econst.c ehead.h

lcalc: libml.a lcalc.o ieee.o econst.o $(OBJS)
	$(CC) $(CFLAGS) -o lcalc lcalc.o ieee.o econst.o libml.a $(LIBS)

lcalc.o: lcalc.c lcalc.h ehead.h

ieee.o: ieee.c ehead.h

# Use $(OBJS) in ar command for libml.a if possible; else *.o
libml.a: $(OBJS) mconf.h
	ar -rv libml.a $(OBJS)
	ranlib libml.a
#	ar -rv libml.a *.o

# Use assembly language sqrt and floor if available.
# sqrtl.c does not give strictly rounded results.
# The floor.c routines are slow.
#floorl.o: floorl.387
#	as -o floorl.o floorl.387
#sqrtl.o: sqrtl.387
#	as -o sqrtl.o sqrtl.387
#setprec.o: setprec.387
#	as -o setprec.o setprec.387

# ELF versions for linux
#floorl.o: floorlelf.387
#	as -o floorl.o floorlelf.387
#sqrtl.o: sqrtlelf.387
#	as -o sqrtl.o sqrtlelf.387
#setprec.o: setprelf.387
#	as -o setprec.o setprelf.387

# 68K routines
# For Sun 3
#floorl.o: floorl.sun
#	as -o floorl.o floorl.sun
#setprec.o: setprec.688
#	as -o setprec.o setprec.688
#sqrtl.o: sqrtl.sun
#	as -o sqrtl.o sqrtl.sun

# For Motorola 68k sysv
# Thanks to Philippe De Muyter <phdm@info.ucl.ac.be>
#floorl.o: floorl.mot
#	as -o floorl.o floorl.mot
#setprec.o: setprec.mot
#	as -o setprec.o setprec.mot
#sqrtl.o: sqrtl.mot
#	as -o sqrtl.o sqrtl.mot

fltestl: fltestl.c libml.a
	$(CC) $(CFLAGS) -o fltestl fltestl.c libml.a

fltestl.o: fltestl.c

flrtstl: flrtstl.c libml.a
	$(CC) $(CFLAGS) -o flrtstl flrtstl.c libml.a

flrtstl.o: flrtstl.c

nantst: nantst.c libml.a
	$(CC) $(CFLAGS) -o nantst nantst.c libml.a

nantst.o: nantst.c

testvect: testvect.o libml.a
	$(CC) $(CFLAGS) -o testvect testvect.o libml.a

testvect.o: testvect.c
	$(CC) -g -c -o testvect.o testvect.c

monotl: monotl.o libml.a
	$(CC) $(CFLAGS) -o monotl monotl.o libml.a

monotl.o: monotl.c
	$(CC) -g -c -o monotl.o monotl.c

acoshl.o: acoshl.c
asinhl.o: asinhl.c
asinl.o: asinl.c
atanhl.o: atanhl.c
atanl.o: atanl.c
bdtrl.o: bdtrl.c
btdtrl.o: btdtrl.c
cbrtl.o: cbrtl.c
chdtrl.o: chdtrl.c
coshl.o: coshl.c
ellpel.o: ellpel.c
ellpkl.o: ellpkl.c
elliel.o: elliel.c
ellikl.o: ellikl.c
ellpjl.o: ellpjl.c
exp10l.o: exp10l.c
exp2l.o: exp2l.c
expl.o: expl.c
fdtrl.o: fdtrl.c
gammal.o: gammal.c
gdtrl.o: gdtrl.c
igamil.o: igamil.c
igaml.o: igaml.c
incbetl.o: incbetl.c
incbil.o: incbil.c
isnanl.o: isnanl.c
j0l.o: j0l.c
j1l.o: j1l.c
jnl.o: jnl.c
ldrand.o: ldrand.c
log10l.o: log10l.c
log2l.o: log2l.c
logl.o: logl.c
nbdtrl.o: nbdtrl.c
ndtril.o: ndtril.c
ndtrl.o: ndtrl.c
pdtrl.o: pdtrl.c
powl.o: powl.c
powil.o: powil.c
sinhl.o: sinhl.c
sinl.o: sinl.c
sqrtl.o: sqrtl.c
stdtrl.o: stdtrl.c
tanhl.o: tanhl.c
tanl.o: tanl.c
unityl.o: unityl.c
ynl.o: ynl.c
floorl.o: floorl.c
polevll.o: polevll.c
unityl.o: unityl.c
mtherr.o: mtherr.c

clean:
	rm -f *.o
	rm -f mtstl
	rm -f lcalc
	rm -f fltestl
	rm -f flrtstl
	rm -f nantst
	rm -f lparanoi
	rm -f testvect
	rm -f libml.a

# Run test programs
check: mtstl fltestl testvect monotl libml.a
	-mtstl
	-fltestl
	-testvect
	-monotl
