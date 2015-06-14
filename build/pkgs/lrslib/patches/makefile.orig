# Makefile for lrslib including lrs and redund 
# 2000.6.14
#
# choose one of the first 4 and gmp if applicable - see README
#
# make all      normal version for 32bit machines with timing/signals handling
# make all64    for machines with 64bit integers, eg. DEC Alpha
# make ansi     ansi standard version for 32bit machines without signals handling
# make nosigs   ansi standard version for 32bit machines without timing/signals handling

# make demo     compile the demo programs (vanilla machine)

# make gmp      uses gmp arithmetic, set paths for include files *first*

# make clean    to clean all executables


#Select one of the following INCLUDE,LIB paths only needed for gmp version

#linux at mcgill with gmp version 3
#INCLUDEDIR = /usr/local/include
#LIBDIR     = /usr/local/lib

#linux at mcgill with gmp version 2
#INCLUDEDIR = /labs/cgm/gmp2/include
#LIBDIR     = /labs/cgm/gmp2/lib


#TRUE64 at mcgill gmp version 3
INCLUDEDIR = /labs/cgm/include
LIBDIR     = /labs/cgm/lib

#TRUE64 at mcgill gmp version 2
#INCLUDEDIR = /labs/cgm/gmp2/include
#LIBDIR     = /labs/cgm/gmp2/lib


all:	lrs.c lrslib.c lrslib.h lrsmp.c lrsmp.h lrslong.c lrslong.h redund.c buffer.c nash.c
	gcc -O3 -DTIMES -DSIGNALS -o lrs  lrs.c lrslib.c lrsmp.c
	gcc -O3 -DTIMES -DSIGNALS -o redund  redund.c lrslib.c lrsmp.c
	gcc -O3 -DTIMES -DSIGNALS -DLONG -o lrs1  lrs.c lrslib.c lrslong.c
	gcc -O3 -DTIMES -DSIGNALS -DLONG -o redund1  redund.c lrslib.c lrslong.c
	gcc -O3 -DLRS_QUIET  -DTIMES -DSIGNALS -o nash nash.c lrslib.c lrsmp.c
	gcc -O3 -o setupnash setupnash.c lrslib.c lrsmp.c
	gcc -O3 -o setupnash2 setupnash2.c lrslib.c lrsmp.c
	gcc -Wall -O3 -o fourier  fourier.c lrslib.c lrsmp.c
	gcc -O3 -o buffer buffer.c

gmp:	lrs.c redund.c lrslib.h lrslib.c lrsgmp.h lrsgmp.c nash.c
	gcc -O3 -static -DTIMES -DSIGNALS  -DGMP -I${INCLUDEDIR} lrs.c lrslib.c lrsgmp.c -L${LIBDIR}  -lgmp -o glrs
	gcc -O3 -static -DTIMES -DSIGNALS -DGMP -I${INCLUDEDIR} redund.c lrslib.c lrsgmp.c -L${LIBDIR} -lgmp -o gredund
	gcc -O3 -static -DLRS_QUIET -DTIMES -DSIGNALS -DGMP -I${INCLUDEDIR} nash.c lrslib.c lrsgmp.c -L${LIBDIR} -lgmp -o gnash
	gcc -O3 -static -DTIMES -DSIGNALS  -DGMP -I${INCLUDEDIR} fourier.c lrslib.c lrsgmp.c -L${LIBDIR}  -lgmp -o gfourier
	gcc -O3 -o buffer buffer.c

gnash:	lrslib.h lrslib.c lrsgmp.h lrsgmp.c nash.c
	gcc -O3 -static -DLRS_QUIET -DTIMES -DSIGNALS -DGMP -I${INCLUDEDIR} nash.c lrslib.c lrsgmp.c -L${LIBDIR} -lgmp -o gnash

all64:	lrs.c lrslib.c lrslib.h lrsmp.c lrsmp.h lrslong.c lrslong.h redund.c buffer.c
	gcc -DTIMES -DSIGNALS -DB64 -O3 -o lrs  lrs.c lrslib.c lrsmp.c
	gcc -DTIMES -DSIGNALS -DB64 -O3 -o redund  redund.c lrslib.c lrsmp.c
	gcc -DTIMES -DSIGNALS -DLONG -DB64 -O3 -o lrs1  lrs.c lrslib.c lrslong.c
	gcc -DTIMES -DSIGNALS -DLONG -DB64 -O3 -o redund1  redund.c lrslib.c lrslong.c
	gcc -O3 -o buffer buffer.c

ansi:	lrs.c lrslib.c lrslib.h lrsmp.c lrsmp.h lrslong.c lrslong.h redund.c buffer.c nash.c
	gcc -ansi -DTIMES   -O3 -o lrs  lrs.c lrslib.c lrsmp.c
	gcc -ansi -DTIMES   -O3 -o redund redund.c lrslib.c lrsmp.c
	gcc -ansi -DTIMES  -DLONG  -O3 -o lrs1  lrs.c lrslib.c lrslong.c
	gcc -ansi -DTIMES  -DLONG  -O3 -o redund1 redund.c lrslib.c lrslong.c
	gcc -O3 -o buffer buffer.c
	gcc -Wall -ansi -O3 -o nash nash.c lrslib.c lrsmp.c

nosigs:	lrs.c lrslib.c lrslib.h lrsmp.c lrsmp.h lrslong.c lrslong.h redund.c buffer.c
	gcc -ansi  -O3 -o lrs  lrs.c lrslib.c lrsmp.c
	gcc -ansi  -O3 -o redund redund.c lrslib.c lrsmp.c
	gcc -ansi  -O3 -DLONG -o lrs1  lrs.c lrslib.c lrslong.c
	gcc -ansi  -O3 -DLONG -o redund1 redund.c lrslib.c lrslong.c
	gcc -ansi -O3 -o buffer buffer.c

lrs:    lrs.c lrslib.c lrslong.c lrsmp.c
	gcc -Wall -ansi -O3 -o lrs  lrs.c lrslib.c lrsmp.c

nash:    setupnash2.c setupnash.c nash.c lrslib.c  lrsmp.c
	gcc -Wall -DTIMES -ansi -O3 -o nash nash.c lrslib.c lrsmp.c
	gcc -Wall -o setupnash setupnash.c lrslib.c lrsmp.c
	gcc -Wall -o setupnash2 setupnash2.c lrslib.c lrsmp.c

fourier:    fourier.c lrslib.c lrslong.c lrsmp.c
	gcc -Wall -O3 -o fourier  fourier.c lrslib.c lrsmp.c
	gcc -O3 -static -DTIMES -DSIGNALS  -DGMP -I${INCLUDEDIR} fourier.c lrslib.c lrsgmp.c -L${LIBDIR}  -lgmp -o gfourier

demo:	lpdemo.c chdemo.c vedemo.c lrslib.c lrslong.c lrsmp.c
	gcc -Wall -ansi -O3 -o lpdemo lpdemo.c lrslib.c lrsmp.c
	gcc -Wall -ansi -O3 -o vedemo vedemo.c lrslib.c lrsmp.c
	gcc -Wall -ansi -O3 -o chdemo chdemo.c lrslib.c lrsmp.c

float: float2rat.c rat2float.c lrsmp.c 
	gcc -DLRSMP -Wall -ansi -o float2rat float2rat.c lrsmp.c
	gcc -DLRSMP -Wall -ansi -o rat2float rat2float.c lrsmp.c

clean:
	rm -rf lrs lrs1 redund redund1 buffer glrs gredund
	rm -rf foo gfoo
	rm -rf lpdemo vedemo chdemo
	rm -rf fourier gfourier                   
	rm -rf nash gnash setupnash setupnash2

foo:	foo.c lrslib.h lrslib.c lrsmp.h lrsmp.c
	gcc -O3 -static -DTIMES -DSIGNALS  foo.c lrslib.c lrsmp.c -L${LIBDIR} -o foo

gfoo:	foo.c lrslib.h lrslib.c lrsgmp.h lrsgmp.c
	gcc -O3 -static -DTIMES -DSIGNALS  -DGMP -I${INCLUDEDIR} foo.c lrslib.c lrsgmp.c -L${LIBDIR}  -lgmp -o gfoo

