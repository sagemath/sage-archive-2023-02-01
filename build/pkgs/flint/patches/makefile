ifndef FLINT_TUNE
	# defaults for sage.math development machine
	FLINT_TUNE = -mtune=opteron -march=opteron

	# for the record, here's what I use on my G5 powerpc:
	# FLINT_TUNE = -m64 -mcpu=970 -mtune=970 -mpowerpc64 -falign-loops=16 -falign-functions=16 -falign-labels=16 -falign-jumps=16

	# and here's for my laptop:
	# FLINT_TUNE =
endif

ifndef FLINT_LINK_OPTIONS
	FLINT_LINK_OPTIONS=
endif

# default GMP directories on sage.math development machine
ifndef FLINT_GMP_INCLUDE_DIR
	FLINT_GMP_INCLUDE_DIR = "/home/dmharvey/gmp/install/include"
endif

ifndef FLINT_GMP_LIB_DIR
	FLINT_GMP_LIB_DIR = "/home/dmharvey/gmp/install/lib"
endif

# qd include and library directories

ifndef FLINT_QD_INCLUDE_DIR
	FLINT_QD_INCLUDE_DIR = "/home/dmharvey/gmp/install/include"
endif

ifndef FLINT_QD_LIB_DIR
	FLINT_QD_LIB_DIR = "qd"
endif

#qdexists := $(shell ls -d qd)
#ifeq ($(qdexists), qd)
#	LIBS = -L$(FLINT_GMP_LIB_DIR) -L$(FLINT_QD_LIB_DIR) $(FLINT_LINK_OPTIONS) -lgmp -lpthread -lm -lqd
#else
        qdexists := 0
	LIBS = -L$(FLINT_GMP_LIB_DIR) -L$(FLINT_QD_LIB_DIR) $(FLINT_LINK_OPTIONS) -lgmp -lpthread -lm
#endif

INCS =  -I$(FLINT_GMP_INCLUDE_DIR) -I$(FLINT_QD_INCLUDE_DIR)

CC = gcc -std=c99

CPP = g++

CFLAGS = $(INCS) -funroll-loops -fexpensive-optimizations $(FLINT_TUNE) -O3

RM = rm -f

HEADERS = \
	Z.h \
	Z_mpn.h \
	Z_mpn_mul-tuning.h \
	ZmodF.h \
	ZmodF_mul-tuning.h \
	ZmodF_mul.h \
	ZmodF_poly.h \
	extras.h \
	flint.h \
	fmpz.h \
	fmpz_poly.h \
	longlong.h \
	longlong_wrapper.h \
	memory-manager.h \
	mpn_extras.h \
	mpz_poly-tuning.h \
	mpz_poly.h \
	profiler-main.h \
	profiler.h \
	test-support.h


####### library object files

FLINTOBJ = \
	mpn_extras.o \
	Z.o \
	memory-manager.o \
	Z_mpn.o \
	ZmodF.o \
	ZmodF_mul.o \
	ZmodF_mul-tuning.o \
	fmpz.o \
	fmpz_poly.o \
	mpz_poly-tuning.o \
	mpz_poly.o \
	ZmodF_poly.o


mpn_extras.o: mpn_extras.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpn_extras.c -o mpn_extras.o

Z.o: Z.c $(HEADERS)
	$(CC) $(CFLAGS) -c Z.c -o Z.o

memory-manager.o: memory-manager.c $(HEADERS)
	$(CC) $(CFLAGS) -c memory-manager.c -o memory-manager.o

Z_mpn.o: Z_mpn.c $(HEADERS)
	$(CC) $(CFLAGS) -c Z_mpn.c -o Z_mpn.o

ZmodF.o: ZmodF.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF.c -o ZmodF.o

ZmodF_mul.o: ZmodF_mul.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul.c -o ZmodF_mul.o

ZmodF_mul-tuning.o: ZmodF_mul-tuning.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul-tuning.c -o ZmodF_mul-tuning.o

fmpz.o: fmpz.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz.c -o fmpz.o

fmpz_poly.o: fmpz_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly.c -o fmpz_poly.o

mpz_poly.o: mpz_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly.c -o mpz_poly.o

mpz_poly-tuning.o: mpz_poly-tuning.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly-tuning.c -o mpz_poly-tuning.o

ZmodF_poly.o: ZmodF_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_poly.c -o ZmodF_poly.o

long_extras.o: long_extras.c long_extras.h
	$(CC) $(CFLAGS) -c long_extras.c -o long_extras.o



####### test program object files

test-support.o: test-support.c $(HEADERS)
	$(CC) $(CFLAGS) -c test-support.c -o test-support.o


fmpz_poly-test.o: fmpz_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly-test.c -o fmpz_poly-test.o

ZmodF-test.o: ZmodF-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF-test.c -o ZmodF-test.o

ZmodF_poly-test.o: ZmodF_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_poly-test.c -o ZmodF_poly-test.o

mpz_poly-test.o: mpz_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly-test.c -o mpz_poly-test.o

Z_mpn-test.o: Z_mpn-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c Z_mpn-test.c -o Z_mpn-test.o

ZmodF_mul-test.o: ZmodF_mul-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul-test.c -o ZmodF_mul-test.o

long_extras-test.o: long_extras-test.c
	$(CC) $(CFLAGS) -c long_extras-test.c -o long_extras-test.o



####### test program targets

Z_mpn-test: Z_mpn-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) Z_mpn-test.o test-support.o -o Z_mpn-test $(FLINTOBJ) $(LIBS)

fmpz_poly-test: fmpz_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) fmpz_poly-test.o test-support.o -o fmpz_poly-test $(FLINTOBJ) $(LIBS)

ZmodF-test: ZmodF-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) ZmodF-test.o test-support.o -o ZmodF-test $(FLINTOBJ) $(LIBS)

ZmodF_poly-test: ZmodF_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) ZmodF_poly-test.o test-support.o -o ZmodF_poly-test $(FLINTOBJ) $(LIBS)

mpz_poly-test: mpz_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) mpz_poly-test.o test-support.o -o mpz_poly-test $(FLINTOBJ) $(LIBS)

ZmodF_mul-test: ZmodF_mul-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) ZmodF_mul-test.o test-support.o -o ZmodF_mul-test $(FLINTOBJ) $(LIBS)

long_extras-test: long_extras.o long_extras-test.o test-support.o memory-manager.o
	$(CC) $(CFLAGS) long_extras.o long_extras-test.o test-support.o memory-manager.o -o long_extras-test $(LIBS)


####### tuning program object files

ZmodF_mul-tune.o: ZmodF_mul-tune.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul-tune.c -o ZmodF_mul-tune.o

mpz_poly-tune.o: mpz_poly-tune.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly-tune.c -o mpz_poly-tune.o



####### tuning program targets

ZmodF_mul-tune: ZmodF_mul-tune.o test-support.o profiler.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) ZmodF_mul-tune.o test-support.o profiler.o -o ZmodF_mul-tune $(FLINTOBJ) $(LIBS)

mpz_poly-tune: mpz_poly-tune.o test-support.o profiler.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) mpz_poly-tune.o test-support.o profiler.o -o mpz_poly-tune $(FLINTOBJ) $(LIBS)


####### profiling object files

profiler.o: profiler.c $(HEADERS)
	$(CC) $(CFLAGS) -c profiler.c -o profiler.o

profiler-main.o: profiler-main.c $(HEADERS)
	$(CC) $(CFLAGS) -c profiler-main.c -o profiler-main.o


fmpz_poly-profile-tables.o: fmpz_poly-profile.c $(HEADERS)
	python make-profile-tables.py fmpz_poly
	$(CC) $(CFLAGS) -c fmpz_poly-profile-tables.c -o fmpz_poly-profile-tables.o
	rm fmpz_poly-profile-tables.c

fmpz_poly-profile.o: fmpz_poly-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly-profile.c -o fmpz_poly-profile.o


mpz_poly-profile-tables.o: mpz_poly-profile.c $(HEADERS)
	python make-profile-tables.py mpz_poly
	$(CC) $(CFLAGS) -c mpz_poly-profile-tables.c -o mpz_poly-profile-tables.o
	rm mpz_poly-profile-tables.c

mpz_poly-profile.o: mpz_poly-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly-profile.c -o mpz_poly-profile.o


ZmodF_poly-profile-tables.o: ZmodF_poly-profile.c $(HEADERS)
	python make-profile-tables.py ZmodF_poly
	$(CC) $(CFLAGS) -c ZmodF_poly-profile-tables.c -o ZmodF_poly-profile-tables.o
	rm ZmodF_poly-profile-tables.c

ZmodF_poly-profile.o: ZmodF_poly-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_poly-profile.c -o ZmodF_poly-profile.o


ZmodF_mul-profile-tables.o: ZmodF_mul-profile.c $(HEADERS)
	python make-profile-tables.py ZmodF_mul
	$(CC) $(CFLAGS) -c ZmodF_mul-profile-tables.c -o ZmodF_mul-profile-tables.o
	rm ZmodF_mul-profile-tables.c

ZmodF_mul-profile.o: ZmodF_mul-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul-profile.c -o ZmodF_mul-profile.o




####### profiling program targets

PROFOBJ = $(FLINTOBJ) profiler.o profiler-main.o

fmpz_poly-profile: fmpz_poly-profile.o fmpz_poly-profile-tables.o test-support.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o fmpz_poly-profile fmpz_poly-profile.o fmpz_poly-profile-tables.o test-support.o $(PROFOBJ) $(LIBS)

mpz_poly-profile: mpz_poly-profile.o mpz_poly-profile-tables.o test-support.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o mpz_poly-profile mpz_poly-profile.o mpz_poly-profile-tables.o test-support.o $(PROFOBJ) $(LIBS)

ZmodF_mul-profile: ZmodF_mul-profile.o ZmodF_mul-profile-tables.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o ZmodF_mul-profile ZmodF_mul-profile.o ZmodF_mul-profile-tables.o $(PROFOBJ) $(LIBS)

ZmodF_poly-profile: ZmodF_poly-profile.o ZmodF_poly-profile-tables.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o ZmodF_poly-profile ZmodF_poly-profile.o ZmodF_poly-profile-tables.o $(PROFOBJ) $(LIBS)


kara-profile: kara-profile.c profiler.o test-support.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o kara-profile kara-profile.c profiler.o test-support.o $(FLINTOBJ) $(LIBS)



####### example programs

delta_qexp.o: delta_qexp.c $(HEADERS)
	$(CC) $(CFLAGS) -c delta_qexp.c -o delta_qexp.o

delta_qexp: delta_qexp.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o delta_qexp delta_qexp.o $(FLINTOBJ) $(LIBS)

BLTcubes: long_extras.o BLTcubes.c
	$(CC) $(CFLAGS) -o BLTcubes BLTcubes.c long_extras.o $(LIBS)

BPTJCubes: long_extras.o memory-manager.o
	$(CC) $(CFLAGS) -o BPTJCubes BPTJCubes.c memory-manager.o long_extras.o $(LIBS)

bernoulli.o: bernoulli.c $(HEADERS)
	$(CC) $(CFLAGS) -c bernoulli.c -o bernoulli.o

bernoulli: bernoulli.o long_extras.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o bernoulli bernoulli.o long_extras.o $(FLINTOBJ) $(LIBS)


####### Quadratic sieve

poly.o: QS/poly.c QS/poly.h
	$(CC) $(CFLAGS) -c QS/poly.c -o poly.o

factor_base.o: QS/factor_base.c QS/factor_base.h
	$(CC) $(CFLAGS) -c QS/factor_base.c -o factor_base.o

sieve.o: QS/sieve.c QS/sieve.h
	$(CC) $(CFLAGS) -c QS/sieve.c -o sieve.o

linear_algebra.o: QS/linear_algebra.c QS/linear_algebra.h
	$(CC) $(CFLAGS) -c QS/linear_algebra.c -o linear_algebra.o

block_lanczos.o: QS/block_lanczos.c QS/block_lanczos.h
	$(CC) $(CFLAGS) -c QS/block_lanczos.c -o block_lanczos.o

tinyQS: QS/tinyQS.c QS/tinyQS.h factor_base.o poly.o sieve.o linear_algebra.o block_lanczos.o long_extras.o memory-manager.o fmpz.o test-support.o
	$(CC) $(CFLAGS) -o tinyQS QS/tinyQS.c factor_base.o poly.o sieve.o linear_algebra.o block_lanczos.o memory-manager.o long_extras.o fmpz.o test-support.o $(LIBS)

mp_sieve.o: QS/mp_sieve.c QS/mp_sieve.h
	$(CC) $(CFLAGS) -c QS/mp_sieve.c -o mp_sieve.o

mp_linear_algebra.o: QS/mp_linear_algebra.c QS/mp_linear_algebra.h
	$(CC) $(CFLAGS) -c QS/mp_linear_algebra.c -o mp_linear_algebra.o

mp_poly.o: QS/mp_poly.c QS/mp_poly.h
	$(CC) $(CFLAGS) -c QS/mp_poly.c -o mp_poly.o

mp_lprels.o: QS/mp_lprels.c QS/mp_lprels.h
	$(CC) $(CFLAGS) -c QS/mp_lprels.c -o mp_lprels.o

mp_factor_base.o: QS/mp_factor_base.c QS/mp_factor_base.h
	$(CC) $(CFLAGS) -c QS/mp_factor_base.c -o mp_factor_base.o

mpQS: QS/mpQS.c QS/mpQS.h mp_factor_base.o mp_poly.o mp_sieve.o mp_linear_algebra.o block_lanczos.o mp_lprels.o long_extras.o memory-manager.o fmpz.o test-support.o mpn_extras.o
	$(CC) $(CFLAGS) -o mpQS QS/mpQS.c mp_factor_base.o mp_poly.o mp_sieve.o mp_linear_algebra.o block_lanczos.o mp_lprels.o memory-manager.o long_extras.o fmpz.o test-support.o mpn_extras.o $(LIBS)

####### Integer multiplication timing

ZMULOBJ = memory-manager.o fmpz.o ZmodF_mul-tuning.o fmpz_poly.o ZmodF_poly.o Z_mpn.o profiler.o ZmodF_mul.o ZmodF.o mpn_extras.o Z_mul_timing.o

Z_mul_timing: $(ZMULOBJ)
	$(CC) $(ZMULOBJ) -o Zmul $(LIBS)

####### Linear Algebra

memory-manager2.o: memory-manager.c $(HEADERS)
	$(CPP) $(CFLAGS) -c memory-manager.c -o memory-manager2.o

vec3d.o: vec3d.c vec3d.h
	$(CPP) $(CFLAGS) -c vec3d.c -o vec3d.o

mat3d.o: mat3d.c mat3d.h
	$(CPP) $(CFLAGS) -c mat3d.c -o mat3d.o

vecmat3d: vecmat3d-driver.c mat3d.o vec3d.o memory-manager2.o
	$(CPP) $(CFLAGS) -o vecmat3d vecmat3d-driver.c mat3d.o vec3d.o memory-manager2.o $(LIBS)

x3y3z3k: x3y3z3k.c mat3d.o vec3d.o memory-manager2.o
	$(CPP) -o x3y3z3k $(CFLAGS) x3y3z3k.c mat3d.o vec3d.o memory-manager2.o $(LIBS)

dd_vecmat3d: dd_vecmat3d-driver.c mat3d.o vec3d.o memory-manager2.o
	$(CPP) $(CFLAGS) -o dd_vecmat3d dd_vecmat3d-driver.c mat3d.o vec3d.o memory-manager2.o $(LIBS)

dd_x3y3z3k: dd_x3y3z3k.c mat3d.o vec3d.o memory-manager2.o
	$(CPP) -o dd_x3y3z3k $(CFLAGS) dd_x3y3z3k.c mat3d.o vec3d.o memory-manager2.o $(LIBS)


