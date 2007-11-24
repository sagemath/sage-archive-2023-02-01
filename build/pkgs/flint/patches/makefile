ifndef FLINT_TUNE
	# defaults for sage.math development machine
	FLINT_TUNE = -march=opteron -mtune=opteron

	# for the record, here's what I use on my G5 powerpc:
	# FLINT_TUNE = -m64 -mcpu=970 -mtune=970 -mpowerpc64 -falign-loops=16 -falign-functions=16 -falign-labels=16 -falign-jumps=16

	# and here's for my laptop:
	# FLINT_TUNE =
endif

ifndef NO_NTL
   NTL = -lntl
endif

ifndef FLINT_LINK_OPTIONS
	FLINT_LINK_OPTIONS =
endif


qdexists := $(ls -d qd 2> /dev/null)
ifeq ($(qdexists), qd)
	LIBS = -L$(FLINT_GMP_LIB_DIR)  -L$(FLINT_NTL_LIB_DIR) -L$(FLINT_QD_LIB_DIR) $(FLINT_LINK_OPTIONS) -lgmp -lpthread -lm -lqd $(NTL)
else
	LIBS = -L$(FLINT_GMP_LIB_DIR)  -L$(FLINT_NTL_LIB_DIR) -L$(FLINT_QD_LIB_DIR) $(FLINT_LINK_OPTIONS) -lgmp -lpthread -lm $(NTL)
endif

all: mpn_extras-test fmpz_poly-test fmpz-test ZmodF-test ZmodF_poly-test mpz_poly-test ZmodF_mul-test long_extras-test ZmodF_mul-tune mpz_poly-tune fmpz_poly-profile mpz_poly-profile ZmodF_mul-profile ZmodF_poly-profile kara-profile delta_qexp BLTcubes BPTJCubes bernoulli F_mpz_mul-timing expmod mpQS tinyQS

INCS = -I$(FLINT_GMP_INCLUDE_DIR) -I$(FLINT_NTL_INCLUDE_DIR) -I$(FLINT_QD_INCLUDE_DIR)

CC = gcc -std=c99

CPP = g++

CFLAGS = $(INCS) -funroll-loops -fexpensive-optimizations $(FLINT_TUNE) -O3

RM = rm -f

HEADERS = \
	mpz_extras.h \
	F_mpn_mul-tuning.h \
	ZmodF.h \
	ZmodF_mul-tuning.h \
	ZmodF_mul.h \
	ZmodF_poly.h \
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
	test-support.h \
	long_extras.h


####### library object files

FLINTOBJ = \
	mpn_extras.o \
	mpz_extras.o \
	memory-manager.o \
	ZmodF.o \
	ZmodF_mul.o \
	ZmodF_mul-tuning.o \
	fmpz.o \
	fmpz_poly.o \
	mpz_poly-tuning.o \
	mpz_poly.o \
	ZmodF_poly.o \
	long_extras.o

mpn_extras.o: mpn_extras.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpn_extras.c -o mpn_extras.o

mpz_extras.o: mpz_extras.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_extras.c -o mpz_extras.o

memory-manager.o: memory-manager.c $(HEADERS)
	$(CC) $(CFLAGS) -c memory-manager.c -o memory-manager.o

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

zmod_poly.o: zmod_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c zmod_poly.c -o zmod_poly.o


####### test program object files

test-support.o: test-support.c $(HEADERS)
	$(CC) $(CFLAGS) -c test-support.c -o test-support.o

fmpz_poly-test.o: fmpz_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly-test.c -o fmpz_poly-test.o

fmpz-test.o: fmpz-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz-test.c -o fmpz-test.o

ZmodF-test.o: ZmodF-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF-test.c -o ZmodF-test.o

ZmodF_poly-test.o: ZmodF_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_poly-test.c -o ZmodF_poly-test.o

mpz_poly-test.o: mpz_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly-test.c -o mpz_poly-test.o

mpn_extras-test.o: mpn_extras-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpn_extras-test.c -o mpn_extras-test.o

ZmodF_mul-test.o: ZmodF_mul-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul-test.c -o ZmodF_mul-test.o

long_extras-test.o: long_extras-test.c
	$(CC) $(CFLAGS) -c long_extras-test.c -o long_extras-test.o

zmod_poly-test.o: zmod_poly-test.c
	$(CC) $(CFLAGS) -c zmod_poly-test.c -o zmod_poly-test.o

####### test program targets

mpn_extras-test: mpn_extras-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) mpn_extras-test.o test-support.o -o mpn_extras-test $(FLINTOBJ) $(LIBS)

fmpz_poly-test: fmpz_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) fmpz_poly-test.o test-support.o -o fmpz_poly-test $(FLINTOBJ) $(LIBS)

fmpz-test: fmpz-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) fmpz-test.o test-support.o -o fmpz-test $(FLINTOBJ) $(LIBS)

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

zmod_poly-test: zmod_poly.o zmod_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) zmod_poly.o zmod_poly-test.o test-support.o -o zmod_poly-test $(FLINTOBJ) $(LIBS)


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

NTL-profile-tables.o: NTL-profile.c $(HEADERS)
	python make-profile-tables.py NTL
	$(CPP) $(CFLAGS) -c NTL-profile-tables.c -o NTL-profile-tables.o

zmod_poly-profile-tables.o: zmod_poly-profile.c $(HEADERS)
	python make-profile-tables.py zmod_poly
	$(CC) $(CFLAGS) -c zmod_poly-profile-tables.c -o zmod_poly-profile-tables.o
	rm zmod_poly-profile-tables.c

zmod_poly-profile.o: zmod_poly-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c zmod_poly-profile.c -o zmod_poly-profile.o

bernoulli-profile-tables.o: bernoulli-profile.c $(HEADERS)
	python make-profile-tables.py bernoulli
	$(CC) $(CFLAGS) -c bernoulli-profile-tables.c -o bernoulli-profile-tables.o
	rm bernoulli-profile-tables.c

bernoulli-profile.o: bernoulli-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c bernoulli-profile.c -o bernoulli-profile.o

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

NTL-profile: NTL-profile.c test-support.o NTL-profile-tables.o $(PROFOBJ)
	$(CPP) $(CFLAGS) -o NTL-profile NTL-profile.c NTL-profile-tables.o test-support.o $(PROFOBJ) $(LIBS)

zmod_poly-profile: zmod_poly-profile.o zmod_poly.o zmod_poly-profile-tables.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o zmod_poly-profile zmod_poly.o zmod_poly-profile.o zmod_poly-profile-tables.o $(PROFOBJ) $(LIBS)

bernoulli-profile: bernoulli-profile.o zmod_poly.o bernoulli-profile-tables.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o bernoulli-profile zmod_poly.o bernoulli-profile.o bernoulli-profile-tables.o $(PROFOBJ) $(LIBS)

####### example programs

delta_qexp.o: delta_qexp.c $(HEADERS)
	$(CC) $(CFLAGS) -c delta_qexp.c -o delta_qexp.o

delta_qexp: delta_qexp.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o delta_qexp delta_qexp.o $(FLINTOBJ) $(LIBS)

BLTcubes: long_extras.o BLTcubes.c memory-manager.o
	$(CC) $(CFLAGS) -o BLTcubes BLTcubes.c long_extras.o memory-manager.o $(LIBS)

BPTJCubes: long_extras.o memory-manager.o
	$(CC) $(CFLAGS) -o BPTJCubes BPTJCubes.c memory-manager.o long_extras.o $(LIBS)

bernoulli.o: bernoulli.c $(HEADERS)
	$(CC) $(CFLAGS) -c bernoulli.c -o bernoulli.o

bernoulli: bernoulli.o long_extras.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o bernoulli bernoulli.o $(FLINTOBJ) $(LIBS)

bernoulli_fmpz.o: bernoulli_fmpz.c $(HEADERS)
	$(CC) $(CFLAGS) -c bernoulli_fmpz.c -o bernoulli_fmpz.o

bernoulli_fmpz: bernoulli_fmpz.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o bernoulli_fmpz bernoulli_fmpz.o $(FLINTOBJ) $(LIBS)

bernoulli_zmod.o: bernoulli_zmod.c $(HEADERS)
	$(CC) $(CFLAGS) -c bernoulli_zmod.c -o bernoulli_zmod.o

bernoulli_zmod: bernoulli_zmod.o zmod_poly.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o bernoulli_zmod bernoulli_zmod.o zmod_poly.o $(FLINTOBJ) $(LIBS)

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

tinyQS: QS/tinyQS.c QS/tinyQS.h factor_base.o poly.o sieve.o linear_algebra.o block_lanczos.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o tinyQS QS/tinyQS.c factor_base.o poly.o sieve.o linear_algebra.o block_lanczos.o $(FLINTOBJ) $(LIBS)

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

mpQS: QS/mpQS.c QS/mpQS.h mp_factor_base.o mp_poly.o mp_sieve.o mp_linear_algebra.o block_lanczos.o mp_lprels.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o mpQS QS/mpQS.c mp_factor_base.o mp_poly.o mp_sieve.o mp_linear_algebra.o block_lanczos.o mp_lprels.o $(FLINTOBJ) $(LIBS)

####### Integer multiplication timing

ZMULOBJ = memory-manager.o fmpz.o ZmodF_mul-tuning.o mpz_poly.o mpz_poly-tuning.o fmpz_poly.o ZmodF_poly.o mpz_extras.o profiler.o ZmodF_mul.o ZmodF.o mpn_extras.o F_mpz_mul-timing.o

F_mpz_mul-timing: $(ZMULOBJ)
	$(CC) $(ZMULOBJ) -o Zmul $(LIBS)

####### Linear Algebra

vec3d.o: vec3d.c vec3d.h
	$(CC) $(CFLAGS) -c vec3d.c -o vec3d.o

mat3d.o: mat3d.c mat3d.h
	$(CC) $(CFLAGS) -c mat3d.c -o mat3d.o

vecmat3d: vecmat3d-driver.c mat3d.o vec3d.o memory-manager.o
	$(CC) $(CFLAGS) -o vecmat3d vecmat3d-driver.c mat3d.o vec3d.o memory-manager.o $(LIBS)

x3y3z3k: x3y3z3k.c mat3d.o vec3d.o memory-manager.o
	$(CC) -o x3y3z3k $(CFLAGS) x3y3z3k.c mat3d.o vec3d.o memory-manager.o $(LIBS)

dd_vecmat3d: dd_vecmat3d-driver.c mat3d.o vec3d.o memory-manager.o
	$(CC) $(CFLAGS) -o dd_vecmat3d dd_vecmat3d-driver.c mat3d.o vec3d.o memory-manager.o $(LIBS)

dd_x3y3z3k: dd_x3y3z3k.c mat3d.o vec3d.o memory-manager.o
	$(CC) -o dd_x3y3z3k $(CFLAGS) dd_x3y3z3k.c mat3d.o vec3d.o memory-manager.o $(LIBS)

expmod: expmod.c $(FLINTOBJ)
	$(CC) $(CFLAGS) -o expmod expmod.c $(FLINTOBJ) $(LIBS)

