#ifndef DGS__H
#define DGS__H

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

/**
 * \brief Macro to help with branch prediction.
 */

#define __DGS_LIKELY(cond)    __builtin_expect ((cond) != 0, 1)

/**
 * \brief Macro to help with branch prediction.
 */

#define __DGS_UNLIKELY(cond)  __builtin_expect ((cond) != 0, 0)

/*
 * \brief balanced Bernoulli distribution, multi-precision version
 *
 * This implementation is faster than dgs_bern_mp_t
 */

/**
 * \brief Number of bits sampled at once
 */

#define DGS_BERN_UNIFORM_DEFAULT_LENGTH (sizeof(unsigned long)*8)

/**
 * \brief Maximum number of bits sampled at once
 */

#define DGS_BERN_UNIFORM_MAX_LENGTH (sizeof(unsigned long)*8)

typedef struct {
  size_t length; //< number of bits we sample in each go
  size_t count; //< number of bits we have consumed yet
  mpz_t tmp; //< we sample to this mpz_t
  uint64_t pool; //< we store the pool of random bits here
} dgs_bern_uniform_mp_t;

dgs_bern_uniform_mp_t* dgs_bern_uniform_mp_init(size_t length);
unsigned long dgs_bern_uniform_mp_call(dgs_bern_uniform_mp_t *self, gmp_randstate_t state);
void dgs_bern_uniform_mp_clear(dgs_bern_uniform_mp_t *self);


/**
 * \brief Bernoulli distribution, multi-precision version
 */

typedef struct {
  mpfr_t p;
  mpfr_t tmp;
} dgs_bern_mp_t;

dgs_bern_mp_t *dgs_bern_mp_init(mpfr_t p);
unsigned long dgs_bern_mp_call(dgs_bern_mp_t *self, gmp_randstate_t state);
void dgs_bern_mp_clear(dgs_bern_mp_t *self);

/*
 * Bernoulli with p = exp(-x/f) for integers x, multi-precision version
 */

#define DGS_BERN_EXP_ALLOC_BLOCK_SIZE 16

typedef struct {
  size_t l;
  mpfr_t *p;
  dgs_bern_mp_t **B;
} dgs_bern_exp_mp_t;

dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l);
unsigned long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x, gmp_randstate_t state);
void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self);

/**
 * Discrete Gaussians
 */

typedef enum {
  DGS_DISC_GAUSS_TABLE   = 0x1,
  DGS_DISC_GAUSS_UNIFORM = 0x2,
} dgs_disc_gauss_alg_t;

struct _dgs_disc_gauss_mp_t;

typedef struct _dgs_disc_gauss_mp_t {
  mpfr_t sigma;
  size_t tailcut;
  dgs_disc_gauss_alg_t algorithm;

  dgs_bern_uniform_mp_t *B;
  
  void (*call)(mpz_t rop, struct _dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

  mpz_t upper_bound;
  mpz_t x;
  mpfr_t y;
  mpfr_t *rho;
  
} dgs_disc_gauss_mp_t;

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, size_t tailcut, dgs_disc_gauss_alg_t algorithm);
void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self);

#endif //DGS__H
