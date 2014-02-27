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

#define DGS_BERN_UNIFORM_DEFAULT_LENGTH (sizeof(unsigned long)*8)
#define DGS_BERN_UNIFORM_MAX_LENGTH (sizeof(unsigned long)*8)

/**
 * uniformly random bits: sample 1 with probability 0.5 (faster than dgs_bern_t)
 */

typedef struct {
  size_t length;
  size_t count;
  mpz_t tmp;
  uint64_t pool;
} dgs_bern_uniform_t;

dgs_bern_uniform_t* dgs_bern_uniform_init(size_t length);
unsigned long dgs_bern_uniform_call(dgs_bern_uniform_t *self, gmp_randstate_t state);
void dgs_bern_uniform_clear(dgs_bern_uniform_t *self);


/**
 * Bernoulli distribution: sample 1 with probability c
 */

typedef struct {
  mpfr_t c;
  mpfr_t tmp;
} dgs_bern_t;

dgs_bern_t *dgs_bern_init(mpfr_t c);
unsigned long dgs_bern_call(dgs_bern_t *self, gmp_randstate_t state);
void dgs_bern_clear(dgs_bern_t *self);

/**
 * Sample 1 with probability exp(-x/f)
 */

#define DGS_BERN_EXP_ALLOC_BLOCK_SIZE 16

typedef struct {
  size_t l;
  mpfr_t *c;
  dgs_bern_t **B;
} dgs_bern_exp_t;

dgs_bern_exp_t* dgs_bern_exp_init(mpfr_t f, size_t l);
unsigned long dgs_bern_exp_call(dgs_bern_exp_t *self, mpz_t x, gmp_randstate_t state);
void dgs_bern_exp_clear(dgs_bern_exp_t *self);

