#include "dgs.h"
#include <assert.h>
#include <stdlib.h>

/* uniformly random bits */

dgs_bern_uniform_t* dgs_bern_uniform_init(size_t length) {
  if (length == 0) {
    length = DGS_BERN_UNIFORM_DEFAULT_LENGTH;
  }
  assert(length <= DGS_BERN_UNIFORM_MAX_LENGTH);

  dgs_bern_uniform_t *self = malloc(sizeof(dgs_bern_uniform_t));
  self->length = length;

  self->count = self->length;
  mpz_init(self->tmp);
  return self;
}

unsigned long dgs_bern_uniform_call(dgs_bern_uniform_t *self, gmp_randstate_t state) {
  assert(self != NULL);
  assert(state != NULL);

  if (__DGS_UNLIKELY(self->count == self->length)) {
    mpz_urandomb(self->tmp, state, self->length);
    self->pool = mpz_get_ui(self->tmp);
    self->count = 0;
  }

  unsigned long b = self->pool & 1;
  self->pool >>= 1;
  self->count++;
  return b;
}

void dgs_bern_uniform_clear(dgs_bern_uniform_t *self) {
  mpz_clear(self->tmp);
  free(self);
}

/* Bernoulli distribution */

dgs_bern_t* dgs_bern_init(mpfr_t c) {
  /* we allow 0 here for low precision */
  assert((mpfr_cmp_d(c, 0.0) >= 0) && (mpfr_cmp_d(c, 1.0) < 0));

  dgs_bern_t *self = malloc(sizeof(dgs_bern_t));

  mpfr_init_set(self->c, c, MPFR_RNDN);
  mpfr_init2(self->tmp, mpfr_get_prec(c));
  return self;
}

unsigned long dgs_bern_call(dgs_bern_t *self, gmp_randstate_t state) {
  int r = mpfr_urandomb(self->tmp, state);
  assert(r == 0);

  if (mpfr_cmp(self->tmp, self->c)<0) {
    return 1;
  } else {
    return 0;
  }
}

void dgs_bern_clear(dgs_bern_t *self) {
  mpfr_clear(self->tmp);
  mpfr_clear(self->c);
  free(self);
}

dgs_bern_exp_t* dgs_bern_exp_init(mpfr_t f, size_t l) {
  dgs_bern_exp_t *self = malloc(sizeof(dgs_bern_exp_t));

  /* l == 0, means we use the precision of f to decide l */
  if (l == 0)
    l = SIZE_MAX;
 
  self->l = DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
  self->c = malloc(sizeof(mpfr_t)*self->l); 
  self->B = malloc(sizeof(dgs_bern_exp_t)*self->l);
 
  mpfr_t tmp, tmp2;
  mpfr_init2(tmp2, mpfr_get_prec(f));
  mpfr_init_set(tmp, f, MPFR_RNDN); // f
  mpfr_pow_si(tmp, tmp, -1, MPFR_RNDN); // 1/f
  mpfr_neg(tmp, tmp, MPFR_RNDN); // -1/f

  for(size_t i=0; i<l; i++) {
    mpfr_exp(tmp2, tmp, MPFR_RNDN);
    if (mpfr_zero_p(tmp2)) {
      self->l = i+1;
      break;
    }
    if (i%DGS_BERN_EXP_ALLOC_BLOCK_SIZE == 0 && i!=0) {
      self->l += DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
      self->l = (l>self->l) ? self->l : l;
      self->c = realloc(self->c, sizeof(mpfr_t)*self->l);
      self->B = realloc(self->B, sizeof(dgs_bern_exp_t)*self->l);
    }
    
    mpfr_init_set(self->c[i], tmp2, MPFR_RNDN);
    self->B[i] = dgs_bern_init(self->c[i]);

    mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
  }
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
  return self;
}

unsigned long dgs_bern_exp_call(dgs_bern_exp_t *self, mpz_t x, gmp_randstate_t state) {
  assert(mpz_sgn(x) == 1);
  long int start = (mpz_sizeinbase(x, 2) < self->l) ? mpz_sizeinbase(x, 2) : self->l;

  for(long int i=start-1; i>=0; i--) {
    if (mpz_tstbit(x, i)) {
      if (dgs_bern_call(self->B[i], state) == 0) {
        return 0;
      }
    }
  }
  return 1;
}

void dgs_bern_exp_clear(dgs_bern_exp_t *self) {
  for(size_t i=0; i<self->l; i++) {
    mpfr_clear(self->c[i]);
    dgs_bern_clear(self->B[i]);
  }
  free(self->c);
  free(self->B);
  free(self);
}


