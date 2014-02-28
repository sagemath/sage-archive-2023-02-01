#include "dgs.h"
#include <assert.h>
#include <stdlib.h>

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, size_t tailcut, dgs_disc_gauss_alg_t algorithm) {
  assert(tailcut > 0);

  dgs_disc_gauss_mp_t *self = (dgs_disc_gauss_mp_t*)malloc(sizeof(dgs_disc_gauss_mp_t));
  if (!self) abort();
  
  mpfr_init_set(self->sigma, sigma, MPFR_RNDN);
  self->tailcut = tailcut;

  if(algorithm == (DGS_DISC_GAUSS_TABLE & DGS_DISC_GAUSS_UNIFORM)) {
    self->call = dgs_disc_gauss_mp_call_uniform_table;
  } else {
    self->call = dgs_disc_gauss_mp_call_uniform_table;
  }

  self->B = dgs_bern_uniform_mp_init(0);
  
  if (self->call == dgs_disc_gauss_mp_call_uniform_table) {
    mpfr_t f;
    mpfr_init_set(f, sigma, MPFR_RNDN); // f = σ
    mpfr_sqr(f, f, MPFR_RNDN); // f = σ^2
    mpfr_mul_ui(f, f, 2, MPFR_RNDN); // f = 2 σ^2
    mpfr_ui_div(f, 1, f, MPFR_RNDN);  // f = 1/(2 σ^2)
    mpfr_neg(f, f, MPFR_RNDN); // f = -1/(2 σ^2)

    mpfr_t tmp;
    mpfr_init(tmp);
    mpz_init(self->upper_bound);
    mpfr_mul_ui(tmp, sigma, tailcut, MPFR_RNDN);
    mpfr_get_z(self->upper_bound, tmp, MPFR_RNDU);
    mpfr_clear(tmp);

    self->rho = (mpfr_t*)malloc(sizeof(mpfr_t)*mpz_get_ui(self->upper_bound));
    if (!self->rho) abort();
    mpfr_t x_;
    mpfr_init(x_);
    for(unsigned long x=0; x<mpz_get_ui(self->upper_bound); x++) {
      mpfr_set_ui(x_, x, MPFR_RNDN);
      mpfr_sqr(x_, x_, MPFR_RNDN);
      mpfr_mul(x_, x_, f, MPFR_RNDN);
      mpfr_exp(x_, x_, MPFR_RNDN);
      mpfr_init_set(self->rho[x], x_, MPFR_RNDN);
    }
    mpfr_div_ui(self->rho[0],self->rho[0],2, MPFR_RNDN);
    mpfr_clear(x_);
    mpfr_clear(f);

    mpz_init(self->x);
    mpfr_init(self->y);
  }
  return self;
}

void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  unsigned long x;
  do {
    mpz_urandomm(self->x, state, self->upper_bound);
    x = mpz_get_ui(self->x);
    mpfr_urandomb(self->y, state);
  } while (mpfr_cmp(self->y, self->rho[x]) >= 0);

  mpz_init_set_ui(rop, x);
  if(dgs_bern_uniform_mp_call(self->B, state))
    mpz_neg(rop, rop);
}


void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self) {
  mpfr_clear(self->sigma);
  if (self->B) dgs_bern_uniform_mp_clear(self->B);
  if (self->x) mpz_clear(self->x);
  if (self->y) mpfr_clear(self->y);
  if (self->rho) {
    for(unsigned long x=0; x<mpz_get_ui(self->upper_bound); x++) {
      mpfr_clear(self->rho[x]);
    }
    free(self->rho);
  }
  free(self);
}
