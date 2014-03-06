#include "dgs.h"
#include <assert.h>
#include <stdlib.h>

dgs_disc_gauss_sigma2p_t *dgs_disc_gauss_sigma2p_init() {
  dgs_disc_gauss_sigma2p_t *self = (dgs_disc_gauss_sigma2p_t*)calloc(sizeof(dgs_disc_gauss_sigma2p_t),1);
  if (!self) abort();
  self->B = dgs_bern_uniform_init(0);
  return self;
}

void dgs_disc_gauss_sigma2p_mp_call(mpz_t rop, dgs_disc_gauss_sigma2p_t *self, gmp_randstate_t state) {
  while(1) {
    if (!dgs_bern_uniform_call(self->B, state)) {
      mpz_set_ui(rop, 0);
      return;
    }
    int dobreak = 0;
    for(unsigned long i=1; ;i++) {
      for(size_t j=0; j<2*i-2; j++) {
        if(dgs_bern_uniform_call(self->B, state)) {
          dobreak = 1;
          break;
        }
      }
      if (__DGS_LIKELY(dobreak))
        break;
      if (!dgs_bern_uniform_call(self->B, state)) {
        mpz_set_ui(rop, i);
        return;
      }
    }
  }
}

long dgs_disc_gauss_sigma2p_dp_call(dgs_disc_gauss_sigma2p_t *self, gmp_randstate_t state) {
  while(1) {
    if (!dgs_bern_uniform_call(self->B, state)) {
      return 0;
    }
    int dobreak = 0;
    for(unsigned long i=1; ;i++) {
      for(size_t j=0; j<2*i-2; j++) {
        if(dgs_bern_uniform_call(self->B, state)) {
          dobreak = 1;
          break;
        }
      }
      if (__DGS_LIKELY(dobreak))
        break;
      if (!dgs_bern_uniform_call(self->B, state)) {
        return i;
      }
    }
  }
}


void dgs_disc_gauss_sigma2p_clear(dgs_disc_gauss_sigma2p_t *self) {
  assert(self != NULL);
  if (self->B) dgs_bern_uniform_clear(self->B);
  free(self);
}


static inline void _dgs_disc_gauss_mp_init_f(mpfr_t f, const mpfr_t sigma) {
  mpfr_init2(f, mpfr_get_prec(sigma));
  mpfr_set(f, sigma, MPFR_RNDN);
  mpfr_sqr(f, f, MPFR_RNDN); // f = σ^2
  mpfr_mul_ui(f, f, 2, MPFR_RNDN); // f = 2 σ^2
  mpfr_ui_div(f, 1, f, MPFR_RNDN);  // f = 1/(2 σ^2)
  mpfr_neg(f, f, MPFR_RNDN); // f = -1/(2 σ^2)
}

static inline void _dgs_disc_gauss_mp_init_upper_bound(mpz_t upper_bound, mpz_t two_upper_bound_plus_one,
                                                       mpfr_t sigma, size_t tailcut) {
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(sigma));
  mpz_init(upper_bound);
  mpz_init(two_upper_bound_plus_one);
  mpfr_mul_ui(tmp, sigma, tailcut, MPFR_RNDN);
  mpfr_get_z(upper_bound, tmp, MPFR_RNDU);
  mpz_mul_ui(two_upper_bound_plus_one, upper_bound, 2);
  mpz_add_ui(two_upper_bound_plus_one, two_upper_bound_plus_one, 1);
  mpfr_clear(tmp);
}

static inline void _dgs_disc_gauss_mp_init_bexp(dgs_disc_gauss_mp_t *self, mpfr_t sigma, mpz_t upper_bound) {
  mpfr_init2(self->f, mpfr_get_prec(sigma));
  mpfr_set(self->f, sigma, MPFR_RNDN); // f = σ
  mpfr_sqr(self->f, self->f, MPFR_RNDN); // f = σ^2
  mpfr_mul_ui(self->f, self->f, 2, MPFR_RNDN); // f = 2 σ^2
  size_t l = 2*mpz_sizeinbase(upper_bound, 2);
  self->Bexp = dgs_bern_exp_mp_init(self->f, l);
}

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, size_t tailcut, dgs_disc_gauss_alg_t algorithm) {
  assert(tailcut > 0);

  mpfr_prec_t prec = mpfr_get_prec(sigma);
  
  dgs_disc_gauss_mp_t *self = (dgs_disc_gauss_mp_t*)calloc(sizeof(dgs_disc_gauss_mp_t),1);
  if (!self) abort();

  mpz_init(self->x);
  mpz_init(self->x2);
  mpfr_init2(self->y, prec);
  mpfr_init2(self->z, prec);
  
  mpfr_init2(self->sigma, prec);
  mpfr_set(self->sigma, sigma, MPFR_RNDN);
  self->tailcut = tailcut;
  
  switch(algorithm) {
  case DGS_DISC_GAUSS_UNIFORM_TABLE: {
    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound, self->two_upper_bound_plus_one,
                                        self->sigma, self->tailcut);

    self->call = dgs_disc_gauss_mp_call_uniform_table;
    self->B = dgs_bern_uniform_init(0);
    _dgs_disc_gauss_mp_init_f(self->f, sigma);
    
    self->rho = (mpfr_t*)malloc(sizeof(mpfr_t)*mpz_get_ui(self->upper_bound));
    if (!self->rho) abort();
    mpfr_t x_;
    mpfr_init2(x_, prec);
    for(unsigned long x=0; x<mpz_get_ui(self->upper_bound); x++) {
      mpfr_set_ui(x_, x, MPFR_RNDN);
      mpfr_sqr(x_, x_, MPFR_RNDN);
      mpfr_mul(x_, x_, self->f, MPFR_RNDN);
      mpfr_exp(x_, x_, MPFR_RNDN);
      mpfr_init2(self->rho[x], prec);
      mpfr_set(self->rho[x], x_, MPFR_RNDN);
    }
    mpfr_div_ui(self->rho[0],self->rho[0],2, MPFR_RNDN);
    mpfr_clear(x_);
    break;
  }

  case DGS_DISC_GAUSS_UNIFORM_ONLINE: {
    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound, self->two_upper_bound_plus_one,
                                        self->sigma, self->tailcut);

    self->call = dgs_disc_gauss_mp_call_uniform_online;
    _dgs_disc_gauss_mp_init_f(self->f, self->sigma);

   break;
  }

  case DGS_DISC_GAUSS_UNIFORM_LOGTABLE: {
    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound, self->two_upper_bound_plus_one,
                                        self->sigma, self->tailcut);
    self->call = dgs_disc_gauss_mp_call_uniform_logtable;
    _dgs_disc_gauss_mp_init_bexp(self, self->sigma, self->upper_bound);
   break;
  }

  case DGS_DISC_GAUSS_SIGMA2_LOGTABLE: {
    self->call = dgs_disc_gauss_mp_call_sigma2_logtable;
    mpfr_t tmp;
    mpfr_init2(tmp, prec);

    mpfr_t sigma2;
    mpfr_init2(sigma2, prec);
    mpfr_set_ui(sigma2, 2, MPFR_RNDN);
    mpfr_log(sigma2, sigma2, MPFR_RNDN);
    mpfr_mul_ui(sigma2, sigma2, 2, MPFR_RNDN);
    mpfr_ui_div(sigma2, 1, sigma2, MPFR_RNDN);
    mpfr_sqrt(sigma2, sigma2, MPFR_RNDN);
    mpfr_div(tmp, sigma, sigma2, MPFR_RNDN);
    mpfr_get_z(self->k, tmp, MPFR_RNDN);
    mpfr_mul_z(self->sigma, sigma2, self->k, MPFR_RNDN);
    mpfr_clear(tmp);
    mpfr_clear(sigma2);

    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound, self->two_upper_bound_plus_one,
                                        self->sigma, self->tailcut);
    _dgs_disc_gauss_mp_init_bexp(self, self->sigma, self->upper_bound);
    self->B = dgs_bern_uniform_init(0);
    self->D2 = dgs_disc_gauss_sigma2p_init();
    break;
  }
    
  default:
    abort();
  }  
  return self;
}


void dgs_disc_gauss_mp_call_sigma2_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  do {
    do {
      dgs_disc_gauss_sigma2p_mp_call(self->x, self->D2, state);
      mpz_urandomm(self->y_z, state, self->k);
      mpz_mul(self->x2, self->k, self->x);
      mpz_mul_ui(self->x2, self->x2, 2);
      mpz_add(self->x2, self->x2, self->y_z);
      mpz_mul(self->x2, self->x2, self->y_z);
    } while (dgs_bern_exp_mp_call(self->Bexp, self->x2, state) == 0);
    mpz_mul(rop, self->k, self->x);
    mpz_add(rop, rop, self->y_z);
    if (mpz_sgn(rop) == 0) {
      if (dgs_bern_uniform_call(self->B, state))
        break;
    } else {
      break;
    }
  } while (1);
  if(dgs_bern_uniform_call(self->B, state))
    mpz_neg(rop, rop);
}

void dgs_disc_gauss_mp_call_uniform_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  do {
    mpz_urandomm(self->x, state, self->two_upper_bound_plus_one);
    mpz_sub(self->x, self->x, self->upper_bound);
    mpz_mul(self->x2, self->x, self->x);
  } while (dgs_bern_exp_mp_call(self->Bexp, self->x2, state) == 0);
  mpz_set(rop, self->x);
}

void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  unsigned long x;
  do {
    mpz_urandomm(self->x, state, self->upper_bound);
    x = mpz_get_ui(self->x);
    mpfr_urandomb(self->y, state);
  } while (mpfr_cmp(self->y, self->rho[x]) >= 0);

  mpz_set_ui(rop, x);
  if(dgs_bern_uniform_call(self->B, state))
    mpz_neg(rop, rop);
}

void dgs_disc_gauss_mp_call_uniform_online(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  do {
    mpz_urandomm(self->x, state, self->two_upper_bound_plus_one);
    mpz_sub(self->x, self->x, self->upper_bound);
    mpz_mul(self->x2, self->x, self->x);
    mpfr_mul_z(self->z, self->f, self->x2, MPFR_RNDN);
    mpfr_exp(self->z, self->z, MPFR_RNDN);
    mpfr_urandomb(self->y, state);
  } while (mpfr_cmp(self->y, self->z) >= 0);

  mpz_set(rop, self->x);
}


void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self) {
  mpfr_clear(self->sigma);
  if (self->B) dgs_bern_uniform_clear(self->B);
  if (self->Bexp) dgs_bern_exp_mp_clear(self->Bexp);
  mpz_clear(self->x);
  mpfr_clear(self->y);
  mpfr_clear(self->f);
  mpfr_clear(self->z);
  if (self->rho) {
    for(unsigned long x=0; x<mpz_get_ui(self->upper_bound); x++) {
      mpfr_clear(self->rho[x]);
    }
    free(self->rho);
  }
  free(self);
}



