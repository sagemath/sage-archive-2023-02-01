/**
 * Discrete Gaussian Samplers
 *
 *
 */

/******************************************************************************
*
*                        Discrete Gaussian Samplers
*
* Copyright (c) 2014, Martin Albrecht  <martinralbrecht+dgs@googlemail.com> 
* All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* 
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer. 
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* The views and conclusions contained in the software and documentation are
* those of the authors and should not be interpreted as representing official
* policies, either expressed or implied, of the FreeBSD Project.
******************************************************************************/


#ifndef DGS__H
#define DGS__H

#include <assert.h>
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

static inline unsigned long dgs_bern_uniform_mp_call(dgs_bern_uniform_mp_t *self, gmp_randstate_t state) {
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

void dgs_bern_uniform_mp_clear(dgs_bern_uniform_mp_t *self);

/**
 * \brief Bernoulli distribution, multi-precision version
 */

typedef struct {
  mpfr_t p; //< probability of 1
  mpfr_t tmp; //< used internally
} dgs_bern_mp_t;

dgs_bern_mp_t *dgs_bern_mp_init(mpfr_t p);
unsigned long dgs_bern_mp_call(dgs_bern_mp_t *self, gmp_randstate_t state);
void dgs_bern_mp_clear(dgs_bern_mp_t *self);

/*
 * Bernoulli with p = exp(-x/f) for positive integers x, multi-precision version
 */

#define DGS_BERN_EXP_ALLOC_BLOCK_SIZE 16

typedef struct {
  size_t l; //< we support positive x up to 2^l
  mpfr_t *p; //< probabilities for sub Bernoulli samplers
  dgs_bern_mp_t **B; //< sub Bernoulli samplers
} dgs_bern_exp_mp_t;

dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l);
unsigned long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x, gmp_randstate_t state);
void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self);

/**
 * Discrete Gaussians
 */

typedef enum {
  DGS_DISC_GAUSS_UNIFORM_ONLINE    = 0x1,
  DGS_DISC_GAUSS_UNIFORM_TABLE     = 0x2,
  DGS_DISC_GAUSS_UNIFORM_LOGTABLE  = 0x3,
  DGS_DISC_GAUSS_SIGMA2_LOGTABLE   = 0x7,
} dgs_disc_gauss_alg_t;

struct _dgs_disc_gauss_sigma2p_mp_t;

typedef struct _dgs_disc_gauss_sigma2p_mp_t{
  dgs_bern_uniform_mp_t *B;
} dgs_disc_gauss_sigma2p_mp_t;

dgs_disc_gauss_sigma2p_mp_t *dgs_disc_gauss_sigma2p_mp_init();
void dgs_disc_gauss_sigma2p_mp_call(mpz_t rop, dgs_disc_gauss_sigma2p_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_sigma2p_mp_clear(dgs_disc_gauss_sigma2p_mp_t *self);

struct _dgs_disc_gauss_mp_t;

typedef struct _dgs_disc_gauss_mp_t {
  mpfr_t sigma;
  size_t tailcut;
  dgs_disc_gauss_alg_t algorithm;

  dgs_bern_uniform_mp_t *B;
  dgs_bern_exp_mp_t *Bexp;
  dgs_disc_gauss_sigma2p_mp_t *D2;
  
  void (*call)(mpz_t rop, struct _dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

  mpz_t upper_bound;
  mpz_t two_upper_bound_plus_one;
  mpz_t x;
  mpz_t y_z;
  mpz_t k;
  mpz_t x2;
  mpfr_t y;
  mpfr_t z;
  mpfr_t f;
  mpfr_t *rho;
  
} dgs_disc_gauss_mp_t;

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, size_t tailcut, dgs_disc_gauss_alg_t algorithm);
void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_call_uniform_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_call_uniform_online(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_call_sigma2_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self);


#endif //DGS__H
