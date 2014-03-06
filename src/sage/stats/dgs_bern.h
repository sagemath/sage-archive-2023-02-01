/**
 * \file Bernoulli samplers
 *
 * \author Martin Albrecht <martinralbrecht+dgs@googlemail.com>
 */

/******************************************************************************
*
*                      DGS - Discrete Gaussian Samplers
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

#ifndef DGS_BERN__H
#define DGS_BERN__H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>

#include "dgs_misc.h"

/**
 * \brief Number of bits sampled at once in dgs_bern_uniform_t
 */

#define DGS_BERN_UNIFORM_DEFAULT_LENGTH (sizeof(unsigned long)*8)

/**
 * \brief Maximum number of bits sampled at once in dgs_bern_uniform_t
 */

#define DGS_BERN_UNIFORM_MAX_LENGTH (sizeof(unsigned long)*8)

/**
 * \brief Number of blocks allocated in one go in dgs_bern_exp_mp_t and dgs_bern_exp_dp_t
 */

#define DGS_BERN_EXP_ALLOC_BLOCK_SIZE 16

/**
 * \brief balanced Bernoulli distribution
 *
 * \note This implementation is faster than dgs_bern_mp_t
 *
 * \ingroup Bernoulli
 */

typedef struct {
  size_t   length; //< number of bits we sample in each go
  size_t   count;  //< number of bits we have consumed yet
  mpz_t    tmp;    //< we sample to this mpz_t
  uint64_t pool;   //< we store the pool of random bits here
} dgs_bern_uniform_t;

/**
 * \brief construct a new uniform bit sampler
 *
 * \param length number of bits to sample from GMP
 *
 * \note clear with dgs_bern_uniform_clear
 * 
 * \ingroup Bernoulli
 */ 

dgs_bern_uniform_t* dgs_bern_uniform_init(size_t length);

/**
 * Sample a new uniformly random bit
 *
 * \param self Bernoulli state
 * \param state GMP randstate used as randomness source
 *
 * \ingroup Bernoulli
 */

static inline unsigned long dgs_bern_uniform_call(dgs_bern_uniform_t *self, gmp_randstate_t state) {
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

/**
 * \brief Clear uniformly random bit sampler
 * 
 * \ingroup Bernoulli
 */

void dgs_bern_uniform_clear(dgs_bern_uniform_t *self);

/**
 * \brief Bernoulli distribution, multi-precision version
 */

typedef struct {
  mpfr_t p; //< probability of 1
  mpfr_t tmp; //< used internally
} dgs_bern_mp_t;

/**
 * \brief Create new Bernoulli sampler
 *
 * \param p sampler will return 1 with probability p and 0 with probability 1-p
 *
 * \note clear with dgs_bern_mp_clear
 * 
 * \ingroup Bernoulli
 */

dgs_bern_mp_t *dgs_bern_mp_init(mpfr_t p);

/**
 * \brief Return 1 with probability p
 *
 * \param self Bernoulli state
 * \param state GMP randstate used as randomness source
 *
 * \ingroup Bernoulli
 */

long dgs_bern_mp_call(dgs_bern_mp_t *self, gmp_randstate_t state);

/**
 * \brief Clear Bernoulli sampler
 *
 * \param self Bernoulli state
 * 
 * \ingroup Bernoulli
 */

void dgs_bern_mp_clear(dgs_bern_mp_t *self);

/**
 * \brief Bernoulli distribution, double-precision version
 */

typedef struct {
  double p; //< probability of 1
  mpfr_t tmp; //< used internally
} dgs_bern_dp_t;

/**
 * \brief
 *
 * \param
 *
 * \note
 * 
 * \ingroup
 */

dgs_bern_dp_t *dgs_bern_dp_init(double p);


/**
 * \brief
 *
 * \param
 *
 * \note
 * 
 * \ingroup
 */

long dgs_bern_dp_call(dgs_bern_dp_t *self, gmp_randstate_t state);

/**
 * \brief
 *
 * \param
 *
 * \note
 * 
 * \ingroup
 */

void dgs_bern_dp_clear(dgs_bern_dp_t *self);

/**
 * \brief Bernoulli with p = exp(-x/f) for positive integers x, multi-precision version
 */

/**
 * \brief Bernoulli with p = exp(-x/f) for positive integers x, multi-precision version
 */

typedef struct {
  size_t l; //< we support positive x up to 2^l
  mpfr_t *p; //< probabilities for sub Bernoulli samplers
  dgs_bern_mp_t **B; //< sub Bernoulli samplers
} dgs_bern_exp_mp_t;

/**
 * \brief Create new family of Bernoulli samples
 *
 * \param f samplers return 1 with probability exp(-x/f)
 * \param x we support x up to 2^l
 *
 * \note clear dgs_bern_exp_mp_clear
 * 
 * \ingroup Bernoulli
 */

dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l);


/**
 * \brief Return 1 with probability -exp(x/f).
 *
 * \param self Bernoulli state
 * \param x integer > 0
 * \param state GMP randstate used as randomness source
 *
 * \ingroup Bernoulli
 */

long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x, gmp_randstate_t state);


/**
 * \brief Clear Bernoulli sampler family
 *
 * \param self Bernoulli state
 * 
 * \ingroup Bernoulli
 */

void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self);

/**
 * \brief Bernoulli with p = exp(-x/f) for positive integers x, double-precision version
 */

typedef struct {
  size_t l; //< we support positive x up to 2^l
  double *p; //< probabilities for sub Bernoulli samplers
  dgs_bern_dp_t **B; //< sub Bernoulli samplers
} dgs_bern_exp_dp_t;

/**
 * \brief Create new family of Bernoulli samples
 *
 * \param f samplers return 1 with probability exp(-x/f)
 * \param x we support x up to 2^l
 *
 * \note clear dgs_bern_exp_mp_clear
 * 
 * \ingroup Bernoulli
 */

dgs_bern_exp_dp_t* dgs_bern_exp_dp_init(double f, size_t l);

/**
 * \brief Return 1 with probability -exp(x/f).
 *
 * \param self Bernoulli state
 * \param x integer > 0
 * \param state GMP randstate used as randomness source
 *
 * \ingroup Bernoulli
 */

long dgs_bern_exp_dp_call(dgs_bern_exp_dp_t *self, long x, gmp_randstate_t state);

/**
 * \brief Clear Bernoulli sampler family
 *
 * \param self Bernoulli state
 * 
 * \ingroup Bernoulli
 */

void dgs_bern_exp_dp_clear(dgs_bern_exp_dp_t *self);

#endif //DGS_BERN__H
