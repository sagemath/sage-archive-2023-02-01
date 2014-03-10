/**
 * \file dgs_gauss.h
 *
 * \author Martin Albrecht <martinralbrecht+dgs@googlemail.com>
 */

/******************************************************************************
*
*                        DGS - Discrete Gaussian Samplers
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

#ifndef DGS_GAUSS__H
#define DGS_GAUSS__H

#include "dgs_bern.h"

/**
 * Discrete Gaussians, shared definitions
 */

typedef enum {
  DGS_DISC_GAUSS_UNIFORM_ONLINE    = 0x1,
  DGS_DISC_GAUSS_UNIFORM_TABLE     = 0x2,
  DGS_DISC_GAUSS_UNIFORM_LOGTABLE  = 0x3,
  DGS_DISC_GAUSS_SIGMA2_LOGTABLE   = 0x7,
} dgs_disc_gauss_alg_t;

typedef struct {
  dgs_bern_uniform_t *B;
} dgs_disc_gauss_sigma2p_t;

dgs_disc_gauss_sigma2p_t *dgs_disc_gauss_sigma2p_init();
void dgs_disc_gauss_sigma2p_mp_call(mpz_t rop, dgs_disc_gauss_sigma2p_t *self, gmp_randstate_t state);
long dgs_disc_gauss_sigma2p_dp_call(dgs_disc_gauss_sigma2p_t *self);
void dgs_disc_gauss_sigma2p_clear(dgs_disc_gauss_sigma2p_t *self);

/**
 * Discrete Gaussians, multi-precision version
 */

struct _dgs_disc_gauss_mp_t;

typedef struct _dgs_disc_gauss_mp_t {
  mpfr_t sigma;
  mpz_t c;
  size_t tau;
  dgs_disc_gauss_alg_t algorithm;

  dgs_bern_uniform_t *B;
  dgs_bern_exp_mp_t *Bexp;
  dgs_disc_gauss_sigma2p_t *D2;
  
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

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, mpz_t c, size_t tau, dgs_disc_gauss_alg_t algorithm);
void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_call_uniform_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_call_uniform_online(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_call_sigma2_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);
void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self);

/**
 * \brief Discrete Gaussians, double-precision version
 *
 * Return integer x with probability
 *
 *    $ρ_{σ,c}(x) = exp(-(x-c)^2/(2σ^2))/exp(-(ZZ-c)^2/(2σ^2))$
 *
 * where $exp(-(ZZ-c)^2/(2σ^2)) ≈ \sum_{i=-τσ}^{τσ} exp(-(i-c)^2/(2σ^2))$ is the
 * probability for all of the integers.
 */

typedef struct _dgs_disc_gauss_dp_t {
  double sigma; //< width parameter σ
  long c; //< center parameter c
  size_t tau; //< samples outside of (c-τσ,...,c+τσ) are considered to have probability 0
  dgs_disc_gauss_alg_t algorithm; //< algorithm to use

  dgs_bern_uniform_t *B; //< source of uniform bits
  dgs_bern_exp_dp_t *Bexp; //< source of exponentially biased bits
  dgs_disc_gauss_sigma2p_t *D2; //< source of integers from disc. gaussian with σ = sqrt(1/(2log(2))) ≈ 0.849
  
  long (*call)(struct _dgs_disc_gauss_dp_t *self); //< call this function to sample

  double f; //< typically -1/(2σ^2)
  long upper_bound; //< τσ
  long two_upper_bound_plus_one;  //< 2τσ+1
  long k;  //< used in dgs_disc_gauss_dp_call_sigma2_logtable
  double *rho; //< used in dgs_disc_gauss_dp_call_uniform_table
} dgs_disc_gauss_dp_t;

dgs_disc_gauss_dp_t *dgs_disc_gauss_dp_init(double sigma, long c, size_t tau, dgs_disc_gauss_alg_t algorithm);

/**
 * \brief Sample from dgs_disc_gauss_dp_t by rejection sampling using the uniform
 * distribution
 *
 * \param self Discrete Gaussian sampler
 */

long dgs_disc_gauss_dp_call_uniform_online(dgs_disc_gauss_dp_t *self);

/**
 * \brief Sample from dgs_disc_gauss_dp_t by rejection sampling using the uniform
 * distribution and tabulated exp() evaluations.
 *
 * \param self Discrete Gaussian sampler
 */

long dgs_disc_gauss_dp_call_uniform_table(dgs_disc_gauss_dp_t *self);

/**
 * \brief Sample from dgs_disc_gauss_dp_t by rejection sampling using the
 * uniform distribution avoiding all exp() calls
 *
 * \param self Discrete Gaussian sampler
 */

long dgs_disc_gauss_dp_call_uniform_logtable(dgs_disc_gauss_dp_t *self);

/**
 * \brief Sample from dgs_disc_gauss_dp_t by rejection sampling using the
 * D_{k,σ2} distribution avoiding all exp() calls
 *
 * \param self Discrete Gaussian sampler
 */


long dgs_disc_gauss_dp_call_sigma2_logtable(dgs_disc_gauss_dp_t *self);

void dgs_disc_gauss_dp_clear(dgs_disc_gauss_dp_t *self);


#endif //DGS_GAUSS__H
