"""
AUTHOR: Martin Albrecht <martinralbrecht@googlemail.com>
"""

from sage.libs.gmp.mpz cimport mpz_t
from sage.libs.gmp.random cimport gmp_randstate_t
from sage.libs.mpfr cimport mpfr_t
from libc.stdint cimport uint64_t

cdef extern from "dgs.h":
    int DGS_BERN_UNIFORM_DEFAULT_LENGTH
    int DGS_BERN_UNIFORM_MAX_LENGTH

    ctypedef struct dgs_bern_uniform_mp_t:
        size_t length
        size_t count
        mpz_t tmp
        uint64_t pool

    dgs_bern_uniform_mp_t *dgs_bern_uniform_mp_init(int length)
    unsigned long dgs_bern_uniform_mp_call(dgs_bern_uniform_mp_t *self, gmp_randstate_t state)
    void dgs_bern_uniform_mp_clear(dgs_bern_uniform_mp_t *self)

    ctypedef struct dgs_bern_mp_t:
        mpfr_t c
        mpfr_t tmp

    dgs_bern_mp_t *dgs_bern_mp_init(mpfr_t c)
    unsigned long dgs_bern_mp_call(dgs_bern_mp_t *self, gmp_randstate_t state)
    void dgs_bern_mp_clear(dgs_bern_mp_t *self)

    ctypedef struct dgs_bern_exp_mp_t:
        size_t l
        mpfr_t *c
        dgs_bern_mp_t **B

    dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l)
    unsigned long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x, gmp_randstate_t state)
    void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self)

    ctypedef int dgs_disc_gauss_alg_t
    cdef int DGS_DISC_GAUSS_UNIFORM_TABLE
    cdef int DGS_DISC_GAUSS_UNIFORM_ONLINE
    cdef int DGS_DISC_GAUSS_UNIFORM_LOGTABLE
    cdef int DGS_DISC_GAUSS_SIGMA2_LOGTABLE

    ctypedef struct dgs_disc_gauss_mp_t:
        mpfr_t sigma
        size_t tailcut
        dgs_disc_gauss_alg_t algorithm
        dgs_bern_mp_t *B
        void call(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state)

        unsigned long upper_bound
        mpz_t x
        mpfr_t y
        mpfr_t *rho
  
    dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, size_t tailcut, dgs_disc_gauss_alg_t algorithm)
    void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state)
    void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self)

    ctypedef struct dgs_disc_gauss_sigma2p_mp_t:
        dgs_bern_uniform_mp_t *B

    dgs_disc_gauss_sigma2p_mp_t *dgs_disc_gauss_sigma2p_mp_init()
    void dgs_disc_gauss_sigma2p_mp_call(mpz_t rop, dgs_disc_gauss_sigma2p_mp_t *self, gmp_randstate_t state)
    void dgs_disc_gauss_sigma2p_mp_clear(dgs_disc_gauss_sigma2p_mp_t *self)
