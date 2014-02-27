from sage.libs.gmp.mpz cimport mpz_t
from sage.libs.gmp.random cimport gmp_randstate_t
from sage.libs.mpfr cimport mpfr_t
from libc.stdint cimport uint64_t
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate
from sage.structure.sage_object cimport SageObject

cdef extern from "dgs.h":
    int DGS_BERN_UNIFORM_DEFAULT_LENGTH
    int DGS_BERN_UNIFORM_MAX_LENGTH

    ctypedef struct dgs_bern_uniform_t:
        size_t length
        size_t count
        mpz_t tmp
        uint64_t pool

    dgs_bern_uniform_t *dgs_bern_uniform_init(int length)
    unsigned long dgs_bern_uniform_call(dgs_bern_uniform_t *self, gmp_randstate_t state)
    void dgs_bern_uniform_clear(dgs_bern_uniform_t *self)

    ctypedef struct dgs_bern_t:
        mpfr_t c
        mpfr_t tmp

    dgs_bern_t *dgs_bern_init(mpfr_t c)
    unsigned long dgs_bern_call(dgs_bern_t *self, gmp_randstate_t state)
    void dgs_bern_clear(dgs_bern_t *self)

    ctypedef struct dgs_bern_exp_t:
        size_t l
        mpfr_t *c
        dgs_bern_t **B

    dgs_bern_exp_t* dgs_bern_exp_init(mpfr_t f, size_t l)
    unsigned long dgs_bern_exp_call(dgs_bern_exp_t *self, mpz_t x, gmp_randstate_t state)
    void dgs_bern_exp_clear(dgs_bern_exp_t *self)


cdef class BernoulliBase(SageObject):
    cdef tuple b
    
cdef class BernoulliUniformSampler(BernoulliBase):
    cdef dgs_bern_uniform_t *_gen
    cpdef int raw_call(self)
    
cdef class BernoulliSampler(BernoulliBase):
    cdef dgs_bern_t *_gen
    cpdef int raw_call(self)

cdef class BernoulliExpSampler(BernoulliBase):
    cdef dgs_bern_exp_t *_gen
    cpdef int raw_call(self, x)

    