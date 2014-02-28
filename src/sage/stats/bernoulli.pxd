from dgs cimport dgs_bern_uniform_mp_t, dgs_bern_mp_t, dgs_bern_exp_mp_t

from sage.structure.sage_object cimport SageObject

cdef class BernoulliBase(SageObject):
    cdef tuple b
    
cdef class BernoulliUniformSampler(BernoulliBase):
    cdef dgs_bern_uniform_mp_t *_gen
    cpdef int raw_call(self)
    
cdef class BernoulliSampler(BernoulliBase):
    cdef dgs_bern_mp_t *_gen
    cpdef int raw_call(self)

cdef class BernoulliExpSampler(BernoulliBase):
    cdef dgs_bern_exp_mp_t *_gen
    cpdef int raw_call(self, x)

    