from dgs cimport dgs_disc_gauss_mp_t

from sage.structure.sage_object cimport SageObject

cdef class DiscreteGaussianSampler(SageObject):
    cdef dgs_disc_gauss_mp_t *_gen
