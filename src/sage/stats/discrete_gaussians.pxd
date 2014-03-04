from dgs cimport dgs_disc_gauss_mp_t, dgs_disc_gauss_sigma2p_mp_t

from sage.structure.sage_object cimport SageObject
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.integer cimport Integer

cdef class DiscreteGaussianSampler(SageObject):
    cdef RealNumber sigma
    cdef Integer tailcut
    cdef object algorithm
    cdef dgs_disc_gauss_mp_t *_gen

cdef class DiscreteGaussianSamplerSigma2Plus(SageObject):
    cdef dgs_disc_gauss_sigma2p_mp_t *_gen
