from dgs cimport dgs_disc_gauss_mp_t, dgs_disc_gauss_dp_t

from sage.structure.sage_object cimport SageObject
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.integer cimport Integer

cdef class DiscreteGaussianSampler(SageObject):
    cdef RealNumber _sigma 
    cdef RealNumber _c
    cdef Integer _tau
    cdef object _algorithm
    cdef dgs_disc_gauss_mp_t *_gen_mp
    cdef dgs_disc_gauss_dp_t *_gen_dp

