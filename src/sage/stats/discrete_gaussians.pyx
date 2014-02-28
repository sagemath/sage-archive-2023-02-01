from sage.rings.real_mpfr cimport RealNumber, RealField
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate

from dgs cimport dgs_disc_gauss_mp_init, dgs_disc_gauss_mp_clear, DGS_DISC_GAUSS_UNIFORM, DGS_DISC_GAUSS_TABLE

cdef class DiscreteGaussianSampler(SageObject):
    def __init__(self, sigma, tailcut=6):
        if not isinstance(sigma, RealNumber):
            RR = RealField()
            sigma = RR(sigma)

        self._gen = dgs_disc_gauss_mp_init((<RealNumber>sigma).value, tailcut, DGS_DISC_GAUSS_TABLE & DGS_DISC_GAUSS_UNIFORM)

    def __clear__(self):
        if self._gen:
            dgs_disc_gauss_mp_clear(self._gen)

    def __call__(self):
        cdef randstate rstate = current_randstate()
        cdef Integer rop = Integer()
        self._gen.call(rop.value, self._gen, rstate.gmp_state)
        return rop
