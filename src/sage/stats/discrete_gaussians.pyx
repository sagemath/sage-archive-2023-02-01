from sage.rings.real_mpfr cimport RealNumber, RealField
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate

from dgs cimport dgs_disc_gauss_mp_init, dgs_disc_gauss_mp_clear, DGS_DISC_GAUSS_UNIFORM_TABLE, DGS_DISC_GAUSS_UNIFORM_ONLINE

cdef class DiscreteGaussianSampler(SageObject):
    def __init__(self, sigma, tailcut=6, algorithm="uniform+table"):
        if not isinstance(sigma, RealNumber):
            RR = RealField()
            sigma = RR(sigma)

        if algorithm == "uniform+table":
            algorithm = DGS_DISC_GAUSS_UNIFORM_TABLE
        elif algorithm == "uniform+online":
            algorithm = DGS_DISC_GAUSS_UNIFORM_ONLINE
        else:
            raise ValueError("Algorithm '%s' not supported by class 'DiscreteGaussianSampler'"%(algorithm))
            
        self._gen = dgs_disc_gauss_mp_init((<RealNumber>sigma).value, tailcut, algorithm)

    def __clear__(self):
        if self._gen:
            dgs_disc_gauss_mp_clear(self._gen)

    def __call__(self):
        cdef randstate rstate = current_randstate()
        cdef Integer rop = Integer()
        self._gen.call(rop.value, self._gen, rstate.gmp_state)
        return rop
