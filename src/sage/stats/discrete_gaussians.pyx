"""
Discrete Gaussian Samplers

AUTHOR: Martin Albrecht <martinralbrecht+dgs@googlemail.com>
"""

from sage.rings.real_mpfr cimport RealNumber, RealField
from sage.libs.mpfr cimport mpfr_set, GMP_RNDN
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate

from dgs cimport dgs_disc_gauss_mp_init, dgs_disc_gauss_mp_clear
from dgs cimport DGS_DISC_GAUSS_UNIFORM_TABLE, DGS_DISC_GAUSS_UNIFORM_ONLINE, DGS_DISC_GAUSS_UNIFORM_LOGTABLE, DGS_DISC_GAUSS_SIGMA2_LOGTABLE

from dgs cimport dgs_disc_gauss_sigma2p_mp_init, dgs_disc_gauss_sigma2p_mp_call, dgs_disc_gauss_sigma2p_mp_clear

cdef class DiscreteGaussianSampler(SageObject):
    """A Discrete Gaussian Sampler using rejection sampling."""
    def __init__(self, sigma, tailcut=6, algorithm="uniform+table"):
        """
        
        INPUT:

        - ``sigma`` - samples are sampled with probability proportional to exp(x^2/(2sigma^2))
        - ``tailcut`` - (default: 6)
        - ``algorithm`` - supported choices are:
          - ``"uniform+table"`` -
          - ``"uniform+logtable"`` -
          - ``"uniform+online"`` - 
          - ``"sigma2+logtable"`` - 
        """
        if not isinstance(sigma, RealNumber):
            RR = RealField()
            sigma = RR(sigma)
            
        if algorithm == "uniform+table":
            algorithm = DGS_DISC_GAUSS_UNIFORM_TABLE
        elif algorithm == "uniform+online":
            algorithm = DGS_DISC_GAUSS_UNIFORM_ONLINE
        elif algorithm == "uniform+logtable":
            algorithm = DGS_DISC_GAUSS_UNIFORM_LOGTABLE
        elif algorithm == "sigma2+logtable":
            algorithm = DGS_DISC_GAUSS_SIGMA2_LOGTABLE
        else:
            raise ValueError("Algorithm '%s' not supported by class 'DiscreteGaussianSampler'"%(algorithm))
            
        self._gen = dgs_disc_gauss_mp_init((<RealNumber>sigma).value, tailcut, algorithm)
        self.sigma = sigma.parent()(0)
        mpfr_set(self.sigma.value, self._gen.sigma, GMP_RNDN)
        self.tailcut = Integer(tailcut)
        self.algorithm = algorithm

    def __clear__(self):
        if self._gen:
            dgs_disc_gauss_mp_clear(self._gen)

    def __call__(self):
        cdef randstate rstate = current_randstate()
        cdef Integer rop = Integer()
        self._gen.call(rop.value, self._gen, rstate.gmp_state)
        return rop

    def _repr_(self):
        return "Discrete Gaussian sampler with sigma = %f"%self.sigma

    def _sage_input_(self):
        return "DiscreteGaussianSampler(sigma=%f,tailcut=%d, algorithm='%s')"%(self.sigma, self.tailcut, self.algorithm)

cdef class DiscreteGaussianSamplerSigma2Plus(SageObject):
    def __init__(self):
        self._gen = dgs_disc_gauss_sigma2p_mp_init()

    def __clear__(self):
        if self._gen:
            dgs_disc_gauss_sigma2p_mp_clear(self._gen)

    def __call__(self):
        cdef randstate rstate = current_randstate()
        cdef Integer rop = Integer()
        dgs_disc_gauss_sigma2p_mp_call(rop.value, self._gen, rstate.gmp_state)
        return rop
