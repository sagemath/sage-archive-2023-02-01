# -*- coding: utf-8 -*-
"""Discrete Gaussian Samplers.

This class realizes oracles which returns integers proportionally to
$exp(-x/(2*sigma^2))$. All oracles are implemented using rejection sampling. See
DiscreteGaussianSampler.__init__ to see which strategies are available.

AUTHOR: Martin Albrecht <martinralbrecht+dgs@googlemail.com>

EXAMPLE:

We sample proportionally to exp(-x^2/(2*sigma^2))::

    sage: from sage.stats.discrete_gaussians import DiscreteGaussianSampler
    sage: sigma = 3.0; n=100000
    sage: D = DiscreteGaussianSampler(sigma=sigma)
    sage: l = [D() for _ in xrange(n)]
    sage: c = sum([exp(-x^2/(2*sigma^2)) for x in xrange(-6*sigma,sigma*6+1)]); c
    7.519...
    sage: x=0; l.count(x), ZZ(round(n*exp(-x^2/(2*sigma^2))/c))
    (13363, 13298)
    sage: x=4; l.count(x), ZZ(round(n*exp(-x^2/(2*sigma^2))/c))
    (5415, 5467)
    sage: x=-10; l.count(x), ZZ(round(n*exp(-x^2/(2*sigma^2))/c))
    (40, 51)

REFERENCES:

.. [DDLL13] Léo Ducas, Alain Durmus, Tancrède Lepoint and Vadim
   Lyubashevsky. *Lattice Signatures and Bimodal Gaussians*; in Advances in
   Cryptology – CRYPTO 2013; Lecture Notes in Computer Science Volume 8042,
   2013, pp 40-56 http://www.di.ens.fr/~lyubash/papers/bimodal.pdf

"""

from sage.rings.real_mpfr cimport RealNumber, RealField
from sage.libs.mpfr cimport mpfr_set, GMP_RNDN
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate

from dgs cimport dgs_disc_gauss_mp_init, dgs_disc_gauss_mp_clear
from dgs cimport dgs_disc_gauss_dp_init, dgs_disc_gauss_dp_clear
from dgs cimport DGS_DISC_GAUSS_UNIFORM_TABLE, DGS_DISC_GAUSS_UNIFORM_ONLINE, DGS_DISC_GAUSS_UNIFORM_LOGTABLE, DGS_DISC_GAUSS_SIGMA2_LOGTABLE

cdef class DiscreteGaussianSampler(SageObject):
    """A Discrete Gaussian Sampler using rejection sampling."""
    def __init__(self, sigma, tailcut=6, algorithm="uniform+table", precision="mp"):
        """Construct a new sampler for a discrete Gaussian distribution.

        ALGORITHMS:

        - "uniform+table" - classical rejection sampling, sampling from the
          uniform distribution and accepted with probability
          $exp(-x^2/(2*sigma^2))$ where $exp(-x^2/(2*sigma^2))$ is precomputed
          and stored in a table.

        - "uniform+logtable" - samples are drawn from a uniform distribution and
          accepted with probability $exp(-x^2/(2*sigma^2))$ where
          $exp(-x^2/(2*sigma^2))$ is computed using logarithmically many calls
          to Bernoulli distributions. See [DDLL23]_ for details.

        - "uniform+online" - samples are drawn from a uniform distribution and
          accepted with probability $exp(-x^2/(2*sigma^2))$ where
          $exp(-x^2/(2*sigma^2))$ is computed in each invocation. Typically this
          is very slow.  See [DDLL23]_ for details.

        - "sigma2+logtable" - samples are drawn from an easily samplable
          distribution k*sigma2 and accepted with probability
          $exp(-x^2/(2*sigma^2))$ where $exp(-x^2/(2*sigma^2))$ is computed
          using logarithmically many calls to Bernoulli distributions.  See
          [DDLL23]_ for details. Note that this sampler adjusts sigma to match
          sigma2*k for some integer k.
        
        INPUT:

        - ``sigma`` - samples are sampled with probability proportional to
          $exp(x^2/(2sigma^2))$

        - ``tailcut`` - samples outside the range
          (-sigma*tailcut,...,sigma*tailcut) are considered to have probability
          zero (default: 6)

        - ``algorithm`` - see list above.

        - ``precision`` - either "mp" for multi-precision where the actual
          precision is taken from sigma or "dp" for double precision. In the
          latter case results are not reproducible across plattforms. (default:
          "mp").

        EXAMPLES::

            sage: from sage.stats.discrete_gaussians import DiscreteGaussianSampler
            sage: DiscreteGaussianSampler(3.0, algorithm="uniform+online")
            Discrete Gaussian sampler with sigma = 3.000000
            sage: DiscreteGaussianSampler(3.0, algorithm="uniform+table")
            Discrete Gaussian sampler with sigma = 3.000000
            sage: DiscreteGaussianSampler(3.0, algorithm="uniform+logtable")
            Discrete Gaussian sampler with sigma = 3.000000

        Note that "sigma2+logtable" adjusts sigma::
        
            sage: DiscreteGaussianSampler(3.0, algorithm="sigma2+logtable")
            Discrete Gaussian sampler with sigma = 3.397287

        TESTS::

            sage: from sage.stats.discrete_gaussians import DiscreteGaussianSampler
            sage: DiscreteGaussianSampler(-3.0)
            Traceback (most recent call last):
            ...
            ValueError: sigma must be > 0.0 but got -3.000000

            sage: DiscreteGaussianSampler(3.0, tailcut=-1)
            Traceback (most recent call last):
            ...
            ValueError: tailcut must be >= 1 but got -1

            sage: DiscreteGaussianSampler(3.0, tailcut=2, algorithm="superfastalgorithmyouneverheardof")        
            Traceback (most recent call last):
            ...
            ValueError: Algorithm 'superfastalgorithmyouneverheardof' not supported by class 'DiscreteGaussianSampler'

        """

        if sigma <= 0.0:
            raise ValueError("sigma must be > 0.0 but got %f"%sigma)

        if tailcut < 1:
            raise ValueError("tailcut must be >= 1 but got %d"%tailcut)
            
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

        if precision == "mp":
            if not isinstance(sigma, RealNumber):
                RR = RealField()
                sigma = RR(sigma)            
            self._gen_mp = dgs_disc_gauss_mp_init((<RealNumber>sigma).value, tailcut, algorithm)
            self._gen_dp = NULL
            self.sigma = sigma.parent()(0)
            mpfr_set(self.sigma.value, self._gen_mp.sigma, GMP_RNDN)
        elif precision == "dp":
            RR = RealField()
            if not isinstance(sigma, RealNumber):
                sigma = RR(sigma)            
            self._gen_dp = dgs_disc_gauss_dp_init(sigma, tailcut, algorithm)
            self._gen_mp = NULL
            self.sigma = RR(sigma)
        else:
            raise ValueError("Parameter precision '%s' not supported."%precision)
            
        self.tailcut = Integer(tailcut)
        self.algorithm = algorithm

    def __clear__(self):
        """
        TESTS::

            sage: from sage.stats.discrete_gaussians import DiscreteGaussianSampler
            sage: D = DiscreteGaussianSampler(3.0, algorithm="uniform+online")
            sage: del D
        """
        if self._gen_mp:
            dgs_disc_gauss_mp_clear(self._gen_mp)
        if self._gen_dp:
            dgs_disc_gauss_dp_clear(self._gen_dp)

    def __call__(self):
        """
        Return a new sample.

        EXAMPLES::

            sage: from sage.stats.discrete_gaussians import DiscreteGaussianSampler
            sage: DiscreteGaussianSampler(3.0, algorithm="uniform+online")()
            -3
            sage: DiscreteGaussianSampler(3.0, algorithm="uniform+table")()
            3

        TESTS::

            sage: from sage.stats.discrete_gaussians import DiscreteGaussianSampler
            sage: DiscreteGaussianSampler(3.0, algorithm="uniform+logtable", precision="dp")() # random output
            13
        """
        cdef randstate rstate
        cdef Integer rop
        if self._gen_mp:
            rstate = current_randstate()
            rop = Integer()
            self._gen_mp.call(rop.value, self._gen_mp, rstate.gmp_state)
            return rop
        else:
            return Integer(self._gen_dp.call(self._gen_dp))
            

    def _repr_(self):
        """
        TESTS::
        
            sage: from sage.stats.discrete_gaussians import DiscreteGaussianSampler
            sage: repr(DiscreteGaussianSampler(3.0))
            'Discrete Gaussian sampler with sigma = 3.000000'
        """
        return "Discrete Gaussian sampler with sigma = %f"%self.sigma

