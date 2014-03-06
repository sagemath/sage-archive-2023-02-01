# -*- coding: utf-8 -*-
"""Discrete Gaussian Samplers.

This class realizes oracles which returns integers proportionally to
$exp(-x/(2*sigma^2))$. All oracles are implemented using Rejection sampling. See
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
from dgs cimport DGS_DISC_GAUSS_UNIFORM_TABLE, DGS_DISC_GAUSS_UNIFORM_ONLINE, DGS_DISC_GAUSS_UNIFORM_LOGTABLE, DGS_DISC_GAUSS_SIGMA2_LOGTABLE

cdef class DiscreteGaussianSampler(SageObject):
    """A Discrete Gaussian Sampler using rejection sampling."""
    def __init__(self, sigma, tailcut=6, algorithm="uniform+table"):
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
        if not isinstance(sigma, RealNumber):
            RR = RealField()
            sigma = RR(sigma)

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
            
        self._gen = dgs_disc_gauss_mp_init((<RealNumber>sigma).value, tailcut, algorithm)
        self.sigma = sigma.parent()(0)
        mpfr_set(self.sigma.value, self._gen.sigma, GMP_RNDN)
        self.tailcut = Integer(tailcut)
        self.algorithm = algorithm

    def __clear__(self):
        if self._gen:
            dgs_disc_gauss_mp_clear(self._gen)

    def __call__(self):
        """
        Return a new sample.
        """
        cdef randstate rstate = current_randstate()
        cdef Integer rop = Integer()
        self._gen.call(rop.value, self._gen, rstate.gmp_state)
        return rop

    def _repr_(self):
        return "Discrete Gaussian sampler with sigma = %f"%self.sigma

    def _sage_input_(self):
        return "DiscreteGaussianSampler(sigma=%f,tailcut=%d, algorithm='%s')"%(self.sigma, self.tailcut, self.algorithm)

