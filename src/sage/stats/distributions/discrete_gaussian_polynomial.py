# -*- coding: utf-8 -*-
"""
Discrete Gaussian Samplers for `\\ZZ[x]`.

This class realizes oracles which returns polynomials in `\\ZZ[x]`
where each coefficient is sampled independently with a probability
proportional to `\exp(-(x-c)^2/(2σ^2))`.

AUTHORS:

- Martin Albrecht
- Robert Fitzpatrick
- Daniel Cabracas
- Florian Göpfert
- Michael Schneider

EXAMPLE::

    sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianPolynomialSampler
    sage: sigma = 3.0; n=1000
    sage: l = [DiscreteGaussianPolynomialSampler(ZZ['x'], 64, sigma)() for _ in xrange(n)]
    sage: l = map(lambda f: vector(f).norm().n(), l)
    sage: mean(l), sqrt(64)*sigma
    (23.83..., 24.0...)

"""

from sage.rings.all import RealField, RR, ZZ
from discrete_gaussian_integer import DiscreteGaussianIntegerSampler
from sage.structure.sage_object import SageObject

class DiscreteGaussianPolynomialSampler(SageObject):
    """
    Discrete Gaussian sampler for polynomials.

    EXAMPLE::

        sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianPolynomialSampler
        sage: DiscreteGaussianPolynomialSampler(ZZ['x'], 8, 3.0)()
        3*x^7 + 3*x^6 - 3*x^5 - x^4 - 5*x^2 + 3
        sage: gs = DiscreteGaussianPolynomialSampler(ZZ['x'], 8, 3.0)
        sage: [gs() for _ in xrange(3)]
        [4*x^7 + 4*x^6 - 2*x^5 + x^4 + 4*x^2 - 2*x + 7, 5*x^7 - 4*x^6 + 3*x^4 - 4*x^3 + x^2 - 4, 2*x^7 + 2*x^6 + x^5 + 2*x^3 - 3*x^2 + x]

    .. automethod:: __init__
    .. automethod:: __call__
    """
    def __init__(self, P, n, sigma):
        """
        Construct a sampler for univariate polynomials of degree ``n-1``
        where coefficients are drawn independently with standard deviation
        ``sigma``.

        INPUT:

        - ``P`` - a univariate polynomial ring over the Integers
        - ``n`` - number of coefficients to be sampled
        - ``sigma`` - coefficients `x` are accepted with probability proportional to
          `\exp(-x^2/(2σ^2))`. If an object of type
          :class:`sage.stats.distributions.discrete_gaussian_integer.DiscreteGaussianIntegerSampler` is passed, then this sampler
          is used to sample coefficients.

        EXAMPLE::

            sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianPolynomialSampler
            sage: DiscreteGaussianPolynomialSampler(ZZ['x'], 8, 3.0)()
            3*x^7 + 3*x^6 - 3*x^5 - x^4 - 5*x^2 + 3
            sage: gs = DiscreteGaussianPolynomialSampler(ZZ['x'], 8, 3.0)
            sage: [gs() for _ in xrange(3)]
            [4*x^7 + 4*x^6 - 2*x^5 + x^4 + 4*x^2 - 2*x + 7, 5*x^7 - 4*x^6 + 3*x^4 - 4*x^3 + x^2 - 4, 2*x^7 + 2*x^6 + x^5 + 2*x^3 - 3*x^2 + x]
        """
        if isinstance(sigma, DiscreteGaussianIntegerSampler):
            self.D = sigma
        else:
            self.D = DiscreteGaussianIntegerSampler(RR(sigma))
        self.n = ZZ(n)
        self.P = P

    def __call__(self):
        """
        Return a new sample.

        EXAMPLE::

            sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianPolynomialSampler
            sage: sampler = DiscreteGaussianPolynomialSampler(ZZ['x'], 8, 12.0)
            sage: sampler()
            8*x^7 - 11*x^5 - 19*x^4 + 6*x^3 - 34*x^2 - 21*x + 9
        """
        coeffs = [self.D() for _ in range(self.n)]
        f = self.P(coeffs)
        return f

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianPolynomialSampler
            sage: DiscreteGaussianPolynomialSampler(ZZ['x'], 8, 3.0)
            Discrete Gaussian sampler for polynomials of degree < 8 with σ=3.000000 in each component
        """
        return "Discrete Gaussian sampler for polynomials of degree < %d with σ=%f in each component"%(self.n, self.D.sigma)

