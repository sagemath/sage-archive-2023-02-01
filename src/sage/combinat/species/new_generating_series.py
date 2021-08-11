r"""
Generating Series

This file makes a number of extensions to lazy power series by
endowing them with some semantic content for how they're to be
interpreted.

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevant section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatse10.html.
One notable difference is that we use power-sum symmetric functions
as the coefficients of our cycle index series.

REFERENCES:

.. [BLL] \F. Bergeron, G. Labelle, and P. Leroux.
   "Combinatorial species and tree-like structures".
   Encyclopedia of Mathematics and its Applications, vol. 67, Cambridge Univ. Press. 1998.
.. [BLL-Intro] Francois Bergeron, Gilbert Labelle, and Pierre Leroux.
   "Introduction to the Theory of Species of Structures", March 14, 2008.
"""

#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.lazy_laurent_series import LazyTaylorSeries
from sage.rings.lazy_laurent_series_ring import LazyTaylorSeriesRing
from sage.rings.all import Integer, RationalField
from sage.arith.all import moebius, gcd, lcm, divisors
from sage.combinat.partition import Partition, Partitions
from functools import partial
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.cachefunc import cached_function
from sage.functions.other import factorial

class OrdinaryGeneratingSeriesRing(LazyTaylorSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.new_generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: R = OrdinaryGeneratingSeriesRing(QQ); R
            Lazy Taylor Series Ring in z over Rational Field
            sage: [R(lambda n: 1)[i] for i in range(4)]
            [1, 1, 1, 1]
            sage: R(lambda n: 1).counts(4)
            [1, 1, 1, 1]
            sage: R == loads(dumps(R))
            True

        TESTS:

        We test to make sure that caching works.

        ::

            sage: R is OrdinaryGeneratingSeriesRing(QQ)
            True

        """
        LazyTaylorSeriesRing.Element = OrdinaryGeneratingSeries
        LazyTaylorSeriesRing.__init__(self, R, 'z')


class OrdinaryGeneratingSeries(LazyTaylorSeries):
    def count(self, n):
        """
        Return the number of structures on a set of size ``n``.

        EXAMPLES::

            sage: from sage.combinat.species.new_generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.count(10)
            10
        """
        return self.coefficient(n)

    def counts(self, n):
        """
        Return the number of structures on a set for size ``i`` for
        each ``i`` in ``range(n)``.

        EXAMPLES::

            sage: from sage.combinat.species.new_generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.counts(10)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return [self.count(i) for i in range(n)]


class ExponentialGeneratingSeriesRing(LazyTaylorSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.new_generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ); R
            Lazy Taylor Series Ring in z over Rational Field
            sage: R == loads(dumps(R))
            True
            sage: [R(lambda n: 1).coefficient(i) for i in range(4)]
            [1, 1, 1, 1]
            sage: R(lambda n: 1).counts(4)
            [1, 1, 2, 6]

        TESTS:

        We test to make sure that caching works.

        ::

            sage: R is ExponentialGeneratingSeriesRing(QQ)
            True

        """
        LazyTaylorSeriesRing.Element = ExponentialGeneratingSeries
        LazyTaylorSeriesRing.__init__(self, R, 'z')


class ExponentialGeneratingSeries(LazyTaylorSeries):
    def count(self, n):
        """
        Return the number of structures of size ``n``.

        EXAMPLES::

            sage: from sage.combinat.species.new_generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: f = R(lambda n: 1)
            sage: [f.count(i) for i in range(7)]
            [1, 1, 2, 6, 24, 120, 720]
        """
        return factorial(n) * self.coefficient(n)

    def counts(self, n):
        """
        Return the number of structures on a set for size ``i`` for
        each ``i`` in ``range(n)``.

        EXAMPLES::

            sage: from sage.combinat.species.new_generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.counts(5)
            [0, 1, 4, 18, 96]
        """
        return [self.count(i) for i in range(n)]

    def functorial_composition(self, y):
        r"""
        Return the exponential generating series which is the functorial
        composition of ``self`` with ``y``.

        If `f = \sum_{n=0}^{\infty} f_n \frac{x^n}{n!}` and
        `g = \sum_{n=0}^{\infty} g_n \frac{x^n}{n!}`, then
        functorial composition `f \Box g` is defined as

        .. MATH::

             f \Box g = \sum_{n=0}^{\infty} f_{g_n} \frac{x^n}{n!}

        REFERENCES:

        - Section 2.2 of [BLL]_.

        EXAMPLES::

            sage: G = species.SimpleGraphSpecies()
            sage: g = G.generating_series()
            sage: g.coefficients(10)
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]
        """
        return self._new(partial(self._functorial_compose_gen, y), lambda a,b: 0, self, y)

    def _functorial_compose_gen(self, y, ao):
        """
        Returns a generator for the coefficients of the functorial
        composition of self with y.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: E2 = E.restricted(min=2, max=3)
            sage: WP = species.SubsetSpecies()
            sage: P2 = E2*E
            sage: g1 = WP.generating_series()
            sage: g2 = P2.generating_series()
            sage: g = g1._functorial_compose_gen(g2, 0)
            sage: [next(g) for i in range(10)]
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]
        """
        n = 0
        while True:
            yield self.count(y.count(n)) / factorial(n)
            n += 1