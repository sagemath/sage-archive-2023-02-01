# -*- coding: utf-8 -*-
r"""
Valuations on polynomial rings based on `\phi`-adic expansions

This file implements a base class for discrete valuations on polynomial rings,
defined by a `\phi`-adic expansion.

AUTHORS:

- Julian Rüth (2013-04-15): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013-2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the mac_lane/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.getcwd())
    sys.path.append(os.path.dirname(os.getcwd()))

from valuation import DiscretePseudoValuation
from sage.misc.abstract_method import abstract_method

from sage.misc.cachefunc import cached_method

class DevelopingValuation(DiscretePseudoValuation):
    r"""
    Abstract base class for a discrete valuation of polynomials defined over
    the polynomial ring ``domain`` by the `\phi`-adic development.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, pAdicValuation(QQ, 7))

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def __init__(self, parent, phi):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 7))
            sage: isinstance(v, DevelopingValuation)
            True

        """
        DiscretePseudoValuation.__init__(self, parent)

        domain = parent.domain()
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if not is_PolynomialRing(domain) or not domain.ngens() == 1:
            raise TypeError("domain must be a univariate polynomial ring but %r is not"%(domain,))

        phi = domain.coerce(phi)
        if phi.is_constant() or not phi.is_monic():
            raise ValueError("phi must be a monic non-constant polynomial but %r is not"%(phi,))

        self._phi = phi

    def phi(self):
        r"""
        Return the polynomial `\phi`, the key polynomial of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.phi()
            (1 + O(2^5))*x

        """
        return self._phi

    def effective_degree(self, f):
        r"""
        Return the effective degree of ``f`` with respect to this valuation.

        The effective degree of `f` is the largest `i` such that the valuation
        of `f` and the valuation of `f_i\phi^i` in the development `f=\sum_j
        f_j\phi^j` coincide (see [ML1936'] p.497.)

        INPUT:

        - ``f`` -- a non-zero polynomial in the domain of this valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.effective_degree(x)
            1
            sage: v.effective_degree(2*x + 1)
            0

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            raise ValueError("the effective degree is only defined for non-zero polynomials")

        v = self(f)
        return [i for i,w in enumerate(self.valuations(f)) if w == v][-1]

    def coefficients(self, f):
        r"""
        Return the `\phi`-adic expansion of ``f``.

        INPUT:

        - ``f`` -- a monic polynomial in the domain of this valuation

        OUTPUT:

        An iterator `[f_0,f_1,\dots]` of polynomials in the domain of this
        valuation such that `f=\sum_i f_i\phi^i`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: list(v.coefficients(f)) # note that these constants are in the polynomial ring
            [(1 + 2 + O(2^5)), (2 + O(2^6)), (1 + O(2^5))]
            sage: v = v.augmentation( x^2 + x + 1, 1)
            sage: list(v.coefficients(f))
            [(1 + O(2^5))*x + (2 + O(2^5)), (1 + O(2^5))]

        """
        f = self.domain().coerce(f)

        if self.phi().degree() == 1:
            from itertools import imap
            for c in imap(f.parent(), f(self.phi().parent().gen() - self.phi()[0]).coefficients(sparse=False)): yield c
        else:
            while f.degree() >= 0:
                f,r = self._quo_rem(f)
                yield r

    def _quo_rem(self, f):
        r"""
        Return the quotient and remainder of ``f`` divided by the key
        polynomial :meth:`phi`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 2))
            sage: v._quo_rem(x^2 + 1)
            (x, 1)

        """
        qr = [ self._quo_rem_monomial(i) for i in range(f.degree()+1) ]
        q = [ f[i]*g for i,(g,_) in enumerate(qr) ]
        r = [ f[i]*h for i,(_,h) in enumerate(qr) ]
        return sum(q), sum(r)

    @cached_method
    def _quo_rem_monomial(self, degree):
        r"""
        Return the quotient and remainder of `x^\mathrm{degree}` divided by the
        key polynomial :meth:`phy`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 2))
            sage: v._quo_rem_monomial(10)
            (x^9, 0)

        """
        f = self.domain().one() << degree
        return f.quo_rem(self.phi())

    def newton_polygon(self, f):
        r"""
        Return the newton polygon the `\phi`-adic development of ``f``.

        INPUT::

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v.newton_polygon(f)
            Finite Newton polygon with 2 vertices: (0, 0), (2, 0)

            sage: v = v.augmentation( x^2 + x + 1, 1)
            sage: v.newton_polygon(f)
            Finite Newton polygon with 2 vertices: (0, 0), (1, 1)
            sage: v.newton_polygon( f * v.phi()^3 )
            Finite Newton polygon with 2 vertices: (3, 3), (4, 4)

        """
        f = self.domain().coerce(f)

        from sage.geometry.newton_polygon import NewtonPolygon
        return NewtonPolygon(enumerate(self.valuations(f)))

    def _call_(self, f):
        r"""
        Evaluate this valuation at ``f``.

        INPUT::

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v(f)
            0

            sage: v = v.augmentation( x^2 + x + 1, 1)
            sage: v(f)
            0
            sage: v(f * v.phi()^3 )
            3
            sage: v(S.zero())
            +Infinity

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            from sage.rings.all import infinity
            return infinity

        return min(self.valuations(f))

    @abstract_method
    def valuations(self, f):
        r"""
        Return the valuations of the `f_i\phi^i` in the expansion `f=\sum f_i\phi^i`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        A list, each entry a rational numbers or infinity, the valuations of
        `f_0, f_1\phi, \dots`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S, pAdicValuation(R))
            sage: f = x^2 + 2*x + 16
            sage: list(v.valuations(f))
            [4, 1, 0]

        """

    def _test_effective_degree(self, **options):
        r"""
        Test the correctness of :meth:`effective_degree`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v._test_effective_degree()
        
        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().base_ring().some_elements())
        for x in S:
            if x == 0:
                continue
            tester.assertEqual(self.effective_degree(x), 0)
