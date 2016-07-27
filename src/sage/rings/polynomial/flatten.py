r"""
Class to flatten polynomial rings over polynomial ring

For example ``QQ['a','b'],['x','y']`` flattens to ``QQ['a','b','x','y']``.


EXAMPLES::

    sage: R = QQ['x']['y']['s','t']['X']
    sage: from sage.rings.polynomial.flatten import FlatteningMorphism
    sage: phi = FlatteningMorphism(R); phi
    Generic morphism:
      From: Univariate Polynomial Ring in X over Multivariate Polynomial Ring in s, t over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
      To:   Multivariate Polynomial Ring in x, y, s, t, X over Rational Field
    sage: phi('x*y*s + t*X').parent()
    Multivariate Polynomial Ring in x, y, s, t, X over Rational Field


Authors:

Vincent Delecroix, Ben Hutz (July 2016): initial implementation
"""

#*****************************************************************************
#                  Copyright (C) 2016
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism
from copy import copy
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring_generic import is_MPolynomialRing

class FlatteningMorphism(Morphism):
    r"""
    EXAMPLES::

        sage: R = QQ['a','b']['x','y','z']['t1','t2']
        sage: from sage.rings.polynomial.flatten import FlatteningMorphism
        sage: f = FlatteningMorphism(R)
        sage: f.codomain()
        Multivariate Polynomial Ring in a, b, x, y, z, t1, t2 over Rational Field
        sage: p = R('(a+b)*x + (a^2-b)*t2*(z+y)')
        sage: p
        ((a^2 - b)*y + (a^2 - b)*z)*t2 + (a + b)*x
        sage: f(p)
        a^2*y*t2 + a^2*z*t2 - b*y*t2 - b*z*t2 + a*x + b*x
        sage: f(p).parent()
        Multivariate Polynomial Ring in a, b, x, y, z, t1, t2 over Rational Field

    Also works when univariate polynomial ring are involved::

        sage: R = QQ['x']['y']['s','t']['X']
        sage: from sage.rings.polynomial.flatten import FlatteningMorphism
        sage: f = FlatteningMorphism(R)
        sage: f.codomain()
        Multivariate Polynomial Ring in x, y, s, t, X over Rational Field
        sage: p = R('((x^2 + 1) + (x+2)*y + x*y^3)*(s+t) + x*y*X')
        sage: p
        x*y*X + (x*y^3 + (x + 2)*y + x^2 + 1)*s + (x*y^3 + (x + 2)*y + x^2 + 1)*t
        sage: f(p)
        x*y^3*s + x*y^3*t + x^2*s + x*y*s + x^2*t + x*y*t + x*y*X + 2*y*s + 2*y*t + s + t
        sage: f(p).parent()
        Multivariate Polynomial Ring in x, y, s, t, X over Rational Field    
    """
    def __init__(self, domain):
        """
        The Python constructor

        EXAMPLES::

            sage: R = ZZ['a', 'b', 'c']['x', 'y', 'z']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R)
            Generic morphism:
              From: Multivariate Polynomial Ring in x, y, z over Multivariate
            Polynomial Ring in a, b, c over Integer Ring
              To:   Multivariate Polynomial Ring in a, b, c, x, y, z over Integer
            Ring

        ::

            sage: R = ZZ['a']['b']['c']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R)
            Generic morphism:
              From: Univariate Polynomial Ring in c over Univariate Polynomial Ring
            in b over Univariate Polynomial Ring in a over Integer Ring
              To:   Multivariate Polynomial Ring in a, b, c over Integer Ring

        ::

            sage: R = ZZ['a']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R)
            Traceback (most recent call last):
            ...
            ValueError: clash in variable names

        ::

            sage: K.<v> = NumberField(x^3 - 2)
            sage: R = K['x','y']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f(R('v*a*x^2 + b^2 + 1/v*y'))
            (v)*x^2*a + b^2 + (1/2*v^2)*y

        ::

            sage: R = QQbar['x','y']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: f(R('QQbar(sqrt(2))*a*x^2 + b^2 + QQbar(I)*y'))
            1.414213562373095?*x^2*a + b^2 + I*y
        """
        if not is_PolynomialRing(domain) and not is_MPolynomialRing(domain):
            raise ValueError("domain should be a polynomial ring")

        ring = domain
        variables = []
        intermediate_rings = []

        while is_PolynomialRing(ring) or is_MPolynomialRing(ring):
            intermediate_rings.append(ring)
            v = ring.variable_names()
            if any(vv in variables for vv in v):
                raise ValueError("clash in variable names")
            variables.extend(reversed(v))
            ring = ring.base_ring()
        self._intermediate_rings = intermediate_rings
        variables.reverse()
        codomain = PolynomialRing(ring, variables)

        Morphism.__init__(self, domain, codomain)

    def _call_(self, p):
        r"""
        Evaluate an flatenning morphism.

        EXAMPLES::

            sage: R = QQ['a','b','c']['x','y','z']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: h = FlatteningMorphism(R)('2*a*x + b*z'); h
            2*a*x + b*z
            sage: h.parent()
            Multivariate Polynomial Ring in a, b, c, x, y, z over Rational Field

         TESTS::

             sage: R = QQ['x']['y']['s','t']
             sage: p = R('s*x + y*t + x^2*s + 1 + t')
             sage: from sage.rings.polynomial.flatten import FlatteningMorphism
             sage: f = FlatteningMorphism(R)
             sage: f._call_(p)
             x^2*s + x*s + y*t + t + 1
         """
        p = {(): p}

        for ring in self._intermediate_rings:
            new_p = {}
            if is_PolynomialRing(ring):
                for mon,pp in p.iteritems():
                    assert pp.parent() == ring
                    for i,j in pp.dict().iteritems():
                        new_p[(i,)+(mon)] = j
            elif is_MPolynomialRing(ring):
                for mon,pp in p.iteritems():
                    assert pp.parent() == ring
                    for mmon,q in pp.dict().iteritems():
                        new_p[tuple(mmon)+mon] = q
            else:
                raise RuntimeError
            p = new_p

        return self.codomain()(p)

    def section(self):
        """
        Inverse of this flattenning morphism.

        EXAMPLES::

            sage: R = QQ['a','b','c']['x','y','z']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: h = FlatteningMorphism(R)
            sage: h.section()
            Generic morphism:
              From: Multivariate Polynomial Ring in a, b, c, x, y, z over Rational
            Field
              To:   Multivariate Polynomial Ring in x, y, z over Multivariate
            Polynomial Ring in a, b, c over Rational Field

        ::

            sage: R = ZZ['a']['b']['c']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: FlatteningMorphism(R).section()
            Generic morphism:
              From: Multivariate Polynomial Ring in a, b, c over Integer Ring
              To:   Univariate Polynomial Ring in c over Univariate Polynomial Ring
            in b over Univariate Polynomial Ring in a over Integer Ring
        """
        phi= UnflatteningMorphism(self.codomain(), self.domain())
        phi._intermediate_rings = copy(self._intermediate_rings)
        phi._intermediate_rings.reverse()
        return phi

class UnflatteningMorphism(Morphism):
    r"""
    Inverses for :class:`FlatteningMorphism`

    EXAMPLES::

        sage: R = QQ['x']['y']['s','t']['X']
        sage: p = R.random_element()
        sage: from sage.rings.polynomial.flatten import FlatteningMorphism
        sage: f = FlatteningMorphism(R)
        sage: g = f.section()
        sage: g(f(p)) == p
        True
    """

    def _call_(self, p):
        """
        Evaluate an unflattening morphism.

        This is slow, but works.

        EXAMPLES::

            sage: R = QQ['x']['y']['a,b,c']
            sage: p = R.random_element()
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: g = f.section()
            sage: g(f(p)).parent()
            Multivariate Polynomial Ring in a, b, c over Univariate Polynomial Ring
            in y over Univariate Polynomial Ring in x over Rational Field

        ::

            sage: R = QQbar['x','y']['a','b']
            sage: from sage.rings.polynomial.flatten import FlatteningMorphism
            sage: f = FlatteningMorphism(R)
            sage: p = R('QQbar(sqrt(2))*a*x^2 + b^2 + QQbar(I)*y')
            sage: f.section()(f(p)) ==p
            True
        """
        p = p.dict()

        vars = [R.gens() for R in self._intermediate_rings]
        num = [len(v) for v in vars]
        f = 0
        for mon,pp in p.iteritems():
            ind = 0
            g = pp
            for i in range(len(num)):
                m = mon[ind:ind+num[i]]
                ind += num[i]
                if is_MPolynomialRing(self._intermediate_rings[i]):
                    g = g * self._intermediate_rings[i](dict({tuple(m):1}))
                else:
                    g = g * self._intermediate_rings[i](dict({m[0]:1}))
            f += g

        return f
