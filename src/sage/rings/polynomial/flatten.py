r"""
Class to flatten polynomial rings over polynomial ring

For example ``QQ['a','b'],['x','y']`` flattens to ``QQ['a','b','c','d']``.


EXAMPLES::

    sage: R = QQ['x']['y']['s','t']['X']
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
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring_generic import is_MPolynomialRing

class FlatteningMorphism(Morphism):
    r"""
    EXAMPLES::

        sage: R = QQ['a','b']['x','y','z']['t1','t2']
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
            sage: FlatteningMorphism(R)
            Generic morphism:
              From: Multivariate Polynomial Ring in x, y, z over Multivariate
            Polynomial Ring in a, b, c over Integer Ring
              To:   Multivariate Polynomial Ring in a, b, c, x, y, z over Integer
            Ring
        
        ::

            sage: R = ZZ['a']['b']['c']
            sage: FlatteningMorphism(R)
            Generic morphism:
              From: Univariate Polynomial Ring in c over Univariate Polynomial Ring
            in b over Univariate Polynomial Ring in a over Integer Ring
              To:   Multivariate Polynomial Ring in a, b, c over Integer Ring
        """
        if not is_PolynomialRing(domain) and not is_MPolynomialRing(domain):
            raise ValueError("domain should be a polynomial ring")

        ring = domain
        variables = []
        while is_PolynomialRing(ring) or is_MPolynomialRing(ring):
            v = ring.variable_names()
            if any(vv in variables for vv in v):
                raise ValueError("clash in variable names")
            variables.extend(reversed(v))
            ring = ring.base_ring()
        variables.reverse()
        codomain = PolynomialRing(ring.base_ring(), variables)

        Morphism.__init__(self, domain, codomain)

    def _call_(self, p):
        """
        Evaluate an flatenning morphism.
        
        
        This is slow, but works.
        
        EXAMPLES::
        
            sage: R = QQ['a','b','c']['x','y','z']
            sage: h = FlatteningMorphism(R)('2*a*x + b*z'); h
            2*a*x + b*z
            sage: h.parent()
            Multivariate Polynomial Ring in a, b, c, x, y, z over Rational Field
            
        """
        return self.codomain()(str(p))

    def section(self):
        """
        Inverse of this flattenning morphism.
        
        EXAMPLES::
        
            sage: R = QQ['a','b','c']['x','y','z']
            sage: h = FlatteningMorphism(R)
            sage: h.section()
            Generic morphism:
              From: Multivariate Polynomial Ring in a, b, c, x, y, z over Rational
            Field
              To:   Multivariate Polynomial Ring in x, y, z over Multivariate
            Polynomial Ring in a, b, c over Rational Field

        ::

            sage: R = ZZ['a']['b']['c']
            sage: FlatteningMorphism(R).section()
            Generic morphism:
              From: Multivariate Polynomial Ring in a, b, c over Integer Ring
              To:   Univariate Polynomial Ring in c over Univariate Polynomial Ring
            in b over Univariate Polynomial Ring in a over Integer Ring
        """
        return UnflatteningMorphism(self.codomain(), self.domain())

class UnflatteningMorphism(Morphism):
    r"""
    Inverses for :class:`FlatteningMorphism`
    
    EXAMPLES::

        sage: R = QQ['x']['y']['s','t']['X']
        sage: p = R.random_element()
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
            sage: f = FlatteningMorphism(R)
            sage: g = f.section()
            sage: g(f(p)).parent()
            Multivariate Polynomial Ring in a, b, c over Univariate Polynomial Ring
            in y over Univariate Polynomial Ring in x over Rational Field
        """
        return self.codomain()(str(p))
