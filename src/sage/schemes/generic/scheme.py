"""
Schemes

AUTHORS:
   -- William Stein
   -- David Kohel
"""

#*******************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#  Copyright (C) 2005 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.structure.parent_base import ParentWithBase
from sage.categories.all import Schemes
from sage.rings.all import (IntegerRing, is_CommutativeRing, is_Field,
                            ZZ, is_RingHomomorphism)

SCH = Schemes()

import homset

import morphism

import spec

def is_Scheme(x):
    """
    Return True if $x$ is a scheme.

    EXAMPLES:
        sage: is_Scheme(5)
        False
        sage: X = Spec(QQ)
        sage: is_Scheme(X)
        True
    """
    return isinstance(x, Scheme)

# If a derived class sets any of the properties, _base_scheme, _base_ring,
# or _base_morphism, it will determine the others.  If none are set,
# the base defaults to Spec(Z) with the canonical morphism.

class Scheme(ParentWithBase):
    """
    A scheme.
    """
    def __init__(self,X):
        if spec.is_Spec(X):
            self._base_ring = X.coordinate_ring()
            ParentWithBase.__init__(self, self._base_ring)
        else:
            self._base_scheme = X
            ParentWithBase.__init__(self, self._base_scheme)

    def __cmp__(self, X):
        if not isinstance(X, self.__class__):
            return -1
        return self._cmp_(X)

    def __add__(self, other):
        return self.union(other)

    def _point_morphism_class(self):
        raise NotImplementedError

    def base_extend(self, Y):
        """
        Y is either a scheme in the same category as self or a ring.
        """
        raise NotImplementedError

    def __call__(self, *args):
        """
        If S is a ring or scheme, return the set $X(S)$ of $S$-valued
        points on $X$.  If $S$ is a list or tuple or just the coordinates,
        return a point in $X(T)$, where $T$ is the base scheme of self.

        EXAMPLES:
            sage: A = AffineSpace(2, QQ)

        We create some point sets:
            sage: A(QQ)
            Set of Rational Points of Affine Space of dimension 2 over Rational Field
            sage: A(RR)
            Set of Rational Points of Affine Space of dimension 2 over Real Field with 53 bits of precision

            Space of dimension 2 over Rational Field
            sage: R.<x> = PolynomialRing(QQ)
            sage: A(NumberField(x^2+1, 'a'))
            Set of Rational Points of Affine Space of dimension 2 over Number Field in a with defining polynomial x^2 + 1
            sage: A(GF(7))
            Traceback (most recent call last):
            ...
            ValueError: No natural map from the base ring (=Rational Field) to S (=Finite Field of size 7)


        We create some points:
            sage: A(QQ)([1,0])
            (1, 0)

        We create the same point by giving the coordinates of the point directly.
            sage: A( 1,0 )
            (1, 0)
        """
        if len(args) == 0:
            raise TypeError, "0-dimensional affine space has no points"

        S = args[0]
        if is_CommutativeRing(S):
            return self.point_homset(S)
        elif is_Scheme(S):
            return S.Hom(self)
        elif isinstance(S, (list, tuple)):
            args = S
        elif isinstance(S, morphism.SchemeMorphism_coordinates):
            if S.codomain() == self:
                return S
        else:
            # TODO: fix circular import resulting from non-multiple inheritance
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            if isinstance(S, EllipticCurvePoint_field):
                if S.codomain() == self:
                    return S
                else:
                    return self.point(S)
        return self.point(args)

    def point_homset(self, R=None):
        if R is None or R == self.base_ring():
            # optimize common case
            try:
                return self.__ring_point_homset
            except AttributeError:
                self.__ring_point_homset = self._homset_class(self,self.base_ring())
                return self.__ring_point_homset
        try:
            return self.__point_homset[R]
        except AttributeError:
            self.__point_homset = {}
        except KeyError:
            pass
        H = self._homset_class(self,R)
        self.__point_homset[R] = H
        return H

    def point(self, v, check=True):
        return self._point_class(self, v, check=check)

    def _point_class(self):
        raise NotImplementedError

    def _homset_class(self):
        raise NotImplementedError


    def __div__(self, Y):
        """
        Return the base extension of self to Y.

        EXAMPLES:
            sage: A = AffineSpace(3, ZZ)
            sage: A
            Affine Space of dimension 3 over Integer Ring
            sage: A/QQ
            Affine Space of dimension 3 over Rational Field
            sage: A/GF(7)
            Affine Space of dimension 3 over Finite Field of size 7
        """
        return self.base_extend(Y)

    def base_ring(self):
        try:
            return self._base_ring
        except AttributeError:
            if hasattr(self, '_base_morphism'):
                self._base_ring = self._base_morphism.codomain().coordinate_ring()
            elif hasattr(self, '_base_scheme'):
                self._base_ring = self._base_scheme.coordinate_ring()
            else:
                self._base_ring = ZZ
            return self._base_ring

    def base_scheme(self):
        try:
            return self._base_scheme
        except AttributeError:
            if hasattr(self, '_base_morphism'):
                self._base_scheme = self._base_morphism.codomain()
            elif hasattr(self, '_base_ring'):
                import spec
                self._base_scheme = spec.Spec(self._base_ring)
            else:
                import spec
                self._base_scheme = spec.SpecZ
            return self._base_scheme

    def base_morphism(self):
        try:
            return self._base_morphism
        except AttributeError:
            if hasattr(self, '_base_scheme'):
                self._base_morphism = self.Hom(self._base_scheme, cat=SCH).natural_map()
            elif hasattr(self, '_base_ring'):
                self._base_morphism = self.Hom(spec.Spec(self._base_ring), cat=SCH).natural_map()
            else:
                self._base_morphism = self.Hom(spec.SpecZ, cat=SCH).natural_map()
            return self._base_morphism

    def structure_morphism(self):
        """
        Same as self.base_morphism().
        """
        return self.base_morphism()

    def category(self):
        return Schemes(self.base_scheme())

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme, if defined.  Otherwise raise
        a ValueError.
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            raise ValueError, "This scheme has no associated coordinated ring (defined)."

    def dimension(self):
        """
        Return the relative dimension of this scheme over its base.
        """
        raise NotImplementedError # override in derived class

    def identity_morphism(self):
        return morphism.SchemeMorphism_id(self)

    def hom(self, x, Y=None):
        """
        Return the scheme morphism from self to Y defined by x.  If x is a scheme,
        try to determine a natural map to x.

        If Y is not given, try to determine Y from context.

        EXAMPLES:
        """
        if Y is None:
            if is_Scheme(x):
                return self.Hom(x).natural_map()
            else:
                raise TypeError, "unable to determine codomain"
        return self.Hom(Y)(x)

    def _Hom_(self, Y, cat=None, check=True):
        """
        Return the set of scheme morphisms from self to Y.
        """
        return homset.SchemeHomset(self, Y, cat=cat, check=check)

    def point_set(self, S):
        """
        Return the set of S-valued points of this scheme.
        """
        return self.point_homset(S)

def is_AffineScheme(x):
    """
    Return True if $x$ is an affine scheme.

    EXAMPLES:
        sage: is_AffineScheme(5)
        False
        sage: E = Spec(QQ)
        sage: is_AffineScheme(E)
        True
    """
    return isinstance(x, AffineScheme)

class AffineScheme(Scheme):
    """
    An abstract affine scheme.
    """
    def hom(self, x, Y=None):
        r"""
        Return the scheme morphism from self to Y defined by x.

        If Y is not given, try to determine from context.

        EXAMPLES:
        We construct the inclusion from $\Spec(\Q)$ into $\Spec(\Z)$
        induced by the inclusion from $\Z$ into $\Q$.
            sage: X = Spec(QQ)
            sage: X.hom(ZZ.hom(QQ))
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field
        """
        if is_Scheme(x):
            return self.Hom(x).natural_map()
        if Y is None:
            if is_RingHomomorphism(x):
                import spec
                Y = spec.Spec(x.domain())
        return Scheme.hom(self, x, Y)


#import morphism
#import spec
