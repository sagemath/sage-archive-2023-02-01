"""
Schemes

AUTHORS:

- William Stein, David Kohel, Kiran Kedlaya (2008): added zeta_series
"""

#*******************************************************************************
#  Copyright (C) 2008 Kiran Kedlaya <kedlaya@mit.edu>
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
                            ZZ, is_RingHomomorphism, GF, PowerSeriesRing,
                            Rationals)

SCH = Schemes()

import homset

import morphism

import spec

def is_Scheme(x):
    """
    Return True if `x` is a scheme.

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import is_Scheme
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
    def __init__(self,X):
        """
        A scheme.

        TESTS::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: RmodI = R.quotient(I)
            sage: X = Spec(RmodI)
            sage: X == loads(dumps(X))
            True
        """
        if spec.is_Spec(X):
            self._base_ring = X.coordinate_ring()
            ParentWithBase.__init__(self, self._base_ring)
        else:
            self._base_scheme = X
            ParentWithBase.__init__(self, self._base_scheme)

    def __cmp__(self, X):
        """
        EXAMPLES::

            sage: X = Spec(QQ);  Y = Spec(QQ)
            sage: X == Y
            True
            sage: X is Y
            False
        """
        if not isinstance(X, self.__class__):
            return -1
        return self._cmp_(X)

    def union(self, X):
        """
        Return the disjoint union of the schemes self and X.

        EXAMPLES::

            sage: S = Spec(QQ)
            sage: X = AffineSpace(1, QQ)
            sage: S.union(X)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    __add__ = union

    def _point_morphism_class(self):
        """
        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X._point_morphism_class()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def base_extend(self, Y):
        """
        Y is either a scheme in the same category as self or a ring.

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X.base_scheme()
            Spectrum of Integer Ring
            sage: X.base_extend(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __call__(self, *args):
        """
        If S is a ring or scheme, return the set `X(S)` of
        `S`-valued points on `X`. If `S` is a list
        or tuple or just the coordinates, return a point in `X(T)`,
        where `T` is the base scheme of self.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)

        We create some point sets::

            sage: A(QQ)
            Set of Rational Points of Affine Space of dimension 2 over Rational Field
            sage: A(RR)
            Set of Rational Points of Affine Space of dimension 2 over Real Field with 53 bits of precision

        Space of dimension 2 over Rational Field

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: A(NumberField(x^2+1, 'a'))
            Set of Rational Points of Affine Space of dimension 2 over Number Field in a with defining polynomial x^2 + 1
            sage: A(GF(7))
            Traceback (most recent call last):
            ...
            ValueError: No natural map from the base ring (=Rational Field) to S (=Finite Field of size 7)

        We create some points::

            sage: A(QQ)([1,0])
            (1, 0)

        We create the same point by giving the coordinates of the point
        directly.

        ::

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

    def point_homset(self, S = None):
        """
        Return the set of S-valued points of this scheme.

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: P.point_homset(ZZ)
            Set of Rational Points of Projective Space of dimension 3 over Integer Ring
            sage: P.point_homset(QQ)
            Set of Rational Points of Projective Space of dimension 3 over Rational Field
            sage: P.point_homset(GF(11))
            Set of Rational Points of Projective Space of dimension 3 over Finite Field of size 11
        """
        if S is None or S == self.base_ring():
            # optimize common case
            try:
                return self.__ring_point_homset
            except AttributeError:
                self.__ring_point_homset = self._homset_class(self,self.base_ring())
                return self.__ring_point_homset
        try:
            return self.__point_homset[S]
        except AttributeError:
            self.__point_homset = {}
        except KeyError:
            pass
        H = self._homset_class(self, S)
        self.__point_homset[S] = H
        return H

    def point(self, v, check=True):
        return self._point_class(self, v, check=check)

    def _point_class(self):
        """
        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X._point_class()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _homset_class(self):
        """
        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X._homset_class()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def __div__(self, Y):
        """
        Return the base extension of self to Y.

        EXAMPLES::

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
        """
        Return the base ring of the scheme self.

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: A.base_ring()
            Rational Field

        ::

            sage: X = Spec(QQ)
            sage: X.base_ring()
            Integer Ring
        """
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
        """
        Return the base scheme of the scheme self.

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: A.base_scheme()
            Spectrum of Rational Field

        ::

            sage: X = Spec(QQ)
            sage: X.base_scheme()
            Spectrum of Integer Ring
        """
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
        """
        Return the structure morphism from the scheme self to its base
        scheme.

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: A.base_morphism()
            Scheme morphism:
              From: Affine Space of dimension 4 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map

        ::

            sage: X = Spec(QQ)
            sage: X.base_morphism()
            Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map
        """
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

    structure_morphism = base_morphism

    def category(self):
        """
        Return the category to which this scheme belongs. This is the
        category of all schemes of the base scheme of self.

        EXAMPLES::

            sage: ProjectiveSpace(4, QQ).category()
            Category of schemes over Spectrum of Rational Field
        """
        return Schemes(self.base_scheme())

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme, if defined. Otherwise
        raise a ValueError.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: X = Spec(R.quotient(I))
            sage: X.coordinate_ring()
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2)
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            raise ValueError, "This scheme has no associated coordinated ring (defined)."

    def dimension(self):
        """
        Return the relative dimension of this scheme over its base.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: X = Spec(R.quotient(I))
            sage: X.dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError # override in derived class

    def identity_morphism(self):
        """
        Return the identity morphism of the scheme self.

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X.identity_morphism()
            Scheme endomorphism of Spectrum of Rational Field
              Defn: Identity map
        """
        return morphism.SchemeMorphism_id(self)

    def hom(self, x, Y=None):
        """
        Return the scheme morphism from self to Y defined by x. If x is a
        scheme, try to determine a natural map to x.

        If Y is not given, try to determine Y from context.

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: P.hom(Spec(ZZ))
            Scheme morphism:
              From: Projective Space of dimension 3 over Integer Ring
              To:   Spectrum of Integer Ring
              Defn: Structure map
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

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: S = Spec(ZZ)
            sage: S._Hom_(P)
            Set of points of Projective Space of dimension 3 over Integer Ring defined over Integer Ring
        """
        return homset.SchemeHomset(self, Y, cat=cat, check=check)

    point_set = point_homset

    def count_points(self, n):
        r"""
        Count points over
        `\mathbf{F}_q, \ldots, \mathbf{F}_{q^n}` on a scheme over
        a finite field `\mathbf{F}_q`.

        .. note::

           This is currently only implemented for schemes over prime
           order finite fields.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: C = HyperellipticCurve(x^3+x^2+1)
            sage: C.count_points(4)
            [6, 12, 18, 96]
            sage: C.base_extend(GF(9,'a')).count_points(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Point counting only implemented for schemes over prime fields
        """

        F = self.base_ring()
        if not F.is_finite():
            raise TypeError, "Point counting only defined for schemes over finite fields"
        q = F.cardinality()
        if not q.is_prime():
            raise NotImplementedError, "Point counting only implemented for schemes over prime fields"
        a = []
        for i in range(1, n+1):
            F1 = GF(q**i, name='z')
            S1 = self.base_extend(F1)
            a.append(len(S1.rational_points()))
        return(a)

    def zeta_series(self, n, t):
        """
        Compute a power series approximation to the zeta function of a
        scheme over a finite field.

        INPUT:


        -  ``n`` - the number of terms of the power series to
           compute

        -  ``t`` - the variable which the series should be
           returned


        OUTPUT: A power series approximating the zeta function of self

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: C = HyperellipticCurve(x^3+x^2+1)
            sage: R.<t> = PowerSeriesRing(Integers())
            sage: C.zeta_series(4,t)
            1 + 6*t + 24*t^2 + 78*t^3 + 240*t^4 + O(t^5)
            sage: (1+2*t+3*t^2)/(1-t)/(1-3*t) + O(t^5)
            1 + 6*t + 24*t^2 + 78*t^3 + 240*t^4 + O(t^5)

        Note that this function depends on count_points, which is only
        defined for prime order fields::

            sage: C.base_extend(GF(9,'a')).zeta_series(4,t)
            Traceback (most recent call last):
            ...
            NotImplementedError: Point counting only implemented for schemes over prime fields
        """

        F = self.base_ring()
        if not F.is_finite():
            raise TypeError, "Zeta functions only defined for schemes over finite fields"
        a = self.count_points(n)
        R = PowerSeriesRing(Rationals(), 'u')
        u = R.gen()
        temp = sum(a[i-1]*(u.O(n+1))**i/i for i in range(1,n+1))
        temp2 = temp.exp()
        return(temp2(t).O(n+1))

def is_AffineScheme(x):
    """
    Return True if `x` is an affine scheme.

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import is_AffineScheme
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

        EXAMPLES: We construct the inclusion from
        `\mathrm{Spec}(\mathbb{Q})` into `\mathrm{Spec}(\mathbb{Z})`
        induced by the inclusion from `\mathbb{Z}` into
        `\mathbb{Q}`.

        ::

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
