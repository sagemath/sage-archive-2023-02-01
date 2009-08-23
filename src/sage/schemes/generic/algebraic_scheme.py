"""
Algebraic schemes

An algebraic scheme must be defined by sets of equations in affine
or projective spaces, perhaps by means of gluing relations.
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

from sage.rings.all import (
    is_Ideal,
    is_MPolynomialRing,
    is_FiniteField,
    is_RationalField,
    ZZ)

from sage.structure.all import Sequence

import ambient_space

import affine_space

import projective_space

import morphism

import scheme

def is_AlgebraicScheme(x):
    """
    Return True if `x` is an algebraic scheme, i.e., a
    subscheme of an ambient space over a ring defined by polynomial
    equations.

    EXAMPLES: Affine space is itself not an algebraic scheme, though
    the closed subscheme defined by no equations is.

    ::

        sage: from sage.schemes.generic.algebraic_scheme import is_AlgebraicScheme
        sage: is_AlgebraicScheme(AffineSpace(10, QQ))
        False
        sage: V = AffineSpace(10, QQ).subscheme([]); V
        Closed subscheme of Affine Space of dimension 10 over
        Rational Field defined by:
          (no equations)
        sage: is_AlgebraicScheme(V)
        True

    We create a more complicated closed subscheme.

    ::

        sage: A, x = AffineSpace(10, QQ).objgens()
        sage: X = A.subscheme([sum(x)]); X
        Closed subscheme of Affine Space of dimension 10 over Rational Field defined by:
        x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
        sage: is_AlgebraicScheme(X)
        True

    ::

        sage: is_AlgebraicScheme(QQ)
        False
        sage: S = Spec(QQ)
        sage: is_AlgebraicScheme(S)
        False
    """
    return isinstance(x, AlgebraicScheme)

class AlgebraicScheme(scheme.Scheme):
    """
    An algebraic scheme presented as a subscheme in an ambient space.
    """
    def __init__(self, A):
        """
        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme
            sage: P = ProjectiveSpace(3, ZZ)
            sage: S = AlgebraicScheme(P); S
            Subscheme of Projective Space of dimension 3 over Integer Ring
        """
        if not ambient_space.is_AmbientSpace(A):
            raise TypeError, "A (=%s) must be an ambient space"
        self.__A = A
        self.__divisor_group = {}

    def coordinate_ring(self):
        """
        Return the coordinate ring of this algebraic scheme.  The
        result is cached.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x-y, x-z])
            sage: S.coordinate_ring()
            Quotient of Multivariate Polynomial Ring in x, y, z over Integer Ring by the ideal (x - y, x - z)
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            R = self.__A.coordinate_ring()
            I = self.defining_ideal()
            Q = R.quotient(I)
            self._coordinate_ring = Q
            return Q

    def ambient_space(self):
        """
        Return the ambient space of this algebraic scheme.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(5))
            sage: S = A.subscheme([])
            sage: S.ambient_space()
            Affine Space of dimension 2 over Finite Field of size 5

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x-y, x-z])
            sage: S.ambient_space() is P
            True
        """
        return self.__A

    def ngens(self):
        """
        Return the number of generators of the ambient space of this
        algebraic scheme.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(5))
            sage: S = A.subscheme([])
            sage: S.ngens()
            2

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x-y, x-z])
            sage: P.ngens()
            3
        """
        return self.__A.ngens()

    def _repr_(self):
        """
        Return a string representation of this algebraic scheme.

        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme
            sage: P = ProjectiveSpace(3, ZZ)
            sage: S = AlgebraicScheme(P); S
            Subscheme of Projective Space of dimension 3 over Integer Ring
            sage: S._repr_()
            'Subscheme of Projective Space of dimension 3 over Integer Ring'
        """
        return "Subscheme of %s"%self.__A

    def _homset_class(self, *args, **kwds):
        return self.__A._homset_class(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return self.__A._point_class(*args, **kwds)


class AlgebraicScheme_quasi(AlgebraicScheme):
    """
    The quasi-affine or quasi-projective scheme `X - Y`, where `X` and `Y`
    are both closed subschemes of a common ambient affine or projective
    space.
    """
    def __init__(self, X, Y):
        """
        EXAMPLES::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_quasi
            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: T.complement(S)
            Quasi-projective scheme X - Y, where:
              X: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              (no equations)
              Y: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              x - y

            sage: AlgebraicScheme_quasi(S, T)
            Quasi-projective scheme X - Y, where:
              X: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              (no equations)
              Y: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              x - y
        """
        self.__X = X
        self.__Y = Y
        if not isinstance(X, AlgebraicScheme_subscheme):
            raise TypeError, "X must be a closed subscheme of an ambient space."
        if not isinstance(Y, AlgebraicScheme_subscheme):
            raise TypeError, "Y must be a closed subscheme of an ambient space."
        if X.ambient_space() != Y.ambient_space():
            raise ValueError, "X and Y must be embedded in the same ambient space."
        A = X.ambient_space()
        self._base_ring = A.base_ring()
        AlgebraicScheme.__init__(self, A)

    def _repr_(self):
        """
        Return a string representation of this algebraic scheme.

        EXAMPLES::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_quasi
            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = AlgebraicScheme_quasi(S, T); U
            Quasi-projective scheme X - Y, where:
              X: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              (no equations)
              Y: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              x - y

            sage: U._repr_()
            'Quasi-projective scheme X - Y, where:\n  X: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:\n  (no equations)\n  Y: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:\n  x - y\n'
        """
        if affine_space.is_AffineSpace(self.ambient_space()):
            t = "affine"
        else:
            t = "projective"
        s =  "Quasi-%s scheme X - Y, where:\n"%t
        s += "  X: %s\n"%self.__X
        s += "  Y: %s\n"%self.__Y
        return s

    def X(self):
        """
        Return the scheme `X` such that self is represented as `X - Y`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U.X() is S
            True
        """
        return self.__X

    def Y(self):
        """
        Return the scheme `Y` such that self is represented as `X - Y`.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U.Y() is T
            True
        """
        return self.__Y

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a TypeError.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([])
            sage: T = P.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U._check_satisfies_equations([1, 2, 0])
            True
            sage: U._check_satisfies_equations([1, 1, 0])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 1, 0] do not define a point on Quasi-projective scheme X - Y, where:
              X: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              (no equations)
              Y: Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              x - y

            sage: U._check_satisfies_equations([1, 4])
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match number of variables in parent


            sage: A.<x, y> = AffineSpace(2, GF(7))
            sage: S = A.subscheme([x^2-y])
            sage: T = A.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U._check_satisfies_equations([2, 4])
            True
            sage: U._check_satisfies_equations([1, 1])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 1] do not define a point on Quasi-affine scheme X - Y, where:
              X: Closed subscheme of Affine Space of dimension 2 over Finite Field of size 7 defined by:
              x^2 - y
              Y: Closed subscheme of Affine Space of dimension 2 over Finite Field of size 7 defined by:
              x - y
            sage: U._check_satisfies_equations([1, 0])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 0] do not define a point on Quasi-affine scheme X - Y, where:
              X: Closed subscheme of Affine Space of dimension 2 over Finite Field of size 7 defined by:
              x^2 - y
              Y: Closed subscheme of Affine Space of dimension 2 over Finite Field of size 7 defined by:
              x - y
        """
        for f in self.__X.defining_polynomials():
            if f(v) != 0:
                raise TypeError, "Coordinates %s do not define a point on %s"%(v,self)
        for f in self.__Y.defining_polynomials():
            if f(v) == 0:
                raise TypeError, "Coordinates %s do not define a point on %s"%(v,self)
        return True

    def rational_points(self, F=None, bound=0):
        """
        Return the set of rational points on this algebraic scheme
        over the field `F`.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(7))
            sage: S = A.subscheme([x^2-y])
            sage: T = A.subscheme([x-y])
            sage: U = T.complement(S)
            sage: U.rational_points()
            [(2, 4), (3, 2), (4, 2), (5, 4), (6, 1)]
            sage: U.rational_points(GF(7^2, 'b'))
            [(2, 4), (3, 2), (4, 2), (5, 4), (6, 1), (b, b + 4), (b + 1, 3*b + 5), (b + 2, 5*b + 1),
            (b + 3, 6), (b + 4, 2*b + 6), (b + 5, 4*b + 1), (b + 6, 6*b + 5), (2*b, 4*b + 2),
            (2*b + 1, b + 3), (2*b + 2, 5*b + 6), (2*b + 3, 2*b + 4), (2*b + 4, 6*b + 4),
            (2*b + 5, 3*b + 6), (2*b + 6, 3), (3*b, 2*b + 1), (3*b + 1, b + 2), (3*b + 2, 5),
            (3*b + 3, 6*b + 3), (3*b + 4, 5*b + 3), (3*b + 5, 4*b + 5), (3*b + 6, 3*b + 2),
            (4*b, 2*b + 1), (4*b + 1, 3*b + 2), (4*b + 2, 4*b + 5), (4*b + 3, 5*b + 3),
            (4*b + 4, 6*b + 3), (4*b + 5, 5), (4*b + 6, b + 2), (5*b, 4*b + 2), (5*b + 1, 3),
            (5*b + 2, 3*b + 6), (5*b + 3, 6*b + 4), (5*b + 4, 2*b + 4), (5*b + 5, 5*b + 6),
            (5*b + 6, b + 3), (6*b, b + 4), (6*b + 1, 6*b + 5), (6*b + 2, 4*b + 1), (6*b + 3, 2*b + 6),
            (6*b + 4, 6), (6*b + 5, 5*b + 1), (6*b + 6, 3*b + 5)]
        """
        if F is None:
            F = self.base_ring()

        if bound == 0:
            if is_RationalField(F):
                raise TypeError, "A positive bound (= %s) must be specified."%bound
            if not is_FiniteField(F):
                raise TypeError, "Argument F (= %s) must be a finite field."%F
        pts = []
        for P in self.ambient_space().rational_points(F):
            try:
                if self._check_satisfies_equations(list(P)):
                    pts.append(P)
            except TypeError:
                pass
        pts.sort()
        return pts


class AlgebraicScheme_subscheme(AlgebraicScheme):
    """
    An algebraic scheme presented as a closed subscheme is defined by
    explicit polynomial equations. This is as opposed to a general
    scheme, which could, e.g., be the Neron model of some object, and
    for which we do not want to give explicit equations.

    INPUT:


    -  ``A`` - ambient space (affine or projective n-space
       over a ring)

    -  ``polys`` - ideal or tuple of defining polynomials
    """
    def __init__(self, A, polys):
        """
        TESTS::

            sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: P.subscheme([x^2-y*z])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
            sage: AlgebraicScheme_subscheme(P, [x^2-y*z])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
        """
        AlgebraicScheme.__init__(self, A)
        self._base_ring = A.base_ring()
        if is_Ideal(polys):
            self.__I = polys
            polys = polys.gens()
        if not isinstance(polys, (list, tuple)):
            polys = (polys, )
        else:
            polys = tuple(polys)
        self.__polys = A._validate(polys)  # now polys is a tuple of defining polynomials.

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a TypeError.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: S = P.subscheme([x^2-y*z])
            sage: S._check_satisfies_equations([1, 1, 1])
            True
            sage: S._check_satisfies_equations([1, 0, 1])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 0, 1] do not define a point on Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
            sage: S._check_satisfies_equations([0, 0, 0])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [0, 0, 0] do not define a point on Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - y*z
        """
        for f in self.defining_polynomials():
            if f(v) != 0:   # it must be "!=0" instead of "if f(v)", e.g.,
                            # because of p-adic base rings.
                raise TypeError, "Coordinates %s do not define a point on %s"%(list(v),self)
        try:
            return self.ambient_space()._check_satisfies_equations(v)
        except TypeError:
            raise TypeError, "Coordinates %s do not define a point on %s"%(list(v),self)

    def base_extend(self, R):
        """
        Return the base change to the ring `R` of this scheme.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(11))
            sage: S = P.subscheme([x^2-y*z])
            sage: S.base_extend(GF(11^2, 'b'))
            Closed subscheme of Projective Space of dimension 2 over Finite Field in b of size 11^2 defined by:
              x^2 - y*z
            sage: S.base_extend(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: No natural map from the base ring (=Finite Field of size 11) to S (=Integer Ring)
        """
        A = self.ambient_space().base_extend(R)
        return A.subscheme(self.__polys)

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, QQ)
            sage: X = A.subscheme([x*y, z])
            sage: X == A.subscheme([z, x*y])
            True
            sage: X == A.subscheme([x*y, z^2])
            False
            sage: B.<u, v, t> = AffineSpace(3, QQ)
            sage: X == B.subscheme([u*v, t])
            False
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            return -1
        A = self.ambient_space()
        if other.ambient_space() != A:
            return -1
        return cmp(self.defining_ideal(), other.defining_ideal())

    def _repr_(self):
        """
        Return a string representation of this scheme.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(11))
            sage: S = P.subscheme([x^2-y*z])
            sage: S
            Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:
              x^2 - y*z
            sage: S._repr_()
            'Closed subscheme of Projective Space of dimension 2 over Finite Field of size 11 defined by:\n  x^2 - y*z'
        """
        polys = '\n  '.join([str(f) for f in self.defining_polynomials()])
        if polys == '':
            polys = '(no equations)'
        return "Closed subscheme of %s defined by:\n  %s"%(self.ambient_space(),
                                                           polys)

    def defining_polynomials(self):
        """
        Return the polynomials that define this scheme as a subscheme
        of its ambient space.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x^2-y*z, x^3+z^3])
            sage: S.defining_polynomials()
            (x^2 - y*z, x^3 + z^3)
        """
        return self.__polys

    def defining_ideal(self):
        """
        Return the ideal that defines this scheme as a subscheme
        of its ambient space.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([x^2-y*z, x^3+z^3])
            sage: S.defining_ideal()
            Ideal (x^2 - y*z, x^3 + z^3) of Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        try:
            return self.__I
        except AttributeError:
            R = self.ambient_space().coordinate_ring()
            self.__I = R.ideal(self.defining_polynomials())
            return self.__I

    def irreducible_components(self):
        r"""
        Return the irreducible components of this algebraic scheme, as
        subschemes of the same ambient space.

        OUTPUT: an immutable sequence of irreducible subschemes of the
        ambient space of this scheme

        The components are cached.

        EXAMPLES:

        We define what is clearly a union of four hypersurfaces in
        `\P^4_{\QQ}` then find the irreducible components.

        ::

            sage: PP.<x,y,z,w,v> = ProjectiveSpace(4,QQ)
            sage: V = PP.subscheme( (x^2 - y^2 - z^2)*(w^5 -  2*v^2*z^3)* w * (v^3 - x^2*z) )
            sage: V.irreducible_components()
            [
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            w,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            x^2 - y^2 - z^2,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            x^2*z - v^3,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
            w^5 - 2*z^3*v^2
            ]
        """
        try:
            return self.__irreducible_components
        except AttributeError:
            pass
        I = self.defining_ideal()
        P = I.associated_primes()
        A = self.ambient_space()
        C = Sequence([A.subscheme(X) for X in P], check=False, cr=True)
        C.sort()
        C.set_immutable()
        self.__irreducible_components = C
        return C

    def reduce(self):
        r"""
        Return the corresponding reduced algebraic space associated to this
        scheme.

        EXAMPLES: First we construct the union of a doubled and tripled
        line in the affine plane over `\QQ`.

        ::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: X = A.subscheme([(x-1)^2*(x-y)^3]); X
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^5 - 3*x^4*y + 3*x^3*y^2 - x^2*y^3 - 2*x^4 + 6*x^3*y - 6*x^2*y^2 + 2*x*y^3 + x^3 - 3*x^2*y + 3*x*y^2 - y^3
            sage: X.dimension()
            1

        Then we compute the corresponding reduced scheme.

        ::

            sage: Y = X.reduce(); Y
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^2 - x*y - x + y

        Finally, we verify that the reduced scheme `Y` is the union
        of those two lines.

        ::

            sage: L1 = A.subscheme([x-1]); L1
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x - 1
            sage: L2 = A.subscheme([x-y]); L2
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x - y
            sage: W = L1.union(L2); W             # taken in ambient space
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x^2 - x*y - x + y
            sage: Y == W
            True
        """
        try:
            return self._reduce
        except AttributeError:
            r = self.defining_ideal().radical()
            A = self.ambient_space()
            V = A.subscheme(r)
            V._reduce = V       # so knows it is already reduced!
            self._reduce = V
            return V

    def union(self, other):
        """
        Return the scheme-theoretic union of self and other in their common
        ambient space.

        EXAMPLES: We construct the union of a line and a tripled-point on
        the line.

        ::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: I = ideal([x,y])^3
            sage: P = A.subscheme(I)
            sage: L = A.subscheme([y-1])
            sage: S = L.union(P); S
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
            y^4 - y^3
            x*y^3 - x*y^2
            x^2*y^2 - x^2*y
            x^3*y - x^3
            sage: S.dimension()
            1
            sage: S.reduce()
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
            y^2 - y
            x*y - x

        We can also use the notation "+" for the union::

            sage: A.subscheme([x]) + A.subscheme([y^2 - (x^3+1)])
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
            -x^4 + x*y^2 - x

        Saving and loading::

            sage: loads(S.dumps()) == S
            True
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError, "other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
        return A.subscheme(self.defining_ideal().intersection(other.defining_ideal()))

    __add__ = union

    def intersection(self, other):
        """
        Return the scheme-theoretic intersection of self and other in their
        common ambient space.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: X = A.subscheme([x^2-y])
            sage: Y = A.subscheme([y])
            sage: X.intersection(Y)
            Closed subscheme of Affine Space of dimension 2 over Integer Ring defined by:
              x^2 - y
              y
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError, "other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
        return A.subscheme(self.defining_ideal() + other.defining_ideal())


    def complement(self, other=None):
        """
        Return the scheme-theoretic complement other - self, where
        self and other are both closed algebraic subschemes of the
        same ambient space.

        If other is unspecified, it is taken to be the ambient space
        of self.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: X = A.subscheme([x+y-z])
            sage: Y = A.subscheme([x-y+z])
            sage: Y.complement(X)
            Quasi-affine scheme X - Y, where:
              X: Closed subscheme of Affine Space of dimension 3 over Integer Ring defined by:
              x + y - z
              Y: Closed subscheme of Affine Space of dimension 3 over Integer Ring defined by:
              x - y + z
            sage: Y.complement()
            Quasi-affine scheme X - Y, where:
              X: Closed subscheme of Affine Space of dimension 3 over Integer Ring defined by:
              (no equations)
              Y: Closed subscheme of Affine Space of dimension 3 over Integer Ring defined by:
              x - y + z

            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: X = P.subscheme([x^2+y^2+z^2])
            sage: Y = P.subscheme([x*y+y*z+z*x])
            sage: Y.complement(X)
            Quasi-projective scheme X - Y, where:
              X: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 + y^2 + z^2
              Y: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*y + x*z + y*z
            sage: Y.complement(P)
            Quasi-projective scheme X - Y, where:
              X: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              (no equations)
              Y: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*y + x*z + y*z
        """
        A = self.ambient_space()
        if other is None:
            other = A.subscheme([])
        elif not isinstance(other, AlgebraicScheme_subscheme):
            if other == A:
                other = A.subscheme([])
            else:
                raise TypeError, \
                      "Argument other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
        return AlgebraicScheme_quasi(other, self)

    def rational_points(self, F=None, bound=0):
        """
        EXAMPLES:

        One can enumerate points up to a given bound on a projective scheme
        over the rationals.

        ::

            sage: E = EllipticCurve('37a')
            sage: E.rational_points(bound=8)
            [(0 : 0 : 1),
             (1 : 0 : 1),
             (-1 : 0 : 1),
             (0 : -1 : 1),
             (1 : -1 : 1),
             (-1 : -1 : 1),
             (2 : 2 : 1),
             (2 : -3 : 1),
             (1/4 : -3/8 : 1),
             (1/4 : -5/8 : 1),
             (0 : 1 : 0)]

        For a small finite field, the complete set of points can be
        enumerated.

        ::

            sage: Etilde = E.base_extend(GF(3))
            sage: Etilde.rational_points()
            [(0 : 0 : 1), (0 : 1 : 0), (0 : 2 : 1), (1 : 0 : 1), (1 : 2 : 1), (2 : 0 : 1), (2 : 2 : 1)]

        The class of hyperelliptic curves does not (yet) support
        desingularization of the places at infinity into two points.

        ::

            sage: FF = FiniteField(7)
            sage: P.<x> = PolynomialRing(FiniteField(7))
            sage: C = HyperellipticCurve(x^8+x+1)
            sage: C.rational_points()
            [(2 : 0 : 1), (4 : 0 : 1), (0 : 1 : 1), (6 : 1 : 1), (0 : 6 : 1), (6 : 6 : 1), (0 : 1 : 0)]

        TODO:

        1. The above algorithms enumerate all projective points and
           test whether they lie on the scheme; Implement a more naive
           sieve at least for covers of the projective line.

        2. Implement Stoll's model in weighted projective space to
           resolve singularities and find two points (1 : 1 : 0) and
           (-1 : 1 : 0) at infinity.
        """
        if F == None:
            F = self.base_ring()
        X = self(F)
        if is_RationalField(F) or F == ZZ:
            if not bound > 0:
                raise TypeError, "A positive bound (= %s) must be specified."%bound
            try:
                return X.points(bound)
            except TypeError:
                raise TypeError, "Unable to enumerate points over %s."%F
        try:
            return X.points()
        except TypeError:
            raise TypeError, "Unable to enumerate points over %s."%F

class AlgebraicScheme_subscheme_affine(AlgebraicScheme_subscheme):
    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_affine_space(*args, **kwds)

    def dimension(self):
        """
        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: A.subscheme([]).dimension()
            2
            sage: A.subscheme([x]).dimension()
            1
            sage: A.subscheme([x^5]).dimension()
            1
            sage: A.subscheme([x^2 + y^2 - 1]).dimension()
            1
            sage: A.subscheme([x*(x-1), y*(y-1)]).dimension()
            0

        Something less obvious

        ::

            sage: A.<x,y,z,w> = AffineSpace(4, QQ)
            sage: X = A.subscheme([x^2, x^2*y^2 + z^2, z^2 - w^2, 10*x^2 + w^2 - z^2])
            sage: X
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
            x^2
            x^2*y^2 + z^2
            z^2 - w^2
            10*x^2 - z^2 + w^2
            sage: X.dimension()
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension()
            return self.__dimension

    def projective_embedding(self, i=None, X=None):
        """
        Returns a morphism from this affine scheme into an ambient
        projective space of the same dimension.

        INPUT:


        -  ``i`` - integer (default: dimension of self = last
           coordinate) determines which projective embedding to compute. The
           embedding is that which has a 1 in the i-th coordinate, numbered
           from 0.


        -  ``X`` - (default: None) projective scheme, i.e., codomain of
           morphism; this is constructed if it is not given.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: S = A.subscheme([x*y-z])
            sage: S.projective_embedding()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Integer Ring defined by:
              x*y - z
              To:   Closed subscheme of Projective Space of dimension 3 over Integer Ring defined by:
              x0*x1 - x2*x3
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x : y : z : 1)
        """
        AA = self.ambient_space()
        n = AA.dimension_relative()
        if i is None:
            try:
                i = self._default_embedding_index
            except AttributeError:
                i = int(n)
        else:
            i = int(i)
        if i < 0 or i > n:
            raise ValueError, \
                  "Argument i (=%s) must be between 0 and %s, inclusive"%(i, n)
        try:
            return self.__projective_embedding[i]
        except AttributeError:
            self.__projective_embedding = {}
        except KeyError:
            pass
        if X is None:
            PP = projective_space.ProjectiveSpace(n, AA.base_ring())
            v = list(PP.gens())
            z = v.pop(i)
            v.append(z)
            polys = self.defining_polynomials()
            X = PP.subscheme([ f.homogenize()(v) for f in polys ])
        R = AA.coordinate_ring()
        v = list(R.gens())
        v.insert(i, R(1))
        phi = self.hom(v, X)
        self.__projective_embedding[i] = phi
        return phi


class AlgebraicScheme_subscheme_projective(AlgebraicScheme_subscheme):
    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_projective_space(*args, **kwds)

    def dimension(self):
        """
        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: A.subscheme([]).dimension()
            2
            sage: A.subscheme([x]).dimension()
            1
            sage: A.subscheme([x^5]).dimension()
            1
            sage: A.subscheme([x^2 + y^2 - 1]).dimension()
            1
            sage: A.subscheme([x*(x-1), y*(y-1)]).dimension()
            0

        Something less obvious

        ::

            sage: A.<x,y,z,w> = AffineSpace(4, QQ)
            sage: X = A.subscheme([x^2, x^2*y^2 + z^2, z^2 - w^2, 10*x^2 + w^2 - z^2])
            sage: X
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
            x^2
            x^2*y^2 + z^2
            z^2 - w^2
            10*x^2 - z^2 + w^2
            sage: X.dimension()
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension() - 1
            return self.__dimension

    def affine_patch(self, i):
        r"""
        Return the `i^{th}` affine patch of this projective scheme.
        This is the intersection with this `i^{th}` affine patch of
        its ambient space.

        INPUT:


        -  ``i`` - integer between 0 and dimension of self,
           inclusive.


        OUTPUT: an affine scheme with fixed projective_embedding map.

        EXAMPLES::

            sage: PP = ProjectiveSpace(2, QQ, names='X,Y,Z')
            sage: X,Y,Z = PP.gens()
            sage: C = PP.subscheme(X^3*Y + Y^3*Z + Z^3*X)
            sage: U = C.affine_patch(0)
            sage: U
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
            x0^3*x1 + x1^3 + x0
            sage: U.projective_embedding()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x0^3*x1 + x1^3 + x0
              To:   Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              X^3*Y + Y^3*Z + X*Z^3
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1 : x0 : x1)
        """
        i = int(i)   # implicit type checking
        PP = self.ambient_space()
        n = PP.dimension()
        if i < 0 or i > n:
            raise ValueError, "Argument i (= %s) must be between 0 and %s."%(i, n)
        try:
            return self.__affine_patches[i]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        AA = PP.affine_patch(i)
        phi = AA.projective_embedding()
        polys = self.defining_polynomials()
        xi = phi.defining_polynomials()
        U = AA.subscheme([ f(xi) for f in polys ])
        U._default_embedding_index = i
        phi = U.projective_embedding(i, self)
        self.__affine_patches[i] = U
        return U


