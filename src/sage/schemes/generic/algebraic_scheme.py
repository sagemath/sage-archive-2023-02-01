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

from sage.rings.all import is_Ideal, is_MPolynomialRing

from sage.structure.all import Sequence

import ambient_space

import morphism

import scheme

def is_AlgebraicScheme(x):
    """
    Return True if $x$ is an algebraic scheme, i.e., a subscheme of an
    ambient space over a ring defined by polynomial equations.

    EXAMPLES:
    Affine space is itself not an algebraic scheme, though the closed subscheme
    defined by no equations is.
        sage: is_AlgebraicScheme(AffineSpace(10, Q))
        False
        sage: V = AffineSpace(10, Q).subscheme([]); V
        Closed subscheme of Affine Space of dimension 10 over
        Rational Field defined by:
          (no equations)
        sage: is_AlgebraicScheme(V)
        True

    We create a more complicated closed subscheme.
        sage: A, x = AffineSpace(10, Q).objgens()
        sage: X = A.subscheme([sum(x)]); X
        Closed subscheme of Affine Space of dimension 10 over Rational Field defined by:
          x_9 + x_8 + x_7 + x_6 + x_5 + x_4 + x_3 + x_2 + x_1 + x_0
        sage: is_AlgebraicScheme(X)
        True

        sage: is_AlgebraicScheme(Q)
        False
        sage: S = Spec(Q)
        sage: is_AlgebraicScheme(S)
        False
    """
    return isinstance(x, AlgebraicScheme)

class AlgebraicScheme(scheme.Scheme):
    """
    An algebraic scheme presented as a subscheme in an ambient space.
    """
    def __init__(self, A):
        if not ambient_space.is_AmbientSpace(A):
            raise TypeError, "A (=%s) must be an ambient space"
        self.__A = A

    def ambient_space(self):
        return self.__A

    def degree(self):
        return self.__A.degree()

    def _repr_(self):
        return "An subscheme of %s"%self.__A

    def _homset_class(self, *args, **kwds):
        return self.__A._homset_class(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return self.__A._point_class(*args, **kwds)


class AlgebraicScheme_quasi(AlgebraicScheme):
    """
    The quasi-affine or quasi-projective scheme X - Y, where X and Y
    are both closed subschemes of a common ambient affine or
    projective space.
    """
    def __init__(self, X, Y):
        self.__X = X
        self.__Y = Y
        if not isinstance(X, AlgebraicScheme_subscheme):
            raise TypeError, "X must be a closed subscheme of an ambient space."
        if not isinstance(Y, AlgebraicScheme_subscheme):
            raise TypeError, "Y must be a closed subscheme of an ambient space."
        if X.ambient_space() != Y.ambient_space():
            raise ValueError, "X and Y must be embedded in the same ambient space."
        AlgebraicScheme.__init__(self, X.ambient_space())

    def _repr_(self):
        if affine_space.is_AffineSpace(A):
            t = "affine"
        else:
            t = "projective"
        s =  "Quasi-%s scheme X - Y, where:\n"%t
        s += "  X: %s\n"%self.__X
        s += "  Y: %s\n"%self.__Y
        return s

    def X(self):
        return self.__X

    def Y(self):
        return self.__Y


class AlgebraicScheme_subscheme(AlgebraicScheme):
    """
    An algebraic scheme presented as a closed subscheme is defined by
    explicit polynomial equations.  This is as opposed to a general
    scheme, which could, e.g., by the Neron model of some object, and
    for which we do not want to give explicit equations.

    INPUT:
        A -- ambient space (affine or projective n-space over a ring)
        G -- ideal or tuple of defining polynomials
    """
    def __init__(self, A, G):
        AlgebraicScheme.__init__(self, A)
        self._base_ring = A.base_ring()
        if is_Ideal(G):
            self.__I = G
            G = G.gens()
        if not isinstance(G, (list, tuple)):
            G = (G, )
        else:
            G = tuple(G)
        if len(G) > 0:
            G = self._validate(G)
        self.__G = G  # now G is a tuple of defining polynomials.

    def _validate(self, G):
        G = Sequence(G)  # make sure they have a common parent
        if not is_MPolynomialRing(G.universe()):
            raise TypeError, "each generator must be a multivariate polynomial"
        return tuple(G)

    def _error_bad_coords(self, v):
        raise TypeError, "coordinates %s do not define a point on %s"%(v,self)

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme,
        or raise a TypeError.
        """
        for f in self.defining_polynomials():
            if f(v) != 0:
                self._error_bad_coords(v)

    def base_extend(self, R):
        A = self.ambient_space().base_extend(R)
        return A.subscheme(self.__G)

    def __cmp__(self, other):
        if not isinstance(other, AlgebraicScheme_subscheme):
            return -1
        A = self.ambient_space()
        if other.ambient_space() != A:
            return -1
        return cmp(self.defining_ideal(), other.defining_ideal())

    def _repr_(self):
        polys = '\n  '.join([str(f) for f in self.defining_polynomials()])
        if polys == '':
            polys = '(no equations)'
        return "Closed subscheme of %s defined by:\n  %s"%(self.ambient_space(),
                                                           polys)

    def coordinate_ring(self):
        try:
            return self._coordinate_ring
        except AttributeError:
            R = self.__A.coordinate_ring()
            if len(self.__X) == 0:
                Q = R
            else:
                I = self.defining_ideal()
                Q = R/I
            self._coordinate_ring = Q
            return Q

    def defining_polynomials(self):
        return self.__G

    def defining_ideal(self):
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

        OUTPUT:
            an immutable sequence of irreducible subschemes of the ambient
            space of this scheme

        The components are cached.

        EXAMPLES:

        We define what is clearly a union of four hypersurfaces in
        $\P^4_{\Q}$ then find the irreducible components.
            sage: P, (x,y,z,w,v) = ProjectiveSpace(4,QQ).objgens()
            sage: V = P.subscheme( (x^2 - y^2 - z^2)*(w^5 -  2*v^2*z^3)* w * (v^3 - x^2*z) )
            sage: V.irreducible_components()
            [Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
                      -1*x_4^3 + x_0^2*x_2,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
                      -1*x_2^2 - x_1^2 + x_0^2,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
                      x_3,
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
                      -1*x_3^5 + 2*x_2^3*x_4^2]
        """
        try:
            return self.__irreducible_components
        except AttributeError:
            pass
        I = self.defining_ideal()
        P = I.associated_primes()
        A = self.ambient_space()
        C = Sequence([A.subscheme(X) for X in P], check=False, immutable=True)
        self.__irreducible_components = C
        return C

    def reduce(self):
        r"""
        Return the corresponding reduced algebraic space associated to
        this scheme.

        EXAMPLES:
        First we construct the union of a doubled and triplled line
        in the affine plane over $\Q$.
            sage: A, (x,y) = AffineSpace(2, Q).objgens('xy')
            sage: X = A.subscheme([(x-1)^2*(x-y)^3]); X
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -1*y^3 + 3*x*y^2 + 2*x*y^3 - 3*x^2*y - 6*x^2*y^2 - x^2*y^3 + x^3 + 6*x^3*y + 3*x^3*y^2 - 2*x^4 - 3*x^4*y + x^5
            sage: X.dimension()
            1

        Then we compute the corresponding reduced scheme.
            sage: Y = X.reduce(); Y
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              y - x - x*y + x^2

        Finally, we verify that the reduced scheme $Y$ is the union of those two lines.
            sage: L1 = A.subscheme([x-1]); L1
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -1 + x
            sage: L2 = A.subscheme([x-y]); L2
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -1*y + x
            sage: W = L1.union(L2); W             # taken in ambient space
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              y - x - x*y + x^2
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

        EXAMPLES:
        We construct the union of a line and a tripled-point on the line.
            sage: A, (x,y) = AffineSpace(2, Q, 'xy').objgens()
            sage: I = ideal([x,y])^3
            sage: P = A.subscheme(I)
            sage: L = A.subscheme([y-1])
            sage: S = L.union(P); S
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -1*x*y^2 + x*y^3
              -1*y^3 + y^4
              -1*x^2*y + x^2*y^2
              -1*x^3 + x^3*y
            sage: S.dimension()
            1
            sage: S.reduce()
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -1*x + x*y
              -1*y + y^2

        We can also use the notation "+" for the union:
            sage: A.subscheme([x]) + A.subscheme([y^2 - (x^3+1)])
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -1*x + x*y^2 - x^4

        Saving and loading:
            sage: loads(S.dumps()) == S
            True
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError, "other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
        return A.subscheme(self.defining_ideal().intersection(other.defining_ideal()))


    def intersection(self, other):
        """
        Return the scheme-theoretic intersection of self and other in their common
        ambient space.

        EXAMPLES:
        """
        if not isinstance(other, AlgebraicScheme_subscheme):
            raise TypeError, "other (=%s) must be a closed algebraic subscheme of an ambient space"%other
        A = self.ambient_space()
        if other.ambient_space() != A:
            raise ValueError, "other (=%s) must be in the same ambient space as self"%other
        return A.subscheme(self.defining_ideal() + other.defining_ideal())


class AlgebraicScheme_subscheme_affine(AlgebraicScheme_subscheme):
    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_affine_space(*args, **kwds)

    def dimension(self):
        """
        EXAMPLES:
            sage: A, (x,y) = AffineSpace(2, Q).objgens('xy')
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
            sage: A, (x,y,z,w) = AffineSpace(4, Q).objgens('xyzw')
            sage: X = A.subscheme([x^2, x^2*y^2 + z^2, z^2 - w^2, 10*x^2 + w^2 - z^2])
            sage: X
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x^2
              z^2 + x^2*y^2
              -1*w^2 + z^2
              w^2 - z^2 + 10*x^2
            sage: X.dimension()
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension()
            return self.__dimension


class AlgebraicScheme_subscheme_projective(AlgebraicScheme_subscheme):
    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_projective_space(*args, **kwds)

    def _validate(self, G):
        G = AlgebraicScheme_subscheme._validate(self, G)
        for f in G:
            if not f.is_homogeneous():
                raise TypeError, \
                      "defining polynomials (= %s) must be homogeneous"%G
        return G

    def dimension(self):
        """
        EXAMPLES:
            sage: A, (x,y) = AffineSpace(2, Q).objgens('xy')
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
            sage: A, (x,y,z,w) = AffineSpace(4, Q).objgens('xyzw')
            sage: X = A.subscheme([x^2, x^2*y^2 + z^2, z^2 - w^2, 10*x^2 + w^2 - z^2])
            sage: X
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x^2
              z^2 + x^2*y^2
              -1*w^2 + z^2
              w^2 - z^2 + 10*x^2
            sage: X.dimension()
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension() - 1
            return self.__dimension

