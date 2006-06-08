"""
Affine $n$ space over a ring.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import (
    is_FiniteField,
    is_RationalField,
    is_Ring,
    MPolynomialRing,
    is_MPolynomialRing,
    ZZ)
from sage.misc.all import latex

import algebraic_scheme

import ambient_space

import homset

import morphism

import projective_space

import scheme

def is_AffineSpace(x):
    r"""
    Returns True if x is an affine space, i.e., an ambient space
    $\A^n_R$, where $R$ is a ring and $n\geq 0$ is an integer.

    EXAMPLES:
        sage: is_AffineSpace(AffineSpace(5))
        True
        sage: is_AffineSpace(AffineSpace(5, GF(9)))
        True
        sage: is_AffineSpace(Spec(ZZ))
        False
    """
    return isinstance(x, AffineSpace_generic)

def AffineSpace(n, R=None, names=None):
    r"""
    Return affine space of dimension $n$ over the ring $R$.

    EXAMPLES:
    The dimension and ring can be given in either order.
        sage: AffineSpace(3, QQ)
        Affine Space of dimension 3 over Rational Field
        sage: AffineSpace(5, QQ)
        Affine Space of dimension 5 over Rational Field
        sage: A = AffineSpace(2, QQ, names='XY'); A
        Affine Space of dimension 2 over Rational Field
        sage: A.coordinate_ring()
        Polynomial Ring in X, Y over Rational Field

    Use the divide operator for base extension.
        sage: AffineSpace(5)/GF(17)
        Affine Space of dimension 5 over Finite Field of size 17

    The default base ring is $\Z$.
        sage: AffineSpace(5)
        Affine Space of dimension 5 over Integer Ring

    There is also an affine space associated each polynomial ring.
        sage: R = GF(7)['x,y,z']
        sage: A = AffineSpace(R); A
        Affine Space of dimension 3 over Finite Field of size 7
        sage: A.coordinate_ring() is R
        True
    """
    if is_MPolynomialRing(n) and R is None:
        A = AffineSpace(n.ngens(), n.base_ring())
        A._coordinate_ring = n
        return A
    if R is None:
        R = ZZ  # default is the integers
    return AffineSpace_generic(n, R, names)

class AffineSpace_generic(ambient_space.AmbientSpace, scheme.AffineScheme):
    """
    Affine space of dimension $n$ over the ring $R$.

    EXAMPLES:
        sage: X = AffineSpace(3, QQ)
        sage: X.base_scheme()
        Spectrum of Rational Field
        sage: X.base_ring()
        Rational Field
        sage: X.structure_morphism ()
        Scheme morphism:
          From: Affine Space of dimension 3 over Rational Field
          To:   Spectrum of Rational Field
          Defn: Structure map

    Loading and saving:
        sage: loads(X.dumps()) == X
        True

    We create several other examples of affine spaces.
        sage: AffineSpace(5, PolynomialRing(QQ, 'z'))
        Affine Space of dimension 5 over Univariate Polynomial Ring in z over Rational Field

        sage: AffineSpace(3, RealField())
        Affine Space of dimension 3 over Real Field with 53 bits of precision

        sage: AffineSpace(2, pAdicField(7))
        Affine Space of dimension 2 over 7-adic Field

    Even 0-dimensional affine spaces are supported.
        sage: AffineSpace(0)
        Affine Space of dimension 0 over Integer Ring
    """
    def __init__(self, n, R, names=None):
        ambient_space.AmbientSpace.__init__(self, n, R)
        self.__names = names

    def __iter__(self):
        """
        Return iterator over the elements of this affine space when defined
        over a finite field.

        EXAMPLES:
            sage: FF = FiniteField(3)
            sage: AA = AffineSpace(0, FF)
            sage: [ x for x in AA ]
            [()]
            sage: AA = AffineSpace(1, FF)
            sage: [ x for x in AA ]
            [(0), (1), (2)]
            sage: AA = AffineSpace(2, FF)
            sage: [ x for x in AA ]
            [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]

        AUTHOR: David Kohel <kohel@maths.usyd.edu.au>
        """
        n = self.dimension()
        R = self.base_ring()
        zero = R(0)
        P = [ zero for _ in range(n) ]
        yield self(P)
        iters = [ iter(R) for _ in range(n) ]
        for x in iters: x.next() # put at zero
        i = 0
        while i < n:
            try:
                P[i] = iters[i].next()
                yield self(P)
                i = 0
            except StopIteration:
                iters[i] = iter(R)  # reset
                iters[i].next() # put at zero
                P[i] = zero
                i += 1

    def rational_points(self, F=None):
        if F == None:
            if not is_FiniteField(self.base_ring()):
                raise TypeError, "Base ring (= %s) must be a finite field."%self.base_ring()
            return [ P for P in self ]
        elif not is_FiniteField(F):
            raise TypeError, "Second argument (= %s) must be a finite field."%F
        return [ P for P in self(F) ]

    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_affine_space(*args, **kwds)

    def __cmp__(self, right):
        """
        EXAMPLES:
            sage: AffineSpace(3, QQ) == AffineSpace(3, ZZ)
            False
            sage: AffineSpace(1, ZZ) == AffineSpace(0, ZZ)
            False
            sage: loads(AffineSpace(1, ZZ).dumps()) == AffineSpace(1, ZZ)
            True
        """
        if not isinstance(right, AffineSpace_generic):
            return -1
        return cmp([self.dimension(), self.coordinate_ring()],
                   [right.dimension(), right.coordinate_ring()])

    def _latex_(self):
        r"""
        EXAMPLES:
            sage: print latex(AffineSpace(1, ZZ))
            {\mathbf A}_{\mbox{\bf{}Z}}^1
        """
        return "{\\mathbf A}_{%s}^%s"%(latex(self.base_ring()), self.dimension())

    def _constructor(self, *args, **kwds):
        return AffineSpace(*args, **kwds)

    def _homset_class(self, *args, **kwds):
        return homset.SchemeHomset_affine_coordinates(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return morphism.SchemeMorphism_affine_coordinates(*args, **kwds)

    def _repr_(self):
        return "Affine Space of dimension %s over %s"%(self.dimension(), self.base_ring())

    def _repr_generic_point(self, polys=None):
        if polys is None:
            polys = self.gens()
        return '(%s)'%(", ".join([str(f) for f in polys]))

    def _latex_generic_point(self, v=None):
        if v is None:
            v = self.gens()
        return '\\left(%s\\right)'%(", ".join([str(latex(f)) for f in v]))

    def __pow__(self, m):
        """
        EXAMPLES:
            sage: A = AffineSpace(1, QQ)
            sage: A^5
            Affine Space of dimension 5 over Rational Field
        """
        return self._constructor(self.dimension() * m, self._base_ring)

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme, if defined.  Otherwise raise
        a ValueError.

        EXAMPLES:
            sage: R = AffineSpace(2, GF(9)).coordinate_ring(); R
            Polynomial Ring in x0, x1 over Finite Field in a of size 3^2
            sage: AffineSpace(3, R).coordinate_ring()
            Polynomial Ring in x0, x1, x2 over Polynomial Ring in x0, x1 over Finite Field in a of size 3^2
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            self._coordinate_ring = MPolynomialRing(self.base_ring(), self.dimension(), names=self.__names)
            return self._coordinate_ring

    def projective_embedding(self, i=None, PP=None):
        """
        Returns a morphism from this space into an ambient projective space of
        the same dimension.

        INPUT:
            i -- integer (default: dimension of self = last coordinate) determines
                 which projective embedding to compute.  The embedding is that
                 which has a 1 in the i-th coordinate, numbered from 0.

            PP -- (default: None) ambient projective space, i.e., codomain of morphism;
                 this is constructed if it is not given.

        EXAMPLES:
            sage: AA = AffineSpace(2, QQ)
            sage: pi = AA.projective_embedding(0); pi
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1 : x0 : x1)
            sage: z = AA(3,4)
            sage: pi(z)
            (1/4 : 3/4 : 1)
            sage: pi(AA(0,2))
            (1/2 : 0 : 1)
            sage: pi = AA.projective_embedding(1); pi
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (x0 : 1 : x1)
            sage: pi(z)
            (3/4 : 1/4 : 1)
            sage: pi = AA.projective_embedding(2)
            sage: pi(z)
            (3 : 4 : 1)
        """
        n = self.dimension()
        if i is None:
            try:
                i = self._default_embedding_index
            except AttributeError:
                i = int(n)
        else:
            i = int(i)
        try:
            return self.__projective_embedding[i]
        except AttributeError:
            self.__projective_embedding = {}
        except KeyError:
            pass
        if PP is None:
            PP = projective_space.ProjectiveSpace(n, self.base_ring())
        R = self.coordinate_ring()
        v = list(R.gens())
        if n < 0 or n >self.dimension():
            raise ValueError, \
                  "Argument i (=%s) must be between 0 and %s, inclusive"%(i,n)
        v.insert(i, R(1))
        phi = self.hom(v, PP)
        self.__projective_embedding[i] = phi
        return phi

    def subscheme(self, X):
        """
        Return the closed subscheme defined by X.

        INPUT:
            X -- a list or tuple of equations

        EXAMPLES:
            sage: A, (x,y) = AffineSpace(2, QQ).objgens('xy')
            sage: X = A.subscheme([x, y^2, x*y^2]); X
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x
              y^2
              x*y^2

            sage: X.defining_polynomials ()
            (x, y^2, x*y^2)
            sage: I = X.defining_ideal(); I
            Ideal (y^2, x, x*y^2) of Polynomial Ring in x, y over Rational Field
            sage: I.groebner_basis()
            [x, y^2]
            sage: X.dimension()
            0
            sage: X.base_ring()
            Rational Field
            sage: X.base_scheme()
            Spectrum of Rational Field
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x
              y^2
              x*y^2
              To:   Spectrum of Rational Field
              Defn: Structure map
            sage: X.dimension()
            0
        """
        return algebraic_scheme.AlgebraicScheme_subscheme_affine(self, X)

    def subscheme_complement(self, X, Y):
        X = self.subscheme(X)
        Y = self.subscheme(Y)
        return algebraic_scheme.AlgebraicScheme_quasi(X, Y)
