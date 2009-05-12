"""
Affine `n` space over a ring.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import (
    is_FiniteField,
    is_RationalField,
    is_Ring,
    PolynomialRing,
    is_MPolynomialRing,
    ZZ,
    Integer)

from sage.misc.all import latex

import algebraic_scheme

import ambient_space

import homset

import morphism

import projective_space

import scheme

from sage.structure.parent_gens import normalize_names

def is_AffineSpace(x):
    r"""
    Returns True if x is an affine space, i.e., an ambient space
    `\mathbb{A}^n_R`, where `R` is a ring and
    `n\geq 0` is an integer.

    EXAMPLES::

        sage: from sage.schemes.generic.affine_space import is_AffineSpace
        sage: is_AffineSpace(AffineSpace(5, names='x'))
        True
        sage: is_AffineSpace(AffineSpace(5, GF(9,'alpha'), names='x'))
        True
        sage: is_AffineSpace(Spec(ZZ))
        False
    """
    return isinstance(x, AffineSpace_generic)

def AffineSpace(n, R=None, names='x'):
    r"""
    Return affine space of dimension `n` over the ring `R`.

    EXAMPLES: The dimension and ring can be given in either order.

    ::

        sage: AffineSpace(3, QQ, 'x')
        Affine Space of dimension 3 over Rational Field
        sage: AffineSpace(5, QQ, 'x')
        Affine Space of dimension 5 over Rational Field
        sage: A = AffineSpace(2, QQ, names='XY'); A
        Affine Space of dimension 2 over Rational Field
        sage: A.coordinate_ring()
        Multivariate Polynomial Ring in X, Y over Rational Field

    Use the divide operator for base extension.

    ::

        sage: AffineSpace(5, names='x')/GF(17)
        Affine Space of dimension 5 over Finite Field of size 17

    The default base ring is `\ZZ`.

    ::

        sage: AffineSpace(5, names='x')
        Affine Space of dimension 5 over Integer Ring

    There is also an affine space associated to each polynomial ring.

    ::

        sage: R = GF(7)['x,y,z']
        sage: A = AffineSpace(R); A
        Affine Space of dimension 3 over Finite Field of size 7
        sage: A.coordinate_ring() is R
        True
    """
    if is_MPolynomialRing(n) and R is None:
        R = n
        A = AffineSpace(R.ngens(), R.base_ring(), R.variable_names())
        A._coordinate_ring = R
        return A
    if isinstance(R, (int, long, Integer)):
        n, R = R, n
    if R is None:
        R = ZZ  # default is the integers
    if names is None:
        if n == 0:
            names = ''
        else:
            raise TypeError, "You must specify the variables names of the coordinate ring."
    return AffineSpace_generic(n, R, names)

class AffineSpace_generic(ambient_space.AmbientSpace, scheme.AffineScheme):
    """
    Affine space of dimension `n` over the ring `R`.

    EXAMPLES::

        sage: X.<x,y,z> = AffineSpace(3, QQ)
        sage: X.base_scheme()
        Spectrum of Rational Field
        sage: X.base_ring()
        Rational Field
        sage: X.structure_morphism ()
        Scheme morphism:
          From: Affine Space of dimension 3 over Rational Field
          To:   Spectrum of Rational Field
          Defn: Structure map

    Loading and saving::

        sage: loads(X.dumps()) == X
        True

    We create several other examples of affine spaces.

    ::

        sage: AffineSpace(5, PolynomialRing(QQ, 'z'), 'Z')
        Affine Space of dimension 5 over Univariate Polynomial Ring in z over Rational Field

    ::

        sage: AffineSpace(RealField(), 3, 'Z')
        Affine Space of dimension 3 over Real Field with 53 bits of precision

    ::

        sage: AffineSpace(Qp(7), 2, 'x')
        Affine Space of dimension 2 over 7-adic Field with capped relative precision 20

    Even 0-dimensional affine spaces are supported.

    ::

        sage: AffineSpace(0)
        Affine Space of dimension 0 over Integer Ring
    """
    def __init__(self, n, R, names):
        """
        EXAMPLES::

            sage: AffineSpace(3, Zp(5), 'y')
            Affine Space of dimension 3 over 5-adic Ring with capped relative precision 20
        """
        names = normalize_names(n, names)
        ambient_space.AmbientSpace.__init__(self, n, R)
        self._assign_names(names)

    def __iter__(self):
        """
        Return iterator over the elements of this affine space when defined
        over a finite field.

        EXAMPLES::

            sage: FF = FiniteField(3)
            sage: AA = AffineSpace(FF, 0)
            sage: [ x for x in AA ]
            [()]
            sage: AA = AffineSpace(FF, 1, 'Z')
            sage: [ x for x in AA ]
            [(0), (1), (2)]
            sage: AA.<z,w> = AffineSpace(FF, 2)
            sage: [ x for x in AA ]
            [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]

        AUTHOR:

        - David Kohel
        """
        n = self.dimension_relative()
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

    def ngens(self):
        """
        Return the number of generators of self, i.e. the number of
        variables in the coordinate ring of self.

        EXAMPLES::

            sage: AffineSpace(3, QQ).ngens()
            3
            sage: AffineSpace(7, ZZ).ngens()
            7
        """
        return self.dimension_relative()

    def rational_points(self, F=None):
        """
        Return the list of `F`-rational points on the affine space self,
        where `F` is a given finite field, or the base ring of self.

        EXAMPLES::

            sage: A = AffineSpace(1, GF(3))
            sage: A.rational_points()
            [(0), (1), (2)]
            sage: A.rational_points(GF(3^2, 'b'))
            [(0), (2*b), (b + 1), (b + 2), (2), (b), (2*b + 2), (2*b + 1), (1)]

            sage: AffineSpace(2, ZZ).rational_points(GF(2))
            [(0, 0), (1, 0), (0, 1), (1, 1)]

        TESTS::

            sage: AffineSpace(2, QQ).rational_points()
            Traceback (most recent call last):
            ...
            TypeError: Base ring (= Rational Field) must be a finite field.
            sage: AffineSpace(1, GF(3)).rational_points(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Second argument (= Integer Ring) must be a finite field.
        """
        if F == None:
            if not is_FiniteField(self.base_ring()):
                raise TypeError, "Base ring (= %s) must be a finite field."%self.base_ring()
            return [ P for P in self ]
        elif not is_FiniteField(F):
            raise TypeError, "Second argument (= %s) must be a finite field."%F
        return [ P for P in self.base_extend(F) ]

    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_affine_space(*args, **kwds)

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: AffineSpace(QQ, 3, 'a') == AffineSpace(ZZ, 3, 'a')
            False
            sage: AffineSpace(ZZ,1, 'a') == AffineSpace(ZZ, 0, 'a')
            False
            sage: loads(AffineSpace(ZZ, 1, 'x').dumps()) == AffineSpace(ZZ, 1, 'x')
            True
        """
        if not isinstance(right, AffineSpace_generic):
            return -1
        return cmp([self.dimension_relative(), self.coordinate_ring()],
                   [right.dimension_relative(), right.coordinate_ring()])

    def _latex_(self):
        r"""
        Return a LaTeX representation of this affine space.

        EXAMPLES::

            sage: print latex(AffineSpace(1, ZZ, 'x'))
            \mathbf{A}_{\Bold{Z}}^1

        TESTS::

            sage: AffineSpace(3, Zp(5), 'y')._latex_()
            '\\mathbf{A}_{\\ZZ_{5}}^3'
        """
        return "\\mathbf{A}_{%s}^%s"%(latex(self.base_ring()), self.dimension_relative())

    def _constructor(self, *args, **kwds):
        return AffineSpace(*args, **kwds)

    def _homset_class(self, *args, **kwds):
        return homset.SchemeHomset_affine_coordinates(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return morphism.SchemeMorphism_affine_coordinates(*args, **kwds)

    def _repr_(self):
        """
        Return a string representation of this affine space.

        EXAMPLES::

            sage: AffineSpace(1, ZZ, 'x')
            Affine Space of dimension 1 over Integer Ring

        TESTS::

            sage: AffineSpace(3, Zp(5), 'y')._repr_()
            'Affine Space of dimension 3 over 5-adic Ring with capped relative precision 20'
        """
        return "Affine Space of dimension %s over %s"%(self.dimension_relative(), self.base_ring())

    def _repr_generic_point(self, polys=None):
        """
        Return a string representation of the generic point
        corresponding to the list of polys on this affine space.

        If polys is None, the representation of the generic point of
        the affine space is returned.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: A._repr_generic_point([y-x^2])
            '(-x^2 + y)'
            sage: A._repr_generic_point()
            '(x, y)'
        """
        if polys is None:
            polys = self.gens()
        return '(%s)'%(", ".join([str(f) for f in polys]))

    def _latex_generic_point(self, v=None):
        """
        Return a LaTeX representation of the generic point
        corresponding to the list of polys on this affine space.

        If polys is None, the representation of the generic point of
        the affine space is returned.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: A._latex_generic_point([y-x^2])
            '\\left(- x^{2} + y\\right)'
            sage: A._latex_generic_point()
            '\\left(x, y\\right)'
        """
        if v is None:
            v = self.gens()
        return '\\left(%s\\right)'%(", ".join([str(latex(f)) for f in v]))

    def _check_satisfies_equations(self, v):
        """
        Return True if `v` defines a point on the scheme self; raise a
        TypeError otherwise.

        EXAMPLES::

            sage: A = AffineSpace(3, ZZ)
            sage: A._check_satisfies_equations([1, 1, 0])
            True
            sage: A._check_satisfies_equations((0, 1, 0))
            True
            sage: A._check_satisfies_equations([0, 0, 0])
            True
            sage: A._check_satisfies_equations([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            TypeError: The list v=[1, 2, 3, 4, 5] must have 3 components
            sage: A._check_satisfies_equations([1/2, 1, 1])
            Traceback (most recent call last):
            ...
            TypeError: The components of v=[1/2, 1, 1] must be elements of Integer Ring
            sage: A._check_satisfies_equations(5)
            Traceback (most recent call last):
            ...
            TypeError: The argument v=5 must be a list or tuple
        """
        if not isinstance(v, (list, tuple)):
            raise TypeError, 'The argument v=%s must be a list or tuple'%v
        n = self.ngens()
        if not len(v) == n:
            raise TypeError, 'The list v=%s must have %s components'%(v, n)
        R = self.base_ring()
        from sage.structure.sequence import Sequence
        if not Sequence(v).universe() == R:
            raise TypeError, 'The components of v=%s must be elements of %s'%(v, R)
        return True

    def __pow__(self, m):
        """
        EXAMPLES::

            sage: A = AffineSpace(1, QQ, 'x')
            sage: A^5
            Affine Space of dimension 5 over Rational Field
        """
        mm = int(m)
        if mm != m:
            raise ValueError, "m must be an integer"
        return self._constructor(self.dimension_relative() * mm, self._base_ring, names=self.variable_names() * mm)

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme, if defined.

        EXAMPLES::

            sage: R = AffineSpace(2, GF(9,'alpha'), 'z').coordinate_ring(); R
            Multivariate Polynomial Ring in z0, z1 over Finite Field in alpha of size 3^2
            sage: AffineSpace(3, R, 'x').coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Multivariate Polynomial Ring in z0, z1 over Finite Field in alpha of size 3^2
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            self._coordinate_ring = PolynomialRing(self.base_ring(), self.dimension_relative(), names=self.variable_names())
            return self._coordinate_ring

    def _validate(self, v):
        """
        Return a valid tuple of polynomial functions on self given by
        `v`.  Raise an error if `v` does not consist of valid
        functions.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: A._validate([x*y-z, 1])
            (x*y - z, 1)
            sage: A._validate([x, y, 1/3*z])
            Traceback (most recent call last):
            ...
            ValueError: The arguments [x, y, 1/3*z] are not valid polynomial functions on this affine space
        """
        R = self.coordinate_ring()
        try:
            return tuple([ R(g) for g in v ])
        except:
            raise ValueError, "The arguments %s are not valid polynomial functions on this affine space"%v

    def projective_embedding(self, i=None, PP=None):
        """
        Returns a morphism from this space into an ambient projective space
        of the same dimension.

        INPUT:


        -  ``i`` - integer (default: dimension of self = last
           coordinate) determines which projective embedding to compute. The
           embedding is that which has a 1 in the i-th coordinate, numbered
           from 0.

        -  ``PP`` - (default: None) ambient projective space, i.e.,
           codomain of morphism; this is constructed if it is not
           given.

        EXAMPLES::

            sage: AA = AffineSpace(2, QQ, 'x')
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
        n = self.dimension_relative()
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
        if n < 0 or n >self.dimension_relative():
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


        -  ``X`` - a list or tuple of equations


        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme([x, y^2, x*y^2]); X
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x
              y^2
              x*y^2

        ::

            sage: X.defining_polynomials ()
            (x, y^2, x*y^2)
            sage: I = X.defining_ideal(); I
            Ideal (x, y^2, x*y^2) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: I.groebner_basis()
            [y^2, x]
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
