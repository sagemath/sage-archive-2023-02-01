r"""
Projective `n` space over a ring.

EXAMPLES: We construct projective space over various rings of
various dimensions.

The simplest projective space::

    sage: ProjectiveSpace(0)
    Projective Space of dimension 0 over Integer Ring

A slightly bigger projective space over `\QQ`::

    sage: X = ProjectiveSpace(1000, QQ); X
    Projective Space of dimension 1000 over Rational Field
    sage: X.dimension()
    1000

We can use "over" notation to create projective spaces over various
base rings.

::

    sage: X = ProjectiveSpace(5)/QQ; X
    Projective Space of dimension 5 over Rational Field
    sage: X/CC
    Projective Space of dimension 5 over Complex Field with 53 bits of precision

The third argument specifies the printing names of the generators
of the homogenous coordinate ring. Using objgens() you can obtain
both the space and the generators as ready to use variables.

::

    sage: P2, (x,y,z) = ProjectiveSpace(2, QQ, 'xyz').objgens()
    sage: P2
    Projective Space of dimension 2 over Rational Field
    sage: x.parent()
    Multivariate Polynomial Ring in x, y, z over Rational Field

For example, we use `x,y,z` to define the intersection of
two lines.

::

    sage: V = P2.subscheme([x+y+z, x+y-z]); V
    Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
     x + y + z
     x + y - z
    sage: V.dimension()
    0
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import (PolynomialRing,
                            is_Field,
                            is_FiniteField,
                            is_RationalField,
                            is_Ring,
                            is_CommutativeRing,
                            is_MPolynomialRing,
                            Integer,
                            ZZ)
from sage.modules.all import VectorSpace

from sage.misc.all import latex

import homset
import morphism

import algebraic_scheme
import ambient_space
import affine_space

from sage.structure.parent_gens import normalize_names

from sage.rings.arith import gcd

from sage.combinat.tuple import Tuples

def is_ProjectiveSpace(x):
    r"""
    Return True if `x` is a projective space, i.e., an ambient space
    `\mathbb{P}^n_R`, where `R` is a ring and `n\geq 0` is an
    integer.

    EXAMPLES::

        sage: from sage.schemes.generic.projective_space import is_ProjectiveSpace
        sage: is_ProjectiveSpace(ProjectiveSpace(5, names='x'))
        True
        sage: is_ProjectiveSpace(ProjectiveSpace(5, GF(9,'alpha'), names='x'))
        True
        sage: is_ProjectiveSpace(Spec(ZZ))
        False
    """
    return isinstance(x, ProjectiveSpace_ring)

def ProjectiveSpace(n, R=None, names='x'):
    r"""
    Return projective space of dimension `n` over the ring `R`.

    EXAMPLES: The dimension and ring can be given in either order.

    ::

        sage: ProjectiveSpace(3, QQ)
        Projective Space of dimension 3 over Rational Field
        sage: ProjectiveSpace(5, QQ)
        Projective Space of dimension 5 over Rational Field
        sage: P = ProjectiveSpace(2, QQ, names='XYZ'); P
        Projective Space of dimension 2 over Rational Field
        sage: P.coordinate_ring()
        Multivariate Polynomial Ring in X, Y, Z over Rational Field

    The divide operator does base extension.

    ::

        sage: ProjectiveSpace(5)/GF(17)
        Projective Space of dimension 5 over Finite Field of size 17

    The default base ring is `\ZZ`.

    ::

        sage: ProjectiveSpace(5)
        Projective Space of dimension 5 over Integer Ring

    There is also an projective space associated each polynomial ring.

    ::

        sage: R = GF(7)['x,y,z']
        sage: P = ProjectiveSpace(R); P
        Projective Space of dimension 2 over Finite Field of size 7
        sage: P.coordinate_ring()
        Multivariate Polynomial Ring in x, y, z over Finite Field of size 7
        sage: P.coordinate_ring() is R
        True

    Projective spaces are not cached, i.e., there can be several with
    the same base ring and dimension (to facilitate glueing
    constructions).
    """
    if is_MPolynomialRing(n) and R is None:
        A = ProjectiveSpace(n.ngens()-1, n.base_ring())
        A._coordinate_ring = n
        return A
    if isinstance(R, (int, long, Integer)):
        n, R = R, n
    if R is None:
        R = ZZ  # default is the integers
    if is_Field(R):
        if is_FiniteField(R):
            return ProjectiveSpace_finite_field(n, R, names)
        if is_RationalField(R):
            return ProjectiveSpace_rational_field(n, R, names)
        else:
            return ProjectiveSpace_field(n, R, names)
    elif is_CommutativeRing(R):
        return ProjectiveSpace_ring(n, R, names)
    else:
        raise TypeError, "R (=%s) must be a commutative ring"%R

class ProjectiveSpace_ring(ambient_space.AmbientSpace):
    """
    Projective space of dimension `n` over the ring
    `R`.

    EXAMPLES::

        sage: X.<x,y,z,w> = ProjectiveSpace(3, QQ)
        sage: X.base_scheme()
        Spectrum of Rational Field
        sage: X.base_ring()
        Rational Field
        sage: X.structure_morphism ()
        Scheme morphism:
          From: Projective Space of dimension 3 over Rational Field
          To:   Spectrum of Rational Field
          Defn: Structure map
        sage: X.coordinate_ring()
        Multivariate Polynomial Ring in x, y, z, w over Rational Field

    Loading and saving::

        sage: loads(X.dumps()) == X
        True
    """
    def __init__(self, n, R=ZZ, names=None):
        """
        EXAMPLES::

            sage: ProjectiveSpace(3, Zp(5), 'y')
            Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20
        """
        names = normalize_names(n+1, names)
        ambient_space.AmbientSpace.__init__(self, n, R)
        self._assign_names(names)

    def ngens(self):
        """
        Return the number of generators of self, i.e. the number of
        variables in the coordinate ring of self.

        EXAMPLES::

            sage: ProjectiveSpace(3, QQ).ngens()
            4
            sage: ProjectiveSpace(7, ZZ).ngens()
            8
        """
        return self.dimension_relative() + 1

    def _check_satisfies_equations(self, v):
        """
        Return True if `v` defines a point on the scheme self; raise a
        TypeError otherwise.

        EXAMPLES::

            sage: P = ProjectiveSpace(2, ZZ)
            sage: P._check_satisfies_equations([1, 1, 0])
            True
            sage: P._check_satisfies_equations((0, 1, 0))
            True
            sage: P._check_satisfies_equations([0, 0, 0])
            Traceback (most recent call last):
            ...
            TypeError: The zero vector is not a point in projective space
            sage: P._check_satisfies_equations([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            TypeError: The list v=[1, 2, 3, 4, 5] must have 3 components
            sage: P._check_satisfies_equations([1/2, 1, 1])
            Traceback (most recent call last):
            ...
            TypeError: The components of v=[1/2, 1, 1] must be elements of Integer Ring
            sage: P._check_satisfies_equations(5)
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
        for coord in v:
            if not coord in R:
                raise TypeError, 'The components of v=%s must be elements of %s'%(v, R)
        zero = [R(0)]*n
        if v == zero:
            raise TypeError, 'The zero vector is not a point in projective space'
        return True

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme.

        EXAMPLES::

            sage: ProjectiveSpace(3, GF(19^2,'alpha'), 'abcd').coordinate_ring()
            Multivariate Polynomial Ring in a, b, c, d over Finite Field in alpha of size 19^2

        ::

            sage: ProjectiveSpace(3).coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring

        ::

            sage: ProjectiveSpace(2, QQ, ['alpha', 'beta', 'gamma']).coordinate_ring()
            Multivariate Polynomial Ring in alpha, beta, gamma over Rational Field
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            self._coordinate_ring = PolynomialRing(self.base_ring(),
                               self.variable_names(), self.dimension_relative()+1)
            return self._coordinate_ring

    def _validate(self, v):
        """
        Return a valid tuple of polynomial functions on self given by
        `v`.  Raise an error if `v` does not consist of valid
        functions (in particular, if they are not homogeneous).

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._validate([x*y-z^2, 1])
            (x*y - z^2, 1)
            sage: P._validate([x, y, 1/3*z])
            Traceback (most recent call last):
            ...
            ValueError: The arguments [x, y, 1/3*z] are not valid polynomial functions on this projective space
            sage: P._validate([x*y-z])
            Traceback (most recent call last):
            ...
            TypeError: The polynomial(s) [x*y - z] must be homogeneous
        """
        R = self.coordinate_ring()
        try:
            tup = tuple([ R(g) for g in v ])
        except:
            raise ValueError, "The arguments %s are not valid polynomial functions on this projective space"%v
        for g in tup:
            if not g.is_homogeneous():
                raise TypeError, \
                      "The polynomial(s) %s must be homogeneous"%v
        return tup

    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_projective_space(*args, **kwds)

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: ProjectiveSpace(QQ, 3, 'a') == ProjectiveSpace(ZZ, 3, 'a')
            False
            sage: ProjectiveSpace(ZZ, 1, 'a') == ProjectiveSpace(ZZ, 0, 'a')
            False
            sage: ProjectiveSpace(ZZ, 2, 'a') == AffineSpace(ZZ, 2, 'a')
            False
            sage: loads(AffineSpace(ZZ, 1, 'x').dumps()) == AffineSpace(ZZ, 1, 'x')
            True
        """
        if not isinstance(right, ProjectiveSpace_ring):
            return -1
        return cmp([self.dimension_relative(), self.coordinate_ring()],
                   [right.dimension_relative(), right.coordinate_ring()])

    def _latex_(self):
        r"""
        Return a LaTeX representation of this projective space.

        EXAMPLES::

            sage: print latex(ProjectiveSpace(1, ZZ, 'x'))
            {\mathbf P}_{\Bold{Z}}^1

        TESTS::

            sage: ProjectiveSpace(3, Zp(5), 'y')._latex_()
            '{\\mathbf P}_{\\ZZ_{5}}^3'
        """
        return "{\\mathbf P}_{%s}^%s"%(latex(self.base_ring()), self.dimension_relative())

    def _constructor(self, *args, **kwds):
        return ProjectiveSpace(*args, **kwds)

    def _homset_class(self, *args, **kwds):
        return homset.SchemeHomset_projective_coordinates_ring(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return morphism.SchemeMorphism_projective_coordinates_ring(*args, **kwds)

    def _repr_(self):
        """
        Return a string representation of this projective space.

        EXAMPLES::

            sage: ProjectiveSpace(1, ZZ, 'x')
            Projective Space of dimension 1 over Integer Ring

        TESTS::

            sage: ProjectiveSpace(3, Zp(5), 'y')._repr_()
            'Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20'
        """
        return "Projective Space of dimension %s over %s"%(self.dimension_relative(), self.base_ring())

    def _repr_generic_point(self, v=None):
        """
        Return a string representation of the generic point
        corresponding to the list of polys on this projective space.

        If polys is None, the representation of the generic point of
        the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._repr_generic_point([z*y-x^2])
            '(-x^2 + y*z)'
            sage: P._repr_generic_point()
            '(x : y : z)'
        """
        if v is None:
            v = self.gens()
        return '(%s)'%(" : ".join([repr(f) for f in v]))

    def _latex_generic_point(self, v=None):
        """
        Return a LaTeX representation of the generic point
        corresponding to the list of polys on this projective space.

        If polys is None, the representation of the generic point of
        the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P._latex_generic_point([z*y-x^2])
            '\\left(- x^{2} + y z\\right)'
            sage: P._latex_generic_point()
            '\\left(x : y : z\\right)'
        """
        if v is None:
            v = self.gens()
        return '\\left(%s\\right)'%(" : ".join([str(latex(f)) for f in v]))


    def subscheme(self, X):
        """
        Return the closed subscheme defined by X.

        INPUT:


        -  ``X`` - a list or tuple of equations


        EXAMPLES::

            sage: A.<x,y,z> = ProjectiveSpace(2, QQ)
            sage: X = A.subscheme([x*z^2, y^2*z, x*y^2]); X
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*z^2
              y^2*z
              x*y^2
            sage: X.defining_polynomials ()
            (x*z^2, y^2*z, x*y^2)
            sage: I = X.defining_ideal(); I
            Ideal (x*z^2, y^2*z, x*y^2) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: I.groebner_basis()
            [x*y^2, y^2*z,  x*z^2]
            sage: X.dimension()
            0
            sage: X.base_ring()
            Rational Field
            sage: X.base_scheme()
            Spectrum of Rational Field
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*z^2
              y^2*z
              x*y^2
              To:   Spectrum of Rational Field
              Defn: Structure map
        """
        return algebraic_scheme.AlgebraicScheme_subscheme_projective(self, X)

    def affine_patch(self, i):
        r"""
        Return the `i^{th}` affine patch of this projective space.
        This is an ambient affine space `\mathbb{A}^n_R,` where
        `R` is the base ring of self, whose "projective embedding"
        map is `1` in the `i^{th}` factor.

        INPUT:


        -  ``i`` - integer between 0 and dimension of self,
           inclusive.


        OUTPUT: an ambient affine space with fixed projective_embedding
        map.

        EXAMPLES::

            sage: PP = ProjectiveSpace(5) / QQ
            sage: AA = PP.affine_patch(2)
            sage: AA
            Affine Space of dimension 5 over Rational Field
            sage: AA.projective_embedding()
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1, x2, x3, x4) to
                    (x0 : x1 : 1 : x2 : x3 : x4)
            sage: AA.projective_embedding(0)
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1, x2, x3, x4) to
                    (1 : x0 : x1 : x2 : x3 : x4)
        """
        i = int(i)   # implicit type checking
        n = self.dimension_relative()
        if i < 0 or i > n:
            raise ValueError, "Argument i (= %s) must be between 0 and %s."%(i, n)
        try:
            return self.__affine_patches[i]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        AA = affine_space.AffineSpace(n, self.base_ring(), names='x')
        AA._default_embedding_index = i
        phi = AA.projective_embedding(i, self)
        self.__affine_patches[i] = AA
        return AA


class ProjectiveSpace_field(ProjectiveSpace_ring):
    def _homset_class(self, *args, **kwds):
        return homset.SchemeHomset_projective_coordinates_field(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return morphism.SchemeMorphism_projective_coordinates_field(*args, **kwds)


class ProjectiveSpace_finite_field(ProjectiveSpace_field):
    def __iter__(self):
        r"""
        Return iterator over the elements of this projective space.

        Note that iteration is over the decomposition
        `\PP^n = \mathbb{A}A^n \cup \PP^n-1`, where
        `\mathbb{A}A^n` is the `n`-th affine patch and
        `\PP^n-1` is the hyperplane at infinity
        `x_n = 0`.

        EXAMPLES::

            sage: FF = FiniteField(3)
            sage: PP = ProjectiveSpace(0,FF)
            sage: [ x for x in PP ]
            [(1)]
            sage: PP = ProjectiveSpace(1,FF)
            sage: [ x for x in PP ]
            [(0 : 1), (1 : 1), (2 : 1), (1 : 0)]
            sage: PP = ProjectiveSpace(2,FF)
            sage: [ x for x in PP ]
            [(0 : 0 : 1),
             (1 : 0 : 1),
             (2 : 0 : 1),
             (0 : 1 : 1),
             (1 : 1 : 1),
             (2 : 1 : 1),
             (0 : 2 : 1),
             (1 : 2 : 1),
             (2 : 2 : 1),
             (0 : 1 : 0),
             (1 : 1 : 0),
             (2 : 1 : 0),
             (1 : 0 : 0)]

        AUTHORS:

        - David Kohel

        TODO: Iteration for point sets over finite fields, and return of
        iter of point set over base field. Note that the point set does not
        know whether this is a projective space or subscheme.
        """
        n = self.dimension_relative()
        R = self.base_ring()
        zero = R(0)
        i = n
        while not i < 0:
            P = [ zero for _ in range(i) ] + [ R(1) ] + [ zero for _ in range(n-i) ]
            yield self(P)
            iters = [ iter(R) for _ in range(i) ]
            for x in iters: x.next() # put at zero
            j = 0
            while j < i:
                try:
                    P[j] = iters[j].next()
                    yield self(P)
                    j = 0
                except StopIteration:
                    iters[j] = iter(R)  # reset
                    iters[j].next() # put at zero
                    P[j] = zero
                    j += 1
            i -= 1

    def rational_points(self, F=None):
        """
        Return the list of `F`-rational points on the affine space self,
        where `F` is a given finite field, or the base ring of self.

        EXAMPLES::

            sage: P = ProjectiveSpace(1, GF(3))
            sage: P.rational_points()
            [(0 : 1), (1 : 1), (2 : 1), (1 : 0)]
            sage: P.rational_points(GF(3^2, 'b'))
            [(0 : 1), (2*b : 1), (b + 1 : 1), (b + 2 : 1), (2 : 1), (b : 1), (2*b + 2 : 1), (2*b + 1 : 1), (1 : 1), (1 : 0)]
        """
        if F == None:
            return [ P for P in self ]
        elif not is_FiniteField(F):
            raise TypeError, "Second argument (= %s) must be a finite field."%F
        return [ P for P in self.base_extend(F) ]

class ProjectiveSpace_rational_field(ProjectiveSpace_field):
    def rational_points(self,bound=0):
        r"""
        Returns the projective points `(x_0:\cdots:x_n)` over
        `\QQ` with `|x_i| \leq` bound.

        INPUT:


        -  ``bound`` - integer


        EXAMPLES::

            sage: PP = ProjectiveSpace(0,QQ)
            sage: PP.rational_points(1)
            [(1)]
            sage: PP = ProjectiveSpace(1,QQ)
            sage: PP.rational_points(2)
            [(-2 : 1), (-1 : 1), (0 : 1), (1 : 1), (2 : 1), (-1/2 : 1), (1/2 : 1), (1 : 0)]
            sage: PP = ProjectiveSpace(2,QQ)
            sage: PP.rational_points(2)
            [(-2 : -2 : 1), (-1 : -2 : 1), (0 : -2 : 1), (1 : -2 : 1), (2 : -2 : 1),
            (-2 : -1 : 1), (-1 : -1 : 1), (0 : -1 : 1), (1 : -1 : 1), (2 : -1 : 1),
            (-2 : 0 : 1), (-1 : 0 : 1), (0 : 0 : 1), (1 : 0 : 1), (2 : 0 : 1), (-2 :
            1 : 1), (-1 : 1 : 1), (0 : 1 : 1), (1 : 1 : 1), (2 : 1 : 1), (-2 : 2 :
            1), (-1 : 2 : 1), (0 : 2 : 1), (1 : 2 : 1), (2 : 2 : 1), (-1/2 : -1 :
            1), (1/2 : -1 : 1), (-1 : -1/2 : 1), (-1/2 : -1/2 : 1), (0 : -1/2 : 1),
            (1/2 : -1/2 : 1), (1 : -1/2 : 1), (-1/2 : 0 : 1), (1/2 : 0 : 1), (-1 :
            1/2 : 1), (-1/2 : 1/2 : 1), (0 : 1/2 : 1), (1/2 : 1/2 : 1), (1 : 1/2 :
            1), (-1/2 : 1 : 1), (1/2 : 1 : 1), (-2 : 1 : 0), (-1 : 1 : 0), (0 : 1 :
            0), (1 : 1 : 0), (2 : 1 : 0), (-1/2 : 1 : 0), (1/2 : 1 : 0), (1 : 0 :
            0)]

        .. note::

           The very simple algorithm works as follows: every point
           `(x_0:\cdots:x_n)` in projective space has a unique
           largest index `i` for which `x_i` is not
           zero. The algorithm then iterates downward on this
           index. We normalize by choosing `x_i` positive. Then,
           the points `x_0,\ldots,x_{i-1}` are the points of
           affine `i`-space that are relatively prime to
           `x_i`. We access these by using the Tuples method.

        AUTHORS:

        - Benjamin Antieau (2008-01-12)
        """
        if not bound > 0:
            raise ValueError, \
                  "Argument bound (= %s) must be a positive integer."

        n = self.dimension_relative()


        Q = [ k-bound for k in range(2*bound+1) ]      # the affine coordinates
        R = [ (k+1) for k in range(bound) ]            # the projective coordinate
        S = [ Tuples(Q,(k+1)) for k in range(n) ]
        pts = []

        i=n

        while i > 0:
            P = [ 0 for _ in range(n+1) ]
            for ai in R:
                P[i]=ai
                for tup in S[i-1]:
                    if gcd([ai]+tup)==1:
                        for j in range(i):
                            P[j]=tup[j]
                        pts.append(self(P))
            i-=1

        # now do i=0; this is treated as a special case so that
        # we don't have all points (1:0),(2,0),(3,0),etc.
        P = [ 0 for _ in range(n+1) ]; P[0]=1
        pts.append(self(P))
        return pts
