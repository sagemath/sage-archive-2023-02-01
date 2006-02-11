r"""
Projective $n$ space over a ring.


EXAMPLES:
We construct projective space over various rings of various dimensions.

The simplest projective space:
    sage: ProjectiveSpace(0)
    Projective Space of dimension 0 over Integer Ring

A slightly bigger projective space over $\Q$:
    sage: X = ProjectiveSpace(1000, Q); X
    Projective Space of dimension 1000 over Rational Field
    sage: X.dimension()
    1000

We can use ``over'' notation to create projective spaces over various base rings.
    sage: X = ProjectiveSpace(5)/Q; X
    Projective Space of dimension 5 over Rational Field
    sage: X/CC
    Projective Space of dimension 5 over Complex Field with 53 bits of precision

The third argument specifies the printing names of the generators of the homogenous
coordinate ring.  Using objgens() you can obtain both the space and the generators
as ready to use variables.
    sage: P2, (x,y,z) = ProjectiveSpace(2, Q, 'xyz').objgens()
    sage: P2
    Projective Space of dimension 2 over Rational Field
    sage: x.parent()
    Polynomial Ring in x, y, z over Rational Field

For example, we use $x,y,z$ to define the intersection of two lines.
    sage: V = P2.subscheme([x+y+z, x+y-z]); V
    Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
      z + y + x
      -1*z + y + x
    sage: V.dimension()
    0
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import (MPolynomialRing,
                            is_Field,
                            is_Ring,
                            is_CommutativeRing,
                            is_MPolynomialRing,
                            Z)

from sage.misc.all import latex

import homset
import morphism

import algebraic_scheme
import ambient_space
import affine_space

def is_ProjectiveSpace(x):
    return isinstance(x, ProjectiveSpace_ring)

def ProjectiveSpace(n, R=None, names=None):
    r"""
    Return projective space of dimension $n$ over the ring $R$.

    EXAMPLES:
    The dimension and ring can be given in either order.
        sage: ProjectiveSpace(3, QQ)
        Projective Space of dimension 3 over Rational Field
        sage: ProjectiveSpace(QQ, 5)
        Projective Space of dimension 5 over Rational Field
        sage: P = ProjectiveSpace(QQ, 2, names='XYZ'); P
        Projective Space of dimension 2 over Rational Field
        sage: P.coordinate_ring()
        Polynomial Ring in X, Y, Z over Rational Field

    You can use \code{PP} for ProjectiveSpace:
        sage: PP(5)
        Projective Space of dimension 5 over Integer Ring
        sage: PP(5)/GF(17)
        Projective Space of dimension 5 over Finite Field of size 17

    The default base ring is $\Z$.
        sage: ProjectiveSpace(5)
        Projective Space of dimension 5 over Integer Ring

    There is also an projective space associated each polynomial ring.
        sage: R = GF(7)['x,y,z']
        sage: P = ProjectiveSpace(R); P
        Projective Space of dimension 2 over Finite Field of size 7
        sage: P.coordinate_ring()
        Polynomial Ring in x, y, z over Finite Field of size 7
        sage: P.coordinate_ring() is R
        True

    Projective spaces are not cached, i.e., there can be several with the
    same base ring and dimension (to facilitate glueing constructions).
    """

    if is_Ring(n):
        if R is None and is_MPolynomialRing(n):
            A = ProjectiveSpace(n.ngens()-1, n.base_ring())
            A._coordinate_ring = n
            return A
        (n, R) = (R, n)   # swap

    if R is None:
        R = Z  # default is the integers
    if is_Field(R):
        return ProjectiveSpace_field(n, R, names)
    elif is_CommutativeRing(R):
        return ProjectiveSpace_ring(n, R, names)
    else:
        raise TypeError, "R (=%s) must be a commutative ring"%R

class ProjectiveSpace_ring(ambient_space.AmbientSpace):
    """
    Projective space of dimension $n$ over the ring $R$.

    EXAMPLES:
        sage: X = ProjectiveSpace(3, Q, 'xyzw')
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
        Polynomial Ring in x, y, z, w over Rational Field

    Loading and saving:
        sage: loads(X.dumps()) == X
        True
    """
    def __init__(self, n, R=Z, names=None):
        ambient_space.AmbientSpace.__init__(self, n, R)
        self.__names = names


    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme,
        or raise a TypeError.
        """
        return True

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme, if defined.  Otherwise raise
        a ValueError.

        EXAMPLES:
            sage: ProjectiveSpace(3, GF(19^2), 'abcd').coordinate_ring()
            Polynomial Ring in a, b, c, d over Finite Field in a of size 19^2

            sage: ProjectiveSpace(3).coordinate_ring()
            Polynomial Ring in x_0, x_1, x_2, x_3 over Integer Ring

            sage: ProjectiveSpace(2, Q, ['alpha', 'beta', 'gamma']).coordinate_ring()
            Polynomial Ring in alpha, beta, gamma over Rational Field
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            self._coordinate_ring = MPolynomialRing(self.base_ring(), self.dimension()+1, names=self.__names)
            return self._coordinate_ring

    def _point_morphism_class(self, *args, **kwds):
        return morphism.SchemeMorphism_on_points_projective_space(*args, **kwds)

    def __cmp__(self, right):
        if not isinstance(right, ProjectiveSpace_ring):
            return -1
        return cmp([self.dimension(), self.coordinate_ring()],
                   [right.dimension(), right.coordinate_ring()])

    def _latex_(self):
        return "{\\mathbf P}_{%s}^%s"%(latex(self.base_ring()), self.dimension())


    def _constructor(self, *args, **kwds):
        return ProjectiveSpace(*args, **kwds)

    def _homset_class(self, *args, **kwds):
        return homset.SchemeHomset_projective_coordinates_ring(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return morphism.SchemeMorphism_projective_coordinates_ring(*args, **kwds)

    def _repr_(self):
        return "Projective Space of dimension %s over %s"%(self.dimension(), self.base_ring())

    def _repr_generic_point(self, v=None):
        if v is None:
            v = self.gens()
        return '(%s)'%(" : ".join([str(f) for f in v]))

    def _latex_generic_point(self, v=None):
        if v is None:
            v = self.gens()
        return '\\left(%s\\right)'%(" : ".join([str(latex(f)) for f in v]))


    def subscheme(self, X):
        """
        Return the closed subscheme defined by X.

        INPUT:
            X -- a list or tuple of equations

        EXAMPLES:
            sage: A, (x,y,z) = ProjectiveSpace(2, Q).objgens('xyz')
            sage: X = A.subscheme([x*z^2, y^2*z, x*y^2]); X
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*z^2
              y^2*z
              x*y^2
            sage: X.defining_polynomials ()
            (x*z^2, y^2*z, x*y^2)
            sage: I = X.defining_ideal(); I
            Ideal (y^2*z, x*z^2, x*y^2) of Polynomial Ring in x, y, z over Rational Field
            sage: I.groebner_basis()
            [y^2*z, x*z^2, x*y^2]
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
        Return the $i$-th affine patch of this projective space.  This
        is an ambient affine space $\A^n_R$, where $R$ is the base
        ring of self, whose \code{projective\_embedding} map is $1$
        in the $i$th factor.

        INPUT:
            i -- integer between 0 and dimension of self, inclusive.

        OUTPUT:
            an ambient affine space with fixed projective_embedding map.

        EXAMPLES:
            sage: P = ProjectiveSpace(5) / QQ
            sage: A = P.affine_patch(2)
            sage: A
            Affine Space of dimension 5 over Rational Field
            sage: A.projective_embedding()
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x_0, x_1, x_2, x_3, x_4) to
                    (x_0 : x_1 : x_2 : x_3 : x_4 : 1)
            sage: A.projective_embedding(0)
            Scheme morphism:
              From: Affine Space of dimension 5 over Rational Field
              To:   Projective Space of dimension 5 over Rational Field
              Defn: Defined on coordinates by sending (x_0, x_1, x_2, x_3, x_4) to
                    (1 : x_0 : x_1 : x_2 : x_3 : x_4)
        """
        i = int(i)   # implicit type checking
        n = self.dimension()
        if i < 0 or i > n:
            raise ValueError, "Argument i (= %s) must be between 0 and %s."%(i, n)
        try:
            return self.__affine_patches[i]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        AA = affine_space.AffineSpace(self.base_ring(), n)
        phi = AA.projective_embedding(i, self)
        self.__affine_patches[i] = AA
        return AA


class ProjectiveSpace_field(ProjectiveSpace_ring):
    def _homset_class(self, *args, **kwds):
        return homset.SchemeHomset_projective_coordinates_field(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return morphism.SchemeMorphism_projective_coordinates_field(*args, **kwds)


