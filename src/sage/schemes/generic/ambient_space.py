"""
Ambient Spaces
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import Integer, is_CommutativeRing, ZZ

from sage.structure.parent_gens import ParentWithGens

import algebraic_scheme
import scheme

def is_AmbientSpace(x):
    """
    Return True if `x` is an ambient space.

    EXAMPLES::

        sage: from sage.schemes.generic.ambient_space import is_AmbientSpace
        sage: is_AmbientSpace(ProjectiveSpace(3, ZZ))
        True
        sage: is_AmbientSpace(AffineSpace(2, QQ))
        True
        sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
        sage: is_AmbientSpace(P.subscheme([x+y+z]))
        False
    """
    return isinstance(x, AmbientSpace)

class AmbientSpace(scheme.Scheme, ParentWithGens):
    """
    Base class for ambient spaces over a ring.

    INPUT:


    -  ``n`` - dimension

    -  ``R`` - ring
    """
    def __init__(self, n, R=ZZ):
        """
        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
        """
        if not is_CommutativeRing(R):
            raise TypeError, "R (=%s) must be a commutative ring"%R
        n = Integer(n)
        if n < 0:
            raise ValueError, "n (=%s) must be nonnegative"%n
        self.__n = n
        self._base_ring = R

    #######################################################################
    # Derived classes must overload all of the following functions
    #######################################################################
    def __cmp__(self, right):
        """
        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A.__cmp__(ProjectiveSpace(2, QQ))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _constructor(self):
        """
        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._constructor()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _latex_(self):
        """
        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._latex_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_(self):
        """
        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_generic_point(self, coords=None):
        """
        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._repr_generic_point([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _latex_generic_point(self, coords=None):
        """
        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._latex_generic_point([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a TypeError.

        TEST::

            sage: from sage.schemes.generic.ambient_space import AmbientSpace
            sage: A = AmbientSpace(5, ZZ)
            sage: A._check_satisfies_equations([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    #######################################################################
    # End overloads
    #######################################################################


    def base_extend(self, S, check=True):
        """
        Return the base change of self to the ring `S`, via the natural
        map from the base ring of self to `S`.

        A ValueError is raised if there is no natural map between the
        two rings.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: PQ = P.base_extend(QQ); PQ
            Projective Space of dimension 2 over Rational Field
            sage: PQ.base_extend(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: No natural map from the base ring (=Rational Field) to S (=Finite Field of size 5)
        """
        if is_CommutativeRing(S):
            R = self.base_ring()
            if S == R:
                return self
            if check:
                if not S.has_coerce_map_from(R):
                    raise ValueError, "No natural map from the base ring (=%s) to S (=%s)"%(R, S)
            return self._constructor(self.__n, S, self.variable_names())
        else:
            raise NotImplementedError

    def ambient_space(self):
        """
        Return the ambient space of the scheme self, in this case self
        itself.

        EXAMPLES::

            sage: P = ProjectiveSpace(4, ZZ)
            sage: P.ambient_space() is P
            True

            sage: A = AffineSpace(2, GF(3))
            sage: A.ambient_space()
            Affine Space of dimension 2 over Finite Field of size 3
        """
        return self

    def defining_polynomials(self):
        """
        Return the defining polynomials of the scheme self.  Since
        self is an ambient space, this is an empty list.

        EXAMPLES::

            sage: ProjectiveSpace(2, QQ).defining_polynomials()
            ()
            sage: AffineSpace(0, ZZ).defining_polynomials()
            ()
        """
        return ()

    ######################################################################
    # Associated MPolynomial ring generators
    ######################################################################

    def gen(self, n=0):
        """
        Return the `n`-th generator of the coordinate ring of the
        scheme self.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: P.gen(1)
            y
        """
        return self.coordinate_ring().gen(n)

    def gens(self):
        """
        Return the generators of the coordinate ring of the scheme
        self.

        EXAMPLES::

            sage: AffineSpace(0, QQ).gens()
            ()

            sage: P.<x, y, z> = ProjectiveSpace(2, GF(5))
            sage: P.gens()
            (x, y, z)
        """
        return self.coordinate_ring().gens()

    def ngens(self):
        """
        Return the number of generators of the coordinate ring of the
        scheme self.

        EXAMPLES::

            sage: AffineSpace(0, QQ).ngens()
            0

            sage: ProjectiveSpace(50, ZZ).ngens()
            51
        """
        return len(self.gens())

##     def assign_names(self, names=None):
##         """
##         EXAMPLES:
##             sage: A = AffineSpace(2, QQ, 'ab'); A
##             Affine Space of dimension 2 over Rational Field
##             sage: A.coordinate_ring()
##             Polynomial Ring in a, b over Rational Field
##             sage: A._assign_names('xy'); A.coordinate_ring()
##             Polynomial Ring in x, y over Rational Field
##         """
##         self.coordinate_ring()._assign_names(names)

    def dimension_absolute(self):
        """
        Return the absolute dimension of this scheme.

        EXAMPLES::

            sage: A2Q = AffineSpace(2, QQ)
            sage: A2Q.dimension_absolute()
            2
            sage: A2Q.dimension()
            2
            sage: A2Z = AffineSpace(2, ZZ)
            sage: A2Z.dimension_absolute()
            3
            sage: A2Z.dimension()
            3
        """
        base = self.base_scheme()
        if base.is_noetherian():
            return self.dimension_relative() + base.dimension()
        raise NotImplementedError, "Cannot compute the dimension of this scheme."

    dimension = dimension_absolute

    def dimension_relative(self):
        """
        Return the relative dimension of this scheme over its base.

        EXAMPLES::

            sage: A2Q = AffineSpace(2, QQ)
            sage: A2Q.dimension_relative()
            2
            sage: A2Z = AffineSpace(2, ZZ)
            sage: A2Z.dimension_relative()
            2
        """
        return self.__n
