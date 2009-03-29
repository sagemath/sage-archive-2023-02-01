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
    return isinstance(x, AmbientSpace)

class AmbientSpace(scheme.Scheme, ParentWithGens):
    """
    Base class for ambient spaces over a ring.

    INPUT:


    -  ``n`` - dimension

    -  ``R`` - ring
    """
    def __init__(self, n, R=ZZ):
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
        raise NotImplementedError

    def _constructor(self):
        raise NotImplementedError

    def _latex_(self):
        raise NotImplementedError

    def _repr_(self):
        raise NotImplementedError

    def _repr_generic_point(self, coords=None):
        raise NotImplementedError

    def _latex_generic_point(self, coords=None):
        raise NotImplementedError

    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v define a point on this scheme, or
        raise a TypeError.
        """
        return True

    #######################################################################
    # End overloads
    #######################################################################


    def base_extend(self, S, check=True):
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
        return self

    def defining_polynomials(self):
        return ()

    ######################################################################
    # Associated MPolynomial ring generators
    ######################################################################

    def gen(self, n=0):
        return self.coordinate_ring().gen(n)

    def gens(self):
        return self.coordinate_ring().gens()

    def ngens(self):
        raise NotImplementedError

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
