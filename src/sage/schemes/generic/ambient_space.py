"""
Ambient Spaces
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import Integer, is_CommutativeRing, Z
from sage.structure.all import Generators

import algebraic_scheme
import scheme

def is_AmbientSpace(x):
    return isinstance(x, AmbientSpace)

class AmbientSpace(scheme.Scheme, Generators):
    """
    Base class for ambient spaces over a ring.

    INPUT:
        n -- dimension
        R -- ring
    """
    def __init__(self, n, R=Z):
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
        Verify that the coordinates of v define a point on this scheme,
        or raise a TypeError.
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
                try:
                    S._coerce_(R(1))  # make sure there is a natural morphism R --> S
                except TypeError:
                    raise ValueError, "No natural map from the base ring (=%s) to S (=%s)"%(R, S)
            return self._constructor(self.__n, S)
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
        return self.coordinate_ring().ngens()

    def assign_names(self, names=None):
        """
        EXAMPLES:
            sage: A = AffineSpace(2, Q, 'ab'); A
            Affine Space of dimension 2 over Rational Field
            sage: A.coordinate_ring()
            Polynomial Ring in a, b over Rational Field
            sage: A.assign_names('xy'); A.coordinate_ring()
            Polynomial Ring in x, y over Rational Field
        """
        self.coordinate_ring().assign_names(names)

    def dimension(self):
        """
        Return the relative dimension of this scheme over its base.
        """
        return self.__n

