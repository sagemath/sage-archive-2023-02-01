"""
Spec of a ring
"""

#*******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.rings.all import is_CommutativeRing, ZZ
import scheme
import point

def is_Spec(X):
    """
    EXAMPLES:
        sage: from sage.schemes.generic.spec import is_Spec
        sage: is_Spec(QQ^3)
        False
        sage: X = Spec(QQ); X
        Spectrum of Rational Field
        sage: is_Spec(X)
        True
    """
    return isinstance(X, Spec)

class Spec(scheme.AffineScheme):
    r"""
    The spectrum of a commutative ring, as a scheme.

    \note{Calling \code{Spec(R)} twice produces two distinct
    (but equal) schemes, which is important for gluing to
    construct more general schemes.}

    EXAMPLES:
        sage: Spec(QQ)
        Spectrum of Rational Field
        sage: Spec(PolynomialRing(QQ, 'x'))
        Spectrum of Univariate Polynomial Ring in x over Rational Field
        sage: Spec(PolynomialRing(QQ, 'x', 3))
        Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        sage: X = Spec(PolynomialRing(GF(49,'a'), 3, 'x')); X
        Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Finite Field in a of size 7^2
        sage: loads(X.dumps()) == X
        True
        sage: A = Spec(ZZ); B = Spec(ZZ)
        sage: A is B
        False
        sage: A == B
        True

    A TypeError is raised if the input is not a CommutativeRing.
        sage: Spec(5)
        Traceback (most recent call last):
        ...
        TypeError: R (=5) must be a commutative ring
        sage: Spec(FreeAlgebra(QQ,2, 'x'))
        Traceback (most recent call last):
        ...
        TypeError: R (=Free Algebra on 2 generators (x0, x1) over Rational Field) must be a commutative ring

    EXAMPLES:
        sage: X = Spec(ZZ)
        sage: X
        Spectrum of Integer Ring
        sage: X.base_scheme()
        Spectrum of Integer Ring
        sage: X.base_ring()
        Integer Ring
        sage: X.dimension()
        1
    """
    def __init__(self, R, S=None, check=True):
        if not is_CommutativeRing(R):
            raise TypeError, "R (=%s) must be a commutative ring"%R
        self.__R = R
        if not S is None:
            if not is_CommutativeRing(S):
                raise TypeError, "S (=%s) must be a commutative ring"%S
            try:
                S.hom(R)
            except TypeError:
                raise ValueError, "There must be a natural map S --> R, but S = %s and R = %s"%(S,R)
            self._base_ring = S

    def _cmp_(self, X):
        """
        Anything that is not a Spec is less than X.
        Spec's are compared with self using comparison
        of the underlying rings.

        EXAMPLES:
            sage: Spec(QQ) == Spec(QQ)
            True
            sage: Spec(QQ) == Spec(ZZ)
            False
            sage: Spec(QQ) == 5
            False
            sage: Spec(GF(5)) < Spec(GF(7))
            True
            sage: Spec(GF(7)) < Spec(GF(5))
            False
        """
        return cmp(self.__R, X.coordinate_ring())

    def _repr_(self):
        return "Spectrum of %s"%self.__R

    def __call__(self, x):
        """
        Create a point of this scheme.
        """
        return point.SchemePoint_spec(self, x)

    def coordinate_ring(self):
        """
        Return the underlying ring of this scheme.

        EXAMPLES:
            sage: Spec(QQ).coordinate_ring()
            Rational Field
            sage: Spec(PolynomialRing(QQ,3, 'x')).coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return self.__R

    def dimension(self):
        """
        Return the relative dimension of this scheme over its base.
        """
        if self.base_ring() == ZZ:
            return self.__R.krull_dimension()
        raise NotImplementedError


SpecZ = Spec(ZZ)
