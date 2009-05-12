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
    EXAMPLES::

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

    .. note::

        Calling ``Spec(R)`` twice produces two distinct (but equal)
        schemes, which is important for gluing to construct more
        general schemes.

    EXAMPLES::

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

    ::

        sage: Spec(5)
        Traceback (most recent call last):
        ...
        TypeError: R (=5) must be a commutative ring
        sage: Spec(FreeAlgebra(QQ,2, 'x'))
        Traceback (most recent call last):
        ...
        TypeError: R (=Free Algebra on 2 generators (x0, x1) over Rational Field) must be a commutative ring

    EXAMPLES::

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
        """
        EXAMPLES::

            sage: Spec(ZZ)
            Spectrum of Integer Ring
        """
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
        Return the result of the comparison of self and X.

        Spec's are compared with self using comparison of the
        underlying rings.  If X is not a Spec, then the result is
        platform-dependent (either self < X or X < self, but never
        self == X).

        EXAMPLES::

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

        TESTS::

            sage: Spec(QQ).__cmp__(Spec(ZZ))
            1
        """
        return cmp(self.__R, X.coordinate_ring())

    def _repr_(self):
        """
        EXAMPLES::

            sage: Spec(PolynomialRing(QQ, 3, 'x'))
            Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field

        TESTS::

            sage: Spec(PolynomialRing(QQ, 3, 'x'))._repr_()
            'Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field'
        """
        return "Spectrum of %s"%self.__R

    def _latex_(self):
        """
        LaTeX representation of this Spec.

        EXAMPLES::

            sage: S = Spec(PolynomialRing(ZZ, 2, 'x'))
            sage: S
            Spectrum of Multivariate Polynomial Ring in x0, x1 over Integer Ring
            sage: S._latex_()
            '\\mathrm{Spec}(\\Bold{Z}[x_{0}, x_{1}])'
        """
        return "\\mathrm{Spec}(%s)" % self.__R._latex_()

    def __call__(self, x):
        """
        Create a point of this Spec.  The argument `x` can be a prime
        ideal of the coordinate ring, or an element (or list of
        elements) of the coordinate ring which generates a prime
        ideal.

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: P = S(3); P
            Point on Spectrum of Integer Ring defined by the Principal ideal (3) of Integer Ring
            sage: type(P)
            <class 'sage.schemes.generic.point.SchemeTopologicalPoint_prime_ideal'>
            sage: S(ZZ.ideal(next_prime(1000000)))
            Point on Spectrum of Integer Ring defined by the Principal ideal (1000003) of Integer Ring

        ::

            sage: R.<x, y, z> = QQ[]
            sage: S = Spec(R)
            sage: P = S(R.ideal(x, y, z)); P
            Point on Spectrum of Multivariate Polynomial Ring in x, y, z over Rational Field defined by the Ideal (x, y, z) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return point.SchemeTopologicalPoint_prime_ideal(self, x)

    def coordinate_ring(self):
        """
        Return the underlying ring of this scheme.

        EXAMPLES::

            sage: Spec(QQ).coordinate_ring()
            Rational Field
            sage: Spec(PolynomialRing(QQ, 3, 'x')).coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return self.__R

    def is_noetherian(self):
        """
        Return True if this scheme is Noetherian.

        EXAMPLES::

            sage: Spec(ZZ).is_noetherian()
            True
        """
        return self.__R.is_noetherian()

    def dimension_absolute(self):
        """
        Return the absolute dimension of this scheme.

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: S.dimension_absolute()
            1
            sage: S.dimension()
            1
        """
        return self.__R.krull_dimension()

    dimension = dimension_absolute

    def dimension_relative(self):
        """
        Return the relative dimension of this scheme over its base.

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: S.dimension_relative()
            0
        """
        return self.__R.krull_dimension() - self.base_ring().krull_dimension()

SpecZ = Spec(ZZ)
