"""
Spectrum of a ring as affine scheme.
"""

#*******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.rings.commutative_ring import is_CommutativeRing
from sage.rings.integer_ring import ZZ

from sage.schemes.generic.scheme import AffineScheme
from sage.schemes.generic.point import SchemeTopologicalPoint_prime_ideal


def is_Spec(X):
    """
    Test whether ``X`` is a Spec.

    INPUT:

    - ``X`` -- anything.

    OUTPUT:

    Boolean.

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

class Spec(AffineScheme):
    r"""
    The spectrum of a commutative ring, as a scheme.

    .. note::

        Calling ``Spec(R)`` twice produces two distinct (but equal)
        schemes, which is important for gluing to construct more
        general schemes.

    INPUT:

    - ``R`` -- a commutative ring.

    - ``S`` -- a commutative ring (optional, default:`\ZZ`). The base
      ring.

    EXAMPLES::

        sage: Spec(QQ)
        Spectrum of Rational Field
        sage: Spec(PolynomialRing(QQ, 'x'))
        Spectrum of Univariate Polynomial Ring in x over Rational Field
        sage: Spec(PolynomialRing(QQ, 'x', 3))
        Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        sage: X = Spec(PolynomialRing(GF(49,'a'), 3, 'x')); X
        Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Finite Field in a of size 7^2
        sage: TestSuite(X).run(skip = ["_test_an_element", "_test_elements",
        ...                            "_test_some_elements"])

        sage: A = Spec(ZZ); B = Spec(ZZ)
        sage: A is B
        False
        sage: A == B
        True

    A ``TypeError`` is raised if the input is not a commutative ring::

        sage: Spec(5)
        Traceback (most recent call last):
        ...
        TypeError: R (=5) must be a commutative ring
        sage: Spec(FreeAlgebra(QQ,2, 'x'))
        Traceback (most recent call last):
        ...
        TypeError: R (=Free Algebra on 2 generators (x0, x1) over
        Rational Field) must be a commutative ring

    TESTS::

        sage: X = Spec(ZZ)
        sage: X
        Spectrum of Integer Ring
        sage: X.base_scheme()
        Spectrum of Integer Ring
        sage: X.base_ring()
        Integer Ring
        sage: X.dimension()
        1
        sage: Spec(QQ,QQ).base_scheme()
        Spectrum of Rational Field
        sage: Spec(RDF,QQ).base_scheme()
        Spectrum of Rational Field
    """
    def __init__(self, R, S=None):
        """
        Construct the spectrum of the ring ``R``.

        See :class:`Spec` for details.

        EXAMPLES::

            sage: Spec(ZZ)
            Spectrum of Integer Ring
        """
        if not is_CommutativeRing(R):
            raise TypeError("R (=%s) must be a commutative ring"%R)
        self.__R = R
        if not S is None:
            if not is_CommutativeRing(S):
                raise TypeError("S (=%s) must be a commutative ring"%S)
            try:
                S.hom(R)
            except TypeError:
                raise ValueError("There must be a natural map S --> R, but S = %s and R = %s"%(S,R))
        AffineScheme.__init__(self, S)

    def _cmp_(self, X):
        """
        Compare ``self`` and ``X``.

        Spec's are compared with self using comparison of the
        underlying rings.  If X is not a Spec, then the result is
        platform-dependent (either self < X or X < self, but never
        self == X).

        INPUT:

        - ``X`` -- anything.

        OUTPUT:

        ``+1``, ``0``, or ``-1``.

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

    def __hash__(self):
        """
        Return the hash value.

        OUTPUT:

        A 32/64-bit hash value, depending on architecture.

        TESTS::

            sage: hash(Spec(ZZ))
            -1667718069                 # 32-bit
            -5659298568736299957        # 64-bit

            sage: hash(Spec(QQ['x','y','z']))
            -804171295                  # 32-bit
            -4893002889606114847        # 64-bit
        """
        # R is the only defining data, but we'd like to avoid collisions with it.
        return hash("Spec") ^ hash(self.__R)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        String.

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

        OUTPUT:

        String.

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
        Call syntax for Spec.

        INPUT/OUTPUT:

        The argument ``x`` must be one of the following:

        - a prime ideal of the coordinate ring; the output will
          be the corresponding point of X

        - an element (or list of elements) of the coordinate ring
          which generates a prime ideal; the output will be the
          corresponding point of X

        - a ring or a scheme S; the output will be the set X(S) of
          S-valued points on X

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: P = S(3); P
            Point on Spectrum of Integer Ring defined by the Principal ideal (3) of Integer Ring
            sage: type(P)
            <class 'sage.schemes.generic.point.SchemeTopologicalPoint_prime_ideal'>
            sage: S(ZZ.ideal(next_prime(1000000)))
            Point on Spectrum of Integer Ring defined by the Principal ideal (1000003) of Integer Ring

            sage: R.<x, y, z> = QQ[]
            sage: S = Spec(R)
            sage: P = S(R.ideal(x, y, z)); P
            Point on Spectrum of Multivariate Polynomial Ring
            in x, y, z over Rational Field defined by the Ideal (x, y, z)
            of Multivariate Polynomial Ring in x, y, z over Rational Field

        This indicates the fix of :trac:`12734`::
            sage: S = Spec(ZZ)
            sage: S(ZZ)
            Set of rational points of Spectrum of Integer Ring
            sage: S(S)
            Set of rational points of Spectrum of Integer Ring
        """
        if is_CommutativeRing(x):
            return self.point_homset(x)
        from sage.schemes.all import is_Scheme
        if is_Scheme(x):
            return x.Hom(self)

        return SchemeTopologicalPoint_prime_ideal(self, x)

    def _an_element_(self):
        r"""
        Return an element of the spectrum of the ring.

        OUTPUT:

        A point of the affine scheme ``self``.

        EXAMPLES::

            sage: Spec(QQ).an_element()
            Point on Spectrum of Rational Field defined by the Principal ideal (0) of Rational Field
            sage: Spec(ZZ).an_element()    # random output
            Point on Spectrum of Integer Ring defined by the Principal ideal (811) of Integer Ring
        """
        if self.coordinate_ring() is ZZ:
            from sage.rings.arith import random_prime
            return self(random_prime(1000))
        return self(0)

    def coordinate_ring(self):
        """
        Return the underlying ring of this scheme.

        OUTPUT:

        A commutative ring.

        EXAMPLES::

            sage: Spec(QQ).coordinate_ring()
            Rational Field
            sage: Spec(PolynomialRing(QQ, 3, 'x')).coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return self.__R

    def is_noetherian(self):
        """
        Test whether ``self`` is Noetherian.

        OUTPUT:

        Boolean. Return True if this scheme is Noetherian.

        EXAMPLES::

            sage: Spec(ZZ).is_noetherian()
            True
        """
        return self.__R.is_noetherian()

    def dimension_absolute(self):
        """
        Return the absolute dimension of this scheme.

        OUTPUT:

        Integer.

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

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: S.dimension_relative()
            0
        """
        return self.__R.krull_dimension() - self.base_ring().krull_dimension()

    def base_extend(self, R):
        """
        Extend the base ring/scheme.

        INPUT:

        - ``R`` -- an affine scheme or a commutative ring.

        EXAMPLES::

            sage: Spec_ZZ = Spec(ZZ);  Spec_ZZ
            Spectrum of Integer Ring
            sage: Spec_ZZ.base_extend(QQ)
            Spectrum of Rational Field
        """
        if is_CommutativeRing(R):
            return Spec(self.coordinate_ring().base_extend(R), self.base_ring())
        if not self.base_scheme() == R.base_scheme():
            raise ValueError('The new base scheme must be a scheme over the old base scheme.')
        return Spec(self.coordinate_ring().base_extend(new_base.coordinate_ring()),
                    self.base_ring())

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for more details.

        EXAMPLES::

            sage: Spec(QQ)._point_homset(Spec(QQ), Spec(ZZ))
            Set of rational points of Spectrum of Integer Ring
        """
        from sage.schemes.affine.affine_homset import SchemeHomset_points_spec
        return SchemeHomset_points_spec(*args, **kwds)


SpecZ = Spec(ZZ)
