"""
The Spec functor

AUTHORS:

- William Stein

- Peter Bruin (2014): rewrite Spec as a functor

"""

#*******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.rings.integer_ring import ZZ

from sage.schemes.generic.scheme import AffineScheme, is_AffineScheme

from sage.misc.superseded import deprecated_function_alias

is_Spec = deprecated_function_alias(7946, is_AffineScheme)

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
    pass

SpecZ = Spec(ZZ)
