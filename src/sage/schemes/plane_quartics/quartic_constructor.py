"""
Quartic curve constructor
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.projective.projective_space import is_ProjectiveSpace, ProjectiveSpace
from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial

from quartic_generic import QuarticCurve_generic

def QuarticCurve(F, PP=None, check=False):
    """
    Returns the quartic curve defined by the polynomial F.

    INPUT:

    - F -- a polynomial in three variables, homogeneous of degree 4

    - PP -- a projective plane (default:None)

    - check -- whether to check for smoothness or not (default:False)

    EXAMPLES::

        sage: x,y,z=PolynomialRing(QQ,['x','y','z']).gens()
        sage: QuarticCurve(x**4+y**4+z**4)
        Quartic Curve over Rational Field defined by x^4 + y^4 + z^4

    TESTS::

        sage: QuarticCurve(x**3+y**3)
        Traceback (most recent call last):
        ...
        ValueError: Argument F (=x^3 + y^3) must be a homogeneous polynomial of degree 4

        sage: QuarticCurve(x**4+y**4+z**3)
        Traceback (most recent call last):
        ...
        ValueError: Argument F (=x^4 + y^4 + z^3) must be a homogeneous polynomial of degree 4

        sage: x,y=PolynomialRing(QQ,['x','y']).gens()
        sage: QuarticCurve(x**4+y**4)
        Traceback (most recent call last):
        ...
        ValueError: Argument F (=x^4 + y^4) must be a polynomial in 3 variables

    """
    if not is_MPolynomial(F):
        raise ValueError("Argument F (=%s) must be a multivariate polynomial"%F)
    P = F.parent()
    if not P.ngens() == 3:
        raise ValueError("Argument F (=%s) must be a polynomial in 3 variables"%F)
    if not(F.is_homogeneous() and F.degree()==4):
        raise ValueError("Argument F (=%s) must be a homogeneous polynomial of degree 4"%F)

    if not PP is None:
        if not is_ProjectiveSpace(PP) and PP.dimension == 2:
            raise ValueError("Argument PP (=%s) must be a projective plane"%PP)
    else:
        PP = ProjectiveSpace(P)

    if check:
        raise NotImplementedError("Argument checking (for nonsingularity) is not implemented.")

    return QuarticCurve_generic(PP, F)
