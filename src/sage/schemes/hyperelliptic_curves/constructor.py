"""
Hyperelliptic curve constructor
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.generic.all import ProjectiveSpace

from hyperelliptic_generic import HyperellipticCurve_generic
from hyperelliptic_finite_field import HyperellipticCurve_finite_field
from hyperelliptic_rational_field import HyperellipticCurve_rational_field
from hyperelliptic_padic_field import HyperellipticCurve_padic_field
from hyperelliptic_g2_generic import HyperellipticCurve_g2_generic
from hyperelliptic_g2_finite_field import HyperellipticCurve_g2_finite_field
from hyperelliptic_g2_rational_field import HyperellipticCurve_g2_rational_field
from hyperelliptic_g2_padic_field import HyperellipticCurve_g2_padic_field

from sage.rings.all import is_FiniteField, is_RationalField, is_Polynomial, is_pAdicField

def HyperellipticCurve(f,h=None,names=None,PP=None):
    r"""
    Returns the hyperelliptic curve `y^2 + h y = f`, for
    univariate polynomials `h` and `f`. If `h`
    is not given, then it defaults to 0.

    INPUT:


    -  ``f`` - univariate polynomial

    -  ``h`` - optional univariate polynomial


    EXAMPLES: A curve with and without the h term::

        sage: R.<x> = QQ[]
        sage: HyperellipticCurve(x^5 + x + 1)
        Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1
        sage: HyperellipticCurve(x^19 + x + 1, x-2)
        Hyperelliptic Curve over Rational Field defined by y^2 + (x - 2)*y = x^19 + x + 1

    A curve over a non-prime finite field::

        sage: k.<a> = GF(9); R.<x> = k[]
        sage: HyperellipticCurve(x^3 + x - 1, x+a)
        Hyperelliptic Curve over Finite Field in a of size 3^2 defined by y^2 + (x + a)*y = x^3 + x + 2

    Here's one where we change the names of the vars in the homogeneous
    polynomial::

        sage: k.<a> = GF(9); R.<x> = k[]
        sage: HyperellipticCurve(x^3 + x - 1, x+a, names=['X','Y'])
        Hyperelliptic Curve over Finite Field in a of size 3^2 defined by Y^2 + (X + a)*Y = X^3 + X + 2
    """
    if not is_Polynomial(f):
        raise TypeError, "Arguments f (=%s) and h (= %s) must be polynomials"%(f,h)
    P = f.parent()
    if h is None:
        h = P(0)
        g = (f.degree()-1)%2
    try:
        h = P(h)
    except TypeError:
        raise TypeError, \
              "Arguments f (=%s) and h (= %s) must be polynomials in the same ring"%(f,h)
    df = f.degree()
    dh_2 = 2*h.degree()
    if dh_2 < df:
        g = (df-1)//2
    elif df < dh_2:
        g = (dh_2-1)//2
    else:
        a0 = f.leading_coefficient()
        b0 = h.leading_coefficient()
        A0 = 4*a0 + b0^2
        if A0 != 0:
            g = (df-1)//2
        else:
            if P(2) == 0:
                raise TypeError, "Arguments define a curve with finite singularity."
            f0 = 4*f + h^2
            d0 = f0.degree()
            g = (d0-1)//2
    R = P.base_ring()
    PP = ProjectiveSpace(2, R)
    if names is None:
        names = ["x","y"]
    if is_FiniteField(R):
        if g == 2:
            return HyperellipticCurve_g2_finite_field(PP, f, h, names=names, genus=g)
        else:
            return HyperellipticCurve_finite_field(PP, f, h, names=names, genus=g)
    elif is_RationalField(R):
        if g == 2:
            return HyperellipticCurve_g2_rational_field(PP, f, h, names=names, genus=g)
        else:
            return HyperellipticCurve_rational_field(PP, f, h, names=names, genus=g)
    elif is_pAdicField(R):
        if g == 2:
            return HyperellipticCurve_g2_padic_field(PP, f, h, names=names, genus=g)
        else:
            return HyperellipticCurve_padic_field(PP, f, h, names=names, genus=g)
    else:
        if g == 2:
            return HyperellipticCurve_g2_generic(PP, f, h, names=names, genus=g)
        else:
            return HyperellipticCurve_generic(PP, f, h, names=names, genus=g)
