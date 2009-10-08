r"""
Isogenies

An isogeny `\varphi: E_1\to E_2` between two elliptic curves `E_1` and `E_2` is a morphism of curves that sends the
origin of `E_1` to the origin of `E_2`. Such a morphism is automatically a morphism of group schemes and the kernel is a finite subgroup scheme of `E_1`.
Such a subscheme can either be given by a list of generators, which have to be torsion points, or by a polynomial in the coordinate `x` of the Weierstrass equation of `E_1`.

AUTHORS:

- Daniel Shumow <shumow@gmail.com>: 2009-04-19: initial version
- Chris Wuthrich : 7/09: changes: add check of input, not the full list is needed.
"""

#*****************************************************************************
#       Copyright (C) 2009 Daniel Shumow <shumow@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import deepcopy, copy

from sage.categories import homset

from sage.categories.morphism import Morphism

from sage.structure.sage_object import SageObject
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.all import is_Polynomial
from sage.schemes.elliptic_curves.all import EllipticCurve
from sage.schemes.elliptic_curves.all import is_EllipticCurve

from sage.rings.number_field.number_field_base import is_NumberField

from sage.rings.rational_field import is_RationalField

from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

from sage.sets.set import Set

#
# Private function for parsing input to determine the type of algorithm
#
def isogeny_determine_algorithm(E, kernel, codomain, degree, model):
    r"""
    Helper function that allows the various isogeny functions to infer the algorithm type from the parameters passed in.

    If kernel is a list of points on the EllipticCurve `E`, then we assume the algorithm to use is Velu.

    If kernel is a list of coefficients or a univariate polynomial we try to use the Kohel's algorithms.

    EXAMPLES:

    This helper function will be implicitly called by the following examples::

        sage: R.<x> = GF(5)[]
        sage: E = EllipticCurve(GF(5), [0,0,0,1,0])
        sage: phi = EllipticCurveIsogeny(E, x+3)
        sage: phi2 = EllipticCurveIsogeny(E, [GF(5)(3),GF(5)(1)])
        sage: phi == phi2
        True
        sage: phi3 = EllipticCurveIsogeny(E,  E((2,0)) )
        sage: phi3 == phi2
        True
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_determine_algorithm
        sage: isogeny_determine_algorithm(E, x+3, None, None, None)
        'kohel'
        sage: isogeny_determine_algorithm(E, [3, 1], None, None, None)
        'kohel'
        sage: isogeny_determine_algorithm(E, E((2,0)), None, None, None)
        'velu'

    """

    kernel_is_list = (type(kernel) == list)

    if not kernel_is_list and kernel in E :
        kernel = [kernel]
        kernel_is_list = True

    if (is_Polynomial(kernel) or ( kernel_is_list) and (kernel[0] in E.base_ring()) ):
        algorithm = "kohel"
    elif (kernel_is_list) and (kernel[0] in E): # note that if kernel[0] is on an extension of E this condition will be false
        algorithm = "velu"
    else:
        raise ValueError, "Invalid Parameters to EllipticCurveIsogeny constructor."

    return algorithm


def isogeny_codomain_from_kernel(E, kernel, degree=None):
    r"""
    This function computes the isogeny codomain given a kernel.

    INPUT:

    - ``E`` - The domain elliptic curve.

    - ``kernel`` - Either a list of points in the kernel of the isogeny, or a kernel polynomial (specified as a either a univariate polynomial or a coefficient list.)

    - ``degree`` - an integer, (default:None)  optionally specified degree of the kernel.


    OUTPUT:

    (elliptic curve) the codomain of the separable normalized isogeny from this kernel

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_codomain_from_kernel
        sage: E = EllipticCurve(GF(7), [1,0,1,0,1])
        sage: R.<x> = GF(7)[]
        sage: isogeny_codomain_from_kernel(E, [4,1], degree=3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 6 over Finite Field of size 7
        sage: EllipticCurveIsogeny(E, [4,1]).codomain() == isogeny_codomain_from_kernel(E, [4,1], degree=3)
        True
        sage: isogeny_codomain_from_kernel(E, x^3 + x^2 + 4*x + 3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + 6 over Finite Field of size 7
        sage: isogeny_codomain_from_kernel(E, x^3 + 2*x^2 + 4*x + 3)
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 5*x + 2 over Finite Field of size 7

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: kernel_list = [E((15,10)), E((10,3)),E((6,5))]
        sage: isogeny_codomain_from_kernel(E, kernel_list)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15 over Finite Field of size 19

    """

    algorithm = isogeny_determine_algorithm(E, kernel, None, degree, None);

    if ("velu"==algorithm):
        # if we are using Velu's formula, just instantiate the isogeny
        # and return the codomain
        codomain = EllipticCurveIsogeny(E, kernel).codomain()
    elif ("kohel"==algorithm):
        codomain = compute_codomain_kohel(E, kernel, degree)

    return codomain


def compute_codomain_formula(E, v, w):
    r"""
    Given parameters `v` and `w` (as in velu / kohel / etc formulas)
    computes the codomain curve.

    EXAMPLES:

    This formula is used by every Isogeny Instantiation::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, E((1,2)) )
        sage: phi.codomain()
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 13 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_codomain_formula
        sage: v = phi._EllipticCurveIsogeny__v
        sage: w = phi._EllipticCurveIsogeny__w
        sage: compute_codomain_formula(E, v, w)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 13 over Finite Field of size 19

    """

    a1,a2,a3,a4,a6 = E.ainvs()

    A4 = a4 - 5*v
    A6 = a6 - (a1**2 + 4*a2)*v - 7*w

    return EllipticCurve([a1, a2, a3, A4, A6])


def compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4):
    r"""
    The formula for computing `v` and `w` using Kohel's formulas for isogenies of degree 2.

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, [9,1])
        sage: phi
        Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 19 to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 8 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_even_deg1
        sage: a1,a2,a3,a4,a6 = E.ainvs()
        sage: x0 = -9
        sage: y0 = -(a1*x0 + a3)/2
        sage: compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4)
        (18, 9)

    """

    v = (3*x0**2 + 2*a2*x0 + a4 - a1*y0)
    w = x0*v

    return (v,w)


def compute_vw_kohel_even_deg3(b2,b4,s1,s2,s3):
    r"""
    The formula for computing `v` and `w` using Kohel's formulas for isogenies of degree 3.

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x^3 + 7*x^2 + 15*x + 12)
        sage: phi
        Isogeny of degree 4 from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 19 to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_even_deg3
        sage: (b2,b4) = (E.b2(), E.b4())
        sage: (s1, s2, s3) = (-7, 15, -12)
        sage: compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3)
        (4, 7)

    """

    temp1 = (s1**2 - 2*s2)
    v = 3*temp1 + b2*s1/2 + 3*b4/2
    w = 3*(s1**3 - 3*s1*s2 + 3*s3) + b2*temp1/2 + b4*s1/2

    return (v,w)


def compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,n):
    r"""
    This function computes the `v` and `w` according to Kohel's formulas.

    EXAMPLES:

    This function will be implicitly called by the following example::

        sage: E = EllipticCurve(GF(19), [18,17,16,15,14])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x^3 + 14*x^2 + 3*x + 11)
        sage: phi
        Isogeny of degree 7 from Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 15*x + 14 over Finite Field of size 19 to Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 18*x + 18 over Finite Field of size 19
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_odd
        sage: (b2,b4,b6) = (E.b2(), E.b4(), E.b6())
        sage: (s1,s2,s3) = (-14,3,-11)
        sage: compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,3)
        (7, 1)

    """

    v = 6*(s1**2 - 2*s2) + b2*s1 + n*b4
    w = 10*(s1**3 - 3*s1*s2 + 3*s3) + 2*b2*(s1**2 - 2*s2) + 3*b4*s1 + n*b6

    return (v,w)


def compute_codomain_kohel(E, kernel, degree):
    r"""

    This function computes the codomain from the kernel polynomial as per Kohel's formulas.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_codomain_kohel
        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, [9,1])
        sage: phi.codomain() == isogeny_codomain_from_kernel(E, [9,1])
        True
        sage: compute_codomain_kohel(E, [9,1], 2)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 9*x + 8 over Finite Field of size 19
        sage: R.<x> = GF(19)[]
        sage: E = EllipticCurve(GF(19), [18,17,16,15,14])
        sage: phi = EllipticCurveIsogeny(E, x^3 + 14*x^2 + 3*x + 11)
        sage: phi.codomain() == isogeny_codomain_from_kernel(E, x^3 + 14*x^2 + 3*x + 11)
        True
        sage: compute_codomain_kohel(E, x^3 + 14*x^2 + 3*x + 11, 7)
        Elliptic Curve defined by y^2 + 18*x*y + 16*y = x^3 + 17*x^2 + 18*x + 18 over Finite Field of size 19
        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: phi = EllipticCurveIsogeny(E, x^3 + 7*x^2 + 15*x + 12)
        sage: isogeny_codomain_from_kernel(E, x^3 + 7*x^2 + 15*x + 12) == phi.codomain()
        True
        sage: compute_codomain_kohel(E, x^3 + 7*x^2 + 15*x + 12,4)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 3*x + 15 over Finite Field of size 19

    NOTES:

        This function uses the formulas of Section 2.4 of [K96].

    REFERENCES:

    - [K96] Kohel, "Endomorphism Rings of Elliptic Curves over Finite Fields"

    """

    # First set up the polynomial ring

    base_field = E.base_ring()

    if (is_Polynomial(kernel)):
        psi = kernel
        kernel_list = psi.list()
        poly_ring = psi.parent()
        x = psi.variables()[0]
    elif (list == type(kernel)) and (kernel[0] in base_field):
        kernel_list = kernel
        poly_ring = base_field.polynomial_ring()
        psi = poly_ring(kernel_list)
        x = poly_ring.gen()
    else:
        raise ValueError, "input not of correct type"


    # next determine the even / odd part of the isogeny
    psi_2tor = two_torsion_part(E, poly_ring, psi, degree)

    if (0 != psi_2tor.degree()): # even degree case

        psi_quo = poly_ring(psi/psi_2tor)

        if (0 != psi_quo.degree()):
            raise ArithmeticError, "For basic Kohel's algorithm, if the kernel degree is even then the kernel must be contained in the two torsion."

        n = psi_2tor.degree()

        if (1 == n):

            a1,a2,a3,a4,a6 = E.ainvs()

            x0 = -psi_2tor.constant_coefficient()

            # determine y0
            if (2 == base_field.characteristic()):
                y0 = (x0**3 + a2*x0**2 + a4*x0 + a6).sqrt()
            else:
                y0 = -(a1*x0 + a3)/2

            (v,w) = compute_vw_kohel_even_deg1(x0,y0,a1,a2,a4)

        elif (3 == n):

            b2 = E.b2()
            b4 = E.b4()

            s = psi_2tor.list()
            s1 = -s[n-1]
            s2 = s[n-2]
            s3 = -s[n-3]

            (v,w) = compute_vw_kohel_even_deg3(b2,b4,s1,s2,s3)

    else:

        n = psi.degree()

        b2 = E.b2()
        b4 = E.b4()
        b6 = E.b6()

        s1 = 0; s2 = 0; s3 = 0

        if (1 <= n):
            s1 = -kernel_list[n-1]

        if (2 <= n):
            s2 = kernel_list[n-2]

        if (3 <= n):
            s3 = -kernel_list[n-3]

        # initializing these allows us to calculate E2.
        (v,w) = compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,n);

    return compute_codomain_formula(E, v, w)


def two_torsion_part(E, poly_ring, psi, degree):
    r"""

    Returns the greatest common divisor of ``psi`` and the 2 torsion polynomial of E.

    EXAMPLES:

    Every function that computes the kernel polynomial via Kohel's formulas will call this function::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x + 13)
        sage: isogeny_codomain_from_kernel(E, x + 13) == phi.codomain()
        True
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import two_torsion_part
        sage: two_torsion_part(E, R, x+13, 2)
        x + 13

    """

    b2 = E.b2()
    b4 = E.b4()
    b6 = E.b6()

    x = poly_ring.gens()[0]

    if (None==degree) or (0 == degree % 2):

        psi_2 = 4*x**3 + b2*x**2 + 2*b4*x + b6
        psi_G = poly_ring(psi.gcd(psi_2))

    else:

        psi_G = poly_ring(1)

    return psi_G


class EllipticCurveIsogeny(Morphism):
    r"""
    Class Implementing Isogenies of Elliptic Curves

    This class implements separable, normalized isogenies of elliptic curves.

    Several different algorithms for computing isogenies are available.  These include:

    - Velu's Formulas: Velu's original formulas for computing
      isogenies.  This algorithm is selected by giving as the ``kernel`` parameter
      a list of points which generate a finite subgroup.

    - Kohel's Formulas:
      Kohel's original formulas for computing isogenies.
      This algorithm is selected by giving as the ``kernel`` parameter a polynomial (or a coefficient list (little endian)) which will define the kernel of the isogeny.

    INPUT:

    - ``E``         - an elliptic curve, the domain of the isogeny to initialize.
    - ``kernel``    - a kernel, either a point in ``E``, a list of points in ``E``, a kernel polynomial, or ``None``.
                      If initiating from a domain/codomain, this must be set to None.
    - ``codomain``  - an elliptic curve (default:None).  If ``kernel`` is None,
                      then this must be the codomain of a separable normalized isogeny,
                      furthermore, ``degree`` must be the degree of the isogeny from ``E`` to ``codomain``.
                      If ``kernel`` is not None, then this must be isomorphic to the codomain of the
                      normalized separable isogeny defined by ``kernel``,
                      in this case, the isogeny is post composed with an isomorphism so that this parameter is the codomain.
    - ``degree``    - an integer (default:None).
                      If ``kernel`` is None, then this is the degree of the isogeny from ``E`` to ``codomain``.
                      If ``kernel`` is not None, then this is used to determine whether or not to skip a gcd
                      of the kernel polynomial with the two torsion polynomial of ``E``.
    - ``model``     - a string (default:None).  Only supported variable is "minimal", in which case if
                      ``E`` is a curve over the rationals, then the codomain is set to be the unique global minimum model.

    EXAMPLES:

    A simple example of creating an isogeny of a field of small characteristic::

        sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
        sage: phi = EllipticCurveIsogeny(E, E((0,0)) ); phi
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
        sage: phi.degree() == 2
        True
        sage: phi.kernel_polynomial()
        x
        sage: phi.rational_maps()
        ((x^2 + 1)/x, (x^2*y - y)/x^2)
        sage: phi == loads(dumps(phi))
        True

    A more complicated example of a characteristic 2 field::

        sage: E = EllipticCurve(GF(2^4,'alpha'), [0,0,1,0,1])
        sage: P = E((1,1))
        sage: phi_v = EllipticCurveIsogeny(E, P); phi_v
        Isogeny of degree 3 from Elliptic Curve defined by y^2 + y = x^3 + 1 over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + y = x^3 over Finite Field in alpha of size 2^4
        sage: phi_ker_poly = phi_v.kernel_polynomial()
        sage: phi_ker_poly
        x + 1
        sage: ker_poly_list = phi_ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(E, ker_poly_list)
        sage: phi_k == phi_v
        True
        sage: phi_k.rational_maps()
        ((x^3 + x + 1)/(x^2 + 1), (x^3*y + x^2*y + x*y + x + y)/(x^3 + x^2 + x + 1))
        sage: phi_v.rational_maps()
        ((x^3 + x + 1)/(x^2 + 1), (x^3*y + x^2*y + x*y + x + y)/(x^3 + x^2 + x + 1))
        sage: phi_k.degree() == phi_v.degree()
        True
        sage: phi_k.degree()
        3
        sage: phi_k.is_separable()
        True
        sage: phi_v(E(0))
        (0 : 1 : 0)
        sage: alpha = E.base_field().gen()
        sage: Q = E((0, alpha*(alpha + 1)))
        sage: phi_v(Q)
        (1 : alpha^2 + alpha : 1)
        sage: phi_v(P) == phi_k(P)
        True
        sage: phi_k(P) == phi_v.codomain()(0)
        True

    We can create an isogeny that has kernel equal to the full 2 torsion::

        sage: E = EllipticCurve(GF(3), [0,0,0,1,1])
        sage: ker_list = E.division_polynomial(2).list()
        sage: phi = EllipticCurveIsogeny(E, ker_list); phi
        Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 3 to Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 3
        sage: phi(E(0))
        (0 : 1 : 0)
        sage: phi(E((0,1)))
        (1 : 0 : 1)
        sage: phi(E((0,2)))
        (1 : 0 : 1)
        sage: phi(E((1,0)))
        (0 : 1 : 0)
        sage: phi.degree()
        4

    We can also create trivial isogenies with the trivial kernel::

        sage: E = EllipticCurve(GF(17), [11, 11, 4, 12, 10])
        sage: phi_v = EllipticCurveIsogeny(E, E(0))
        sage: phi_v.degree()
        1
        sage: phi_v.rational_maps()
        (x, y)
        sage: E == phi_v.codomain()
        True
        sage: P = E.random_point()
        sage: phi_v(P) == P
        True

        sage: E = EllipticCurve(GF(31), [23, 1, 22, 7, 18])
        sage: phi_k = EllipticCurveIsogeny(E, [1])
        sage: phi_k
        Isogeny of degree 1 from Elliptic Curve defined by y^2 + 23*x*y + 22*y = x^3 + x^2 + 7*x + 18 over Finite Field of size 31 to Elliptic Curve defined by y^2 + 23*x*y + 22*y = x^3 + x^2 + 7*x + 18 over Finite Field of size 31
        sage: phi_k.degree()
        1
        sage: phi_k.rational_maps()
        (x, y)
        sage: phi_k.codomain() == E
        True
        sage: phi_k.kernel_polynomial()
        1
        sage: P = E.random_point(); P == phi_k(P)
        True

    Velu and Kohel also work in characteristic 0::

        sage: E = EllipticCurve(QQ, [0,0,0,3,4])
        sage: P_list = E.torsion_points()
        sage: phi = EllipticCurveIsogeny(E, P_list)
        sage: phi
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 3*x + 4 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 27*x + 46 over Rational Field
        sage: P = E((0,2))
        sage: phi(P)
        (6 : -10 : 1)
        sage: phi_ker_poly = phi.kernel_polynomial()
        sage: phi_ker_poly
        x + 1
        sage: ker_poly_list = phi_ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(E, ker_poly_list); phi_k
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 3*x + 4 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 27*x + 46 over Rational Field
        sage: phi_k(P) == phi(P)
        True
        sage: phi_k == phi
        True
        sage: phi_k.degree()
        2
        sage: phi_k.is_separable()
        True

    A more complicated example over the rationals (with odd degree isogeny)::

        sage: E = EllipticCurve('11a1')
        sage: P_list = E.torsion_points()
        sage: phi_v = EllipticCurveIsogeny(E, P_list); phi_v
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: P = E((16,-61))
        sage: phi_v(P)
        (0 : 1 : 0)
        sage: ker_poly = phi_v.kernel_polynomial(); ker_poly
        x^2 - 21*x + 80
        sage: ker_poly_list = ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(E, ker_poly_list); phi_k
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: phi_k == phi_v
        True
        sage: phi_v(P) == phi_k(P)
        True
        sage: phi_k.is_separable()
        True

    We can also do this same example over the number field defined by the irreducible two torsion polynomial of `E`::

        sage: E = EllipticCurve('11a1')
        sage: P_list = E.torsion_points()
        sage: K.<alpha> = NumberField(x^3 - 2* x^2 - 40*x - 158)
        sage: EK = E.change_ring(K)
        sage: P_list = [EK(P) for P in P_list]
        sage: phi_v = EllipticCurveIsogeny(EK, P_list); phi_v
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158 to Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-7820)*x + (-263580) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
        sage: P = EK((alpha/2,-1/2))
        sage: phi_v(P)
        (122/121*alpha^2 + 1633/242*alpha - 3920/121 : -1/2 : 1)
        sage: ker_poly = phi_v.kernel_polynomial()
        sage: ker_poly
        x^2 - 21*x + 80
        sage: ker_poly_list = ker_poly.list()
        sage: phi_k = EllipticCurveIsogeny(EK, ker_poly_list)
        sage: phi_k
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158 to Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-7820)*x + (-263580) over Number Field in alpha with defining polynomial x^3 - 2*x^2 - 40*x - 158
        sage: phi_v == phi_k
        True
        sage: phi_k(P) == phi_v(P)
        True
        sage: phi_k == phi_v
        True
        sage: phi_k.degree()
        5
        sage: phi_v.is_separable()
        True

    The following example shows how to specify an isogeny from domain and codomain::

        sage: E = EllipticCurve('11a1')
        sage: R.<x> = QQ[]
        sage: f = x^2 - 21*x + 80
        sage: phi = E.isogeny(f)
        sage: E2 = phi.codomain()
        sage: phi_s = EllipticCurveIsogeny(E, None, E2, 5)
        sage: phi_s
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
        sage: phi_s == phi
        True
        sage: phi_s.rational_maps() == phi.rational_maps()
        True


    """

    ####################
    # member variables
    ####################

    #
    # variables common to all algorithms
    #
    __E1 = None  # domain curve
    __E2 = None  # codomain curve

    __degree = None

    __separable = True # This class only implements separable isogenies (for now.)

    __algorithm = None

    __this_hash = None

    __check = None
    #
    # pre isomorphism
    #
    __intermediate_domain = None
    __pre_isomorphism = None
    __prei_x_coord_ratl_map = None
    __prei_y_coord_ratl_map = None

    #
    # post isomorphism
    #

    __intermediate_codomain = None
    __post_isomorphism = None
    __posti_x_coord_ratl_map = None
    __posti_y_coord_ratl_map = None

    #
    # algebraic structs
    #
    __base_field = None
    __poly_ring = None
    __x_var = None
    __y_var = None

    #
    # Rational Maps
    #
    __rational_maps_initialized = False
    __X_coord_rational_map = None
    __Y_coord_rational_map = None

    #
    # The dual
    #
    __dual = None

    #
    # Kernel Data
    #

    __kernel_list = None  # list of elements in the kernel

    __kernel_polynomial_list = None # polynomial stored as a little endian list of coefficients

    __kernel_polynomial = None # polynomial with roots at x values for x-coordinate of points in the kernel

    __inner_kernel_polynomial = None # the inner kernel polynomial (ignoring preisomorphism)

    __n = None


    #
    # member variables common to Velu's formula
    #

    # we keep track of the 2 torsion and non2torsion separately
    __kernel_2tor = None
    __kernel_non2tor = None

    # variables used in Velu's formula (as well as Kohel's variant)
    __v = None
    __w = None


    #
    # member variables specific to Kohel's algorithm.
    #

    __psi = None # psi polynomial
    __phi = None # phi polynomial
    __omega = None # omega polynomial


    #
    # Python Special Functions
    #

    def __init__(self, E, kernel, codomain=None, degree=None, model=None, check=True):
        r"""
        Constructor for EllipticCurveIsogeny class.

        EXAMPLES::

            sage: E = EllipticCurve(GF(2), [0,0,1,0,1])
            sage: phi = EllipticCurveIsogeny(E, [1,1])
            sage: phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 + y = x^3 + 1 over Finite Field of size 2 to Elliptic Curve defined by y^2 + y = x^3 over Finite Field of size 2

            sage: E = EllipticCurve(GF(31), [0,0,0,1,0])
            sage: P = E((2,17))
            sage: phi = EllipticCurveIsogeny(E, P)
            sage: phi
            Isogeny of degree 8 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 31 to Elliptic Curve defined by y^2 = x^3 + 10*x + 28 over Finite Field of size 31

            sage: E = EllipticCurve('17a1')
            sage: phi = EllipticCurveIsogeny(E, [41/3, -55, -1, -1, 1])
            sage: phi
            Isogeny of degree 9 from Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - x - 14 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 56*x - 10124 over Rational Field

            sage: E = EllipticCurve('37a1')
            sage: triv = EllipticCurveIsogeny(E, E(0))
            sage: triv
            Isogeny of degree 1 from Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: triv.rational_maps()
            (x, y)

            sage: E = EllipticCurve('49a3')
            sage: R.<X> = QQ[]
            sage: EllipticCurveIsogeny(E,X^3-13*X^2-58*X+503,check=False)
            Isogeny of degree 7 from Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 107*x + 552 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 5252*x - 178837 over Rational Field

       """

        if not is_EllipticCurve(E):
            raise ValueError, "E parameter must be an EllipticCurve."

        if type(kernel) != type([1,1]) and kernel in E : # a single point was given, we put it in a list
                                                    # the first condition assures that [1,1] is treated as x+1
            kernel = [kernel]

        # if the kernel is None and the codomain isn't
        # calculate the kernel polynomial
        pre_isom = None
        post_isom = None

        self.__check = check

        if (None == kernel) and (None != codomain):

            if (None == degree):
                raise ValueError, "If specifying isogeny by domain and codomain, degree parameter must be set."

            # save the domain/codomain: not really used now
            old_domain = E
            old_codomain = codomain

            (pre_isom, post_isom, E, codomain, kernel) = compute_sequence_of_maps(E, codomain, degree)

        self.__init_algebraic_structs(E)

        algorithm = isogeny_determine_algorithm(E, kernel, codomain, degree, model);

        self.__algorithm = algorithm

        if ("velu"==algorithm):
            self.__init_from_kernel_list(kernel)
        elif ("kohel"==algorithm):
            self.__init_from_kernel_polynomial(kernel, degree)

        self.__compute_E2()

        self.__setup_post_isomorphism(codomain, model)

        if (None != pre_isom):
            self.set_pre_isomorphism(pre_isom)

        if (None != post_isom):
            self.set_post_isomorphism(post_isom)

        # Inheritance house keeping

        self.__perform_inheritance_housekeeping()

        return


    def __call__(self, P, output_base_ring=None):
        r"""
        Function that implements the call-ability of elliptic curve isogenies.

        EXAMPLES::

            sage: E = EllipticCurve(GF(17), [1, 9, 5, 4, 3])
            sage: phi = EllipticCurveIsogeny(E, [6,13,1])
            sage: phi(E((1,0)))
            (15 : 13 : 1)

            sage: E = EllipticCurve(GF(23), [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi(E((1,5)))
            (2 : 0 : 1)
            sage: phi(E(15,20), output_base_ring=GF(23^2,'alpha'))
            (12 : 1 : 1)

            sage: E = EllipticCurve(QQ, [0,0,0,3,0])
            sage: P = E((1,2))
            sage: phi = EllipticCurveIsogeny(E, [0,1])
            sage: phi(P)
            (4 : -4 : 1)
            sage: phi(-P)
            (4 : 4 : 1)


        """
        E1 = self.__E1
        E_P = P.curve()

        change_output_ring = False

        # check that the parent curve of the input point is this curve
        # or that the point is on the same curve but whose base ring
        # is an extension of this ring
        if (E1 != E_P):
            if (E1.a_invariants() != E_P.a_invariants()) :
                raise ValueError, "P must be on a curve with same Weierstrass model as the domain curve of this isogeny."
            change_output_ring = True


        if(P.is_zero()):
            return self.__E2(0)

        (xP, yP) = P.xy()

        if not self.__E1.is_on_curve(xP,yP):
            raise InputError, "Input point must be on the domain curve of this isogeny."

        # if there is a pre isomorphism, apply it
        if (None != self.__pre_isomorphism):
            temp_xP = self.__prei_x_coord_ratl_map(x=xP, y=yP)
            temp_yP = self.__prei_y_coord_ratl_map(x=xP, y=yP)
            (xP, yP) = (temp_xP, temp_yP)

        if ("velu" == self.__algorithm):
            outP = self.__compute_via_velu_numeric(xP, yP)
        elif ("kohel" == self.__algorithm):
            outP = self.__compute_via_kohel_numeric(xP,yP)

        # the intermediate functions return the point at infinity
        # if the input point is in the kernel
        if (outP == self.__intermediate_codomain(0)):
            return self.__E2(0)

        # if there is a post isomorphism, apply it
        if (None != self.__post_isomorphism):
            tempX = self.__posti_x_coord_ratl_map(x=outP[0], y=outP[1])
            tempY = self.__posti_y_coord_ratl_map(x=outP[0], y=outP[1])
            outP = [tempX, tempY]

        if change_output_ring:
            if (None == output_base_ring):
                output_base_ring = E_P.base_ring()
            outE2 = self.__E2.change_ring(output_base_ring)
        else:
            output_base_ring = self.__E2.base_ring()
            outE2 = self.__E2

        R = output_base_ring

        return outE2.point([R(outP[0]), R(outP[1]), R(1)], check=False)



    def __hash__(self):
        r"""
        Function that implements the hash ability of Isogeny objects.

        This hashes the underlying kernel polynomial so that equal isogeny objects have the same hash value.
        Also, this hashes the base field, and domain and codomain curves as well, so that isogenies with
        the same kernel polynomial (over different base fields / curves) hash to different values.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi_v = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi_k = EllipticCurveIsogeny(E, [0,1])
            sage: phi_k.__hash__() == phi_v.__hash__()
            True
            sage: E_F17 = EllipticCurve(GF(17), [0,0,0,1,1])
            sage: phi_p = EllipticCurveIsogeny(E_F17, E_F17([0,1]))
            sage: phi_p.__hash__() == phi_v.__hash__()
            False

        """

        if (None != self.__this_hash):
            return self.__this_hash

        ker_poly_list = self.__kernel_polynomial_list

        if (None == ker_poly_list):
            ker_poly_list = self.__init_kernel_polynomial()

        this_hash = 0

        for a in ker_poly_list:
            this_hash = this_hash.__xor__(hash(a))

        this_hash = this_hash.__xor__(hash(self.__E1))
        this_hash = this_hash.__xor__(hash(self.__E2))
        this_hash = this_hash.__xor__(hash(self.__base_field))

        self.__this_hash = this_hash

        return self.__this_hash


    def __cmp__(self, other):
        r"""
        Function that implements comparisons between isogeny objects.

        This function works by comparing the underlying kernel objects.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi_v = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi_k = EllipticCurveIsogeny(E, [0,1])
            sage: phi_k == phi_v
            True
            sage: E_F17 = EllipticCurve(GF(17), [0,0,0,1,0])
            sage: phi_p = EllipticCurveIsogeny(E_F17, [0,1])
            sage: phi_p == phi_v
            False


        """
        if (not isinstance(other, EllipticCurveIsogeny)):
            return -1

        if (None == self.__kernel_polynomial):
            self.__init_kernel_polynomial()

        return cmp(self.__kernel_polynomial, other.kernel_polynomial())


    def __neg__(self):
        r"""
        Function to implement unary negation (-) operator on isogenies.
        Returns a copy of this isogeny that has been negated.

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(17)(0))
            sage: phi = EllipticCurveIsogeny(E,  E((-1,0)) )
            sage: negphi = -phi
            sage: phi(E((0,1))) + negphi(E((0,1)))
            (0 : 1 : 0)

            sage: E = EllipticCurve(j=GF(19)(1728))
            sage: R.<x> = GF(19)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: negphi = -phi
            sage: phi(E((3,7))) + negphi(E((3,12))) == phi(2*E((3,7)))
            True
            sage: negphi(E((18,6)))
            (17 : 0 : 1)

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve('17a1')
            sage: R.<x> = QQ[]
            sage: f = x - 11/4
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: negphi = -phi
            sage: phi.rational_maps()[0] == negphi.rational_maps()[0]
            True
            sage: P = E((7,13))
            sage: phi(P) + negphi(P)
            (0 : 1 : 0)

        """
        # save off the kernel lists
        kernel_list = self.__kernel_list
        self.__kernel_list = None

        output = deepcopy(self)

        # reset the kernel lists
        output.__kernel_list = copy(kernel_list)
        self.__kernel_list = kernel_list

        output.switch_sign()
        return output



    #
    # Sage Special Functions
    #

    def _repr_(self):
        r"""
        Special sage specific function that implement the functionality to display the isogeny self as a string.

        EXAMPLES::

            sage: E = EllipticCurve(GF(31), [1,0,1,1,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)) )
            sage: phi._repr_()
            'Isogeny of degree 17 from Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 31 to Elliptic Curve defined by y^2 + x*y + y = x^3 + 14*x + 9 over Finite Field of size 31'

            sage: E = EllipticCurve(QQ, [1,0,0,1,9])
            sage: phi = EllipticCurveIsogeny(E, [2,1])
            sage: phi._repr_()
            'Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y = x^3 + x + 9 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 - 59*x + 165 over Rational Field'

        """
        return 'Isogeny of degree ' + self.__degree.__repr__() + ' from ' + \
                 self.__E1.__repr__() + ' to ' + self.__E2.__repr__()


    def _latex_(self):
        r"""
        Special sage specific function that implements functionality to display an isogeny object as a latex string.

        This function returns a latex string representing the isogeny self as the `x` and `y` coordinate rational functions.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, E(0))
            sage: phi._latex_()
            '\\left( x , y \\right)'

            sage: E = EllipticCurve(GF(17), [0,0,0,1,-1])
            sage: R.<X> = GF(17)[]
            sage: phi = EllipticCurveIsogeny(E, X+11)
            sage: phi._latex_()
            '\\left( \\frac{x^{2} + 11 x + 7}{x + 11} , \\frac{x^{2} y + 5 x y + 12 y}{x^{2} + 5 x + 2} \\right)'


        """
        ratl_maps = self.rational_maps()
        return '\\left( %s , %s \\right)' % (ratl_maps[0]._latex_(), ratl_maps[1]._latex_())


    ###########################
    # Private Common Functions
    ###########################

    # delete the hash value
    def __clear_cached_values(self):
        r"""
        A private function to clear the hash if the codomain has been modified by a pre or post isomorphism.

        EXAMPLES::

            sage: F = GF(7);
            sage: E = EllipticCurve(j=F(0))
            sage: phi = EllipticCurveIsogeny(E, [E((0,-1)), E((0,1))])
            sage: old_hash = hash(phi)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (-1,2,-3,4)))
            sage: hash(phi) == old_hash
            False

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: old_ratl_maps = phi.rational_maps()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (-1,0,0,0)))
            sage: old_ratl_maps == phi.rational_maps()
            False
            sage: old_ratl_maps[1] == -phi.rational_maps()[1]
            True

            sage: F = GF(127); R.<x> = F[]
            sage: E = EllipticCurve(j=F(1728))
            sage: f = x^5 + 43*x^4 + 97*x^3 + 81*x^2 + 42*x + 82
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: old_hash = hash(phi)
            sage: old_ratl_maps = phi.rational_maps()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (-13,13,-13,13)))
            sage: old_hash == hash(phi)
            False
            sage: old_ratl_maps == phi.rational_maps()
            False
            sage: phi._EllipticCurveIsogeny__clear_cached_values()

        """
        self.__this_hash = None
        self.__rational_maps_initialized = False
        self.__X_coord_rational_map = None
        self.__Y_coord_rational_map = None
        self.__dual


    # performs the inheritance house keeping
    def __perform_inheritance_housekeeping(self):
        r"""
        Internal helper function, sets values on the super classes of this class.

        EXAMPLES:

        The following examples will implicitly exercise this function::

            sage: E = EllipticCurve(GF(43), [2,3,5,7,11])
            sage: R.<x> = GF(43)[]; f = x + 42
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi._EllipticCurveIsogeny__perform_inheritance_housekeeping()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E2 = phi.codomain()
            sage: post_isom = WeierstrassIsomorphism(E2, (41, 37, 31, 29))
            sage: phi.set_post_isomorphism(post_isom)
            sage: E1pr = WeierstrassIsomorphism(E, (-1, 2, -3, 4)).codomain().codomain()
            sage: pre_isom = E1pr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(pre_isom)

        """

        # one of the superclasses uses these fields
        self._domain = self.__E1
        self._codomain = self.__E2

        # sets up the parent
        parent = homset.Hom(self.__E1(0).parent(), self.__E2(0).parent())
        Morphism.__init__(self, parent)

        return


    # initializes the base field
    def __init_algebraic_structs(self, E):
        r"""
        An internal function for EllipticCurveIsogeny objects that sets up the member variables necessary for algebra.

        EXAMPLES:

        The following tests inherently exercise this function::

            sage: E = EllipticCurve(j=GF(17)(0))
            sage: phi = EllipticCurveIsogeny(E,  E((-1,0)))
            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)

            sage: F = GF(19); R.<x> = F[]
            sage: E = EllipticCurve(j=GF(19)(0))
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi._EllipticCurveIsogeny__init_algebraic_structs(E)

        """
        self.__E1 = E
        self.__base_field = E.base_ring()

        poly_ring = self.__poly_ring = PolynomialRing(self.__base_field, ['x','y'])

        self.__x_var = poly_ring('x')
        self.__y_var = poly_ring('y')

        self.__intermediate_domain = E

        return


    def __compute_E2(self):
        r"""
        Private function that computes and sets the isogeny codomain.

        EXAMPLES:

        These examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2()

            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2()

        """

        if ("velu" == self.__algorithm):
            E2 = self.__compute_E2_via_velu()
        elif ("kohel" == self.__algorithm):
            E2 = self.__compute_E2_via_kohel()

        self.__E2 = E2
        self.__intermediate_codomain = E2

        return


    # initializes the rational maps fields
    def __initialize_rational_maps(self):
        r"""
        Private function that computes and initializes the rational maps.

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps()
            sage: phi.rational_maps()
            ((x^2 + 1)/x, (x^2*y - y)/x^2)

            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.rational_maps()
            ((x^2 + 1)/x, (x^2*y - y)/x^2)
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps()

        """

        if (self.__rational_maps_initialized):
            return

        if ("velu"==self.__algorithm):
            (X_map, Y_map) = self.__initialize_rational_maps_via_velu()
        if ("kohel"==self.__algorithm):
            (X_map, Y_map) = self.__initialize_rational_maps_via_kohel()

        if (None != self.__prei_x_coord_ratl_map):
            prei_X_map = self.__prei_x_coord_ratl_map
            prei_Y_map = self.__prei_y_coord_ratl_map

            X_map = X_map.subs(x=prei_X_map, y=prei_Y_map)
            Y_map = Y_map.subs(x=prei_X_map, y=prei_Y_map)

        if (None != self.__posti_x_coord_ratl_map):
            X_map = self.__posti_x_coord_ratl_map.subs(x=X_map, y=Y_map)
            Y_map = self.__posti_y_coord_ratl_map.subs(x=X_map, y=Y_map)

        self.__X_coord_rational_map = X_map
        self.__Y_coord_rational_map = Y_map

        self.__rational_maps_initialized = True

        return


    def __init_kernel_polynomial(self):
        r"""
        Private function that initializes the kernel polynomial
        (if the algorithm does not take it as a parameter.)

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.kernel_polynomial()
            x
            sage: phi._EllipticCurveIsogeny__init_kernel_polynomial()
            [0, 1]

        """

        if (None != self.__kernel_polynomial_list):
            return self.__kernel_polynomial_list

        if ("velu" == self.__algorithm):
            ker_poly_list = self.__init_kernel_polynomial_velu()
        else:
            raise InputError, "The kernel polynomial should already be defined!"

        return ker_poly_list


    def __set_pre_isomorphism(self, domain, isomorphism):
        r"""
        Private function to set the pre isomorphism and domain (and keep track of the domain of the isogeny.)

        EXAMPLES::

            sage: E = EllipticCurve(GF(43), [2,3,5,7,11])
            sage: R.<x> = GF(43)[]; f = x + 42
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi._EllipticCurveIsogeny__perform_inheritance_housekeeping()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E1pr = WeierstrassIsomorphism(E, (-1, 2, -3, 4)).codomain().codomain()
            sage: pre_isom = E1pr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(pre_isom)
            sage: phi._EllipticCurveIsogeny__set_pre_isomorphism(E, WeierstrassIsomorphism(E, (-1, 3, -3, 4)))
            sage: E == phi.domain()
            True

        """

        self.__E1 = domain

        # set the isomorphism
        self.__pre_isomorphism = isomorphism

        # calculate the isomorphism as a rational map.

        (u, r, s, t) = isomorphism.tuple()

        x = self.__x_var;
        y = self.__y_var;

        self.__prei_x_coord_ratl_map = (x - r)/u**2
        self.__prei_y_coord_ratl_map = (y - s*(x-r) - t)/u**3

        if (None != self.__kernel_polynomial):
            ker_poly = self.__kernel_polynomial
            ker_poly = ker_poly.subs(x=self.__prei_x_coord_ratl_map)
            kp_lc = ker_poly.univariate_polynomial().leading_coefficient()
            ker_poly = (1/kp_lc)*ker_poly
            self.__kernel_polynomial = ker_poly

        self.__perform_inheritance_housekeeping()

        return;


    def __set_post_isomorphism(self, codomain, isomorphism):
        r"""
        Private function to set the post isomorphism and codomain (and keep track of the codomain of the isogeny.)

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E2 = phi.codomain()
            sage: isom = WeierstrassIsomorphism(E2, (-1,2,-3,4))
            sage: phi.set_post_isomorphism(isom)
            sage: phi._EllipticCurveIsogeny__set_post_isomorphism(E2, WeierstrassIsomorphism(phi.codomain(), (1,-2,3,-4)))
            sage: E2 == phi.codomain()
            True

        """

        # set the codomains
        self.__E2 = codomain

        # set the isomorphism
        self.__post_isomorphism = isomorphism

        # calculate the isomorphism as a rational map.

        (u, r, s, t) = isomorphism.tuple()

        x = self.__x_var;
        y = self.__y_var;

        self.__posti_x_coord_ratl_map = (x - r)/u**2
        self.__posti_y_coord_ratl_map = (y - s*(x-r) - t)/u**3

        self.__perform_inheritance_housekeeping()

        return;


    def __setup_post_isomorphism(self, codomain, model):
        r"""
        Private function to set up the post isomorphism given the codomain.

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: E2 = EllipticCurve(GF(7), [0,0,0,5,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)), E2); phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 7
            sage: E3 = EllipticCurve(GF(7), [0,0,0,6,0])
            sage: phi._EllipticCurveIsogeny__setup_post_isomorphism(E3, None)
            sage: phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7

            sage: R.<x> = QQ[]
            sage: E = EllipticCurve(j=1728)
            sage: f = x^3 - x
            sage: phi = EllipticCurveIsogeny(E, f, model='minimal'); phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 - x over Rational Field to Elliptic Curve defined by y^2 = x^3 - x over Rational Field

            sage: phi = EllipticCurveIsogeny(E, f, model=None)
            sage: phi._EllipticCurveIsogeny__setup_post_isomorphism(None, 'minimal')
            sage: phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 - x over Rational Field to Elliptic Curve defined by y^2 = x^3 - x over Rational Field

        """

        # TODO: add checks to make sure that
        # codomain and model parameters are consistent with the
        # algorithm used.

        post_isom = None
        newE2 = None

        oldE2 = self.__E2

        if (None != model):

            if (None != codomain):
                raise ValueError, "Cannot specify a codomain and model flag simultaneously."

            if ('minimal' == model):

                if (not is_RationalField(oldE2.base_field())):
                    raise ValueError, "specifying minimal for model flag only valid with curves over the rational numbers."

                newE2 = oldE2.minimal_model()
                post_isom = oldE2.isomorphism_to(newE2)

            else:
                raise ValueError, "Unknown value of model flag."

        elif (None != codomain):
            if (not is_EllipticCurve(codomain)):
                raise ValueError,  "Codomain parameter must be an elliptic curve."

            if (not oldE2.is_isomorphic(codomain)):
                raise ValueError, "Codomain parameter must be isomorphic to computed codomain isogeny"

            newE2 = codomain
            post_isom = oldE2.isomorphism_to(newE2)

        if (None != post_isom):
            self.__set_post_isomorphism(newE2, post_isom)

        return


    ###########################
    # Velu's Formula Functions
    ###########################

    #
    # Setup function for Velu's formula
    #

    def __init_from_kernel_list(self, kernel_gens):
        r"""
        Private function that initializes the isogeny from a list of points which generate the kernel (For Velu's formulas.)

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0))); phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__init_from_kernel_list([E(0), E((0,0))])

        """
        if self.__check :
            for P in kernel_gens:
                if not P.has_finite_order():
                    raise ValueError, "The points in the kernel must be of finite order."
        # work around the current implementation of torsion points. When they are done better this could be
        # reduced but it won't speed things up too much.
        kernel_list = Set([self.__E1(0)])
        for P in kernel_gens:
            points_to_add = []
            for j in range(P.order()):
                for Q in kernel_list:
                    points_to_add.append(j*P+Q)
            kernel_list += Set(points_to_add)

        self.__kernel_list = kernel_list.list()
        self.__kernel_2tor = {}
        self.__kernel_non2tor = {}

        self.__degree = len(kernel_list)

        self.__sort_kernel_list()

        return


    #
    # Precompute the values in Velu's Formula.
    #
    def __sort_kernel_list(self):
        r"""
        Private function that sorts the list of points in the kernel (For Velu's formulas.)
        sorts out the 2 torsion points, and puts them in a dictionary.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P); phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__kernel_2tor = {}
            sage: phi._EllipticCurveIsogeny__kernel_non2tor = {}
            sage: phi._EllipticCurveIsogeny__sort_kernel_list()

        """

        a1,a2,a3,a4,a6 = self.__E1.ainvs()

        v = 0
        w = 0

        for Q in self.__kernel_list:

            if Q.is_zero():
                continue

            (xQ,yQ) = Q.xy()

            gxQ = 3*xQ**2 + 2*a2*xQ + a4 - a1*yQ
            gyQ = -2*yQ - a1*xQ - a3

            uQ = gyQ**2

            # sort torsion points:
            if (2*yQ == -a1*xQ - a3): # Q is 2-torsion
                vQ = gxQ
                self.__kernel_2tor[xQ] = (xQ,yQ,gxQ,gyQ,vQ,uQ)
                v = v + vQ
                w = w + (uQ + xQ*vQ)
            elif (not self.__kernel_non2tor.has_key(xQ)): # Q is not a 2-torsion
                vQ = 2*gxQ - a1*gyQ
                self.__kernel_non2tor[xQ] = (xQ,yQ,gxQ,gyQ,vQ,uQ)
                v = v + vQ
                w = w + (uQ + xQ*vQ)

        self.__v = v
        self.__w = w

        return


    #
    # Velu's formula computing the codomain curve
    #
    def __compute_E2_via_velu(self):
        r"""
        Private function that computes the codomain via Velu's formulas.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P);
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2_via_velu()
            Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 7

        """
        v = self.__v
        w = self.__w

        return compute_codomain_formula(self.__E1, v,w)


    def __velu_sum_helper(self, Qvalues, a1, a3, x, y):
        r"""
        Private function for Velu's formulas, helper function to help perform the summation.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P);
            sage: Q = E((0,0)); phi(Q)
            (0 : 0 : 1)
            sage: phi.rational_maps()
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2),
             (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))

            sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)) )
            sage: Qvals = phi._EllipticCurveIsogeny__kernel_2tor[0]
            sage: phi._EllipticCurveIsogeny__velu_sum_helper(Qvals, 0, 0, 5, 5)
            (3, 3)
            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__velu_sum_helper(Qvals, 0, 0, x, y)
            (1/x, y/x^2)

        """
        xQ = Qvalues[0]
        yQ = Qvalues[1]
        gxQ = Qvalues[2]
        gyQ = Qvalues[3]
        vQ = Qvalues[4]
        uQ = Qvalues[5]

        t1 = x - xQ
        inv_t1 = t1**-1
        inv_t1_2 = inv_t1**2
        inv_t1_3 = inv_t1_2*inv_t1

        tX = vQ*inv_t1 + uQ*(inv_t1_2)

        tY0 = uQ*(2*y + a1*x + a3)
        tY1 = vQ*(a1*t1 + y - yQ)
        tY2 = a1*uQ - gxQ*gyQ

        tY =  ( tY0*inv_t1_3 + (tY1 + tY2)*inv_t1_2 )

        return (tX, tY)


    def __compute_via_velu_numeric(self, xP, yP):
        r"""
        Private function that sorts the list of points in the kernel (For Velu's formulas.)
        sorts out the 2 torsion points, and puts them in a diction

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P);
            sage: Q = E((0,0)); phi(Q)
            (0 : 0 : 1)
            sage: Q = E((-1,0)); phi(Q)
            (0 : 0 : 1)
            sage: phi._EllipticCurveIsogeny__compute_via_velu_numeric(0, 0)
            (0, 0)
            sage: phi._EllipticCurveIsogeny__compute_via_velu_numeric(-1, 0)
            (0, 0)

        """
        E2 = self.__E2

        # first check if the point is in the kernel
        if ( self.__kernel_2tor.has_key(xP) or self.__kernel_non2tor.has_key(xP) ) :
            return E2(0)

        outP = self.__compute_via_velu(xP,yP)

        return outP


    def __compute_via_velu(self, xP, yP):
        r"""
        Private function for Velu's formulas, to perform the summation.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P);
            sage: Q = E((0,0)); phi(Q)
            (0 : 0 : 1)
            sage: phi.rational_maps()
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2),
             (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))
            sage: phi._EllipticCurveIsogeny__compute_via_velu(0, 0)
            (0, 0)
            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__compute_via_velu(x, y)
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2),
             (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))

        """

        ker_2tor = self.__kernel_2tor
        ker_non2tor = self.__kernel_non2tor

        X = 0
        Y = 0

        a1 = self.__E1.a1()
        a3 = self.__E1.a3()

        # next iterate over the 2torsion points of the kernel
        for Qvalues in ker_2tor.itervalues():
            (tX, tY) = self.__velu_sum_helper(Qvalues, a1, a3, xP, yP)
            X = X + tX
            Y = Y + tY

        for Qvalues in ker_non2tor.itervalues():
            (tX, tY) = self.__velu_sum_helper(Qvalues, a1, a3, xP, yP)
            X = X + tX
            Y = Y + tY

        X = xP + X
        Y = yP - Y

        return (X,Y)


    def __initialize_rational_maps_via_velu(self):
        r"""
        Private function for Velu's formulas, helper function to initialize the rational maps.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P);
            sage: phi.rational_maps()
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2),
             (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps_via_velu()
            ((x^4 - 2*x^3 + x^2 - 3*x)/(x^3 - 2*x^2 + 3*x - 2),
             (x^5*y - 2*x^3*y - x^2*y - 2*x*y + 2*y)/(x^5 + 3*x^3 + 3*x^2 + x - 1))

        """

        x = self.__x_var
        y = self.__y_var

        return self.__compute_via_velu(x,y)


    def __init_kernel_polynomial_velu(self):
        r"""
        Private function for Velu's formulas, helper function to initialize the rational maps.

        EXAMPLES:

        The following example inherently exercises this function::

            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: P = E((4,2))
            sage: phi = EllipticCurveIsogeny(E, P);
            sage: phi.kernel_polynomial()
            x^2 + 2*x + 4
            sage: phi._EllipticCurveIsogeny__init_kernel_polynomial_velu()
            [4, 2, 1]

        """

        poly_ring = self.__poly_ring
        x = self.__x_var

        invX = 0

        if (None != self.__pre_isomorphism):
            pre_isom = self.__pre_isomorphism
            u = pre_isom.u
            r = pre_isom.r
            invX = (u**2)*x + r
        else:
            invX = x

        psi = poly_ring(1)

        for Qvalues in self.__kernel_2tor.itervalues():
            xQ = invX(x=Qvalues[0])
            psi = psi*(x - xQ)

        for Qvalues in self.__kernel_non2tor.itervalues():
            xQ = invX(x=Qvalues[0])
            psi = psi*(x - xQ)

        ker_poly_list = psi.univariate_polynomial().list()

        self.__kernel_polynomial_list = ker_poly_list
        self.__kernel_polynomial = psi

        return ker_poly_list



    ###################################
    # Kohel's Variant of Velu's Formula
    ###################################

    def __init_from_kernel_polynomial(self, kernel_polynomial, degree=None):
        r"""
        Private function that initializes the isogeny from a kernel polynomial.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,0,0,-1,0])
            sage: phi = EllipticCurveIsogeny(E, x);phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7

            sage: phi._EllipticCurveIsogeny__init_from_kernel_polynomial(x)

            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: phi._EllipticCurveIsogeny__init_from_kernel_polynomial(x+6, degree=3)

        """

        poly_ring = self.__poly_ring
        x = self.__x_var

        E = self.__E1

        if(is_Polynomial(kernel_polynomial)):
            kernel_polynomial = kernel_polynomial.list()

        n = len(kernel_polynomial)-1

        self.__kernel_polynomial_list = kernel_polynomial

        psi = 0
        for j in xrange(len(kernel_polynomial)):
            psi = psi*x + kernel_polynomial[n-j]

        #
        # Determine if kernel polynomial is entirely a two torsion
        #
        psi_G = two_torsion_part(E, poly_ring, psi, degree);

        # force this polynomial to be monic:
        psi_G = psi_G/psi_G.univariate_polynomial().leading_coefficient()

        if (0 != psi_G.degree()): # even degree case

            psi_quo = poly_ring(psi/psi_G)

            if (0 != psi_quo.degree()):
                raise NotImplementedError, "For basic Kohel's algorithm, if the kernel degree is even then the kernel must be contained in the two torsion."

            (phi, omega, v, w, n, d) = self.__init_even_kernel_polynomial(E, psi_G)

        else: # odd degree case

            (phi, omega, v, w, n, d) = self.__init_odd_kernel_polynomial(E, psi)


        #
        # Set up the necessary instance variables
        #

        self.__kernel_polynomial = psi
        self.__inner_kernel_polynomial = psi

        self.__n = n
        self.__degree = d

        self.__psi = psi
        self.__phi = phi
        self.__omega = omega

        self.__v = v
        self.__w = w

        return


    def __init_even_kernel_polynomial(self, E, psi_G):
        r"""
        Private function that initializes the isogeny from a kernel polynomial,
        for Kohel's algorithm in the even degree case.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [-1,0])
            sage: phi = EllipticCurveIsogeny(E, x);phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 7

            sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import two_torsion_part
            sage: psig = two_torsion_part(E,R,x,None)(phi._EllipticCurveIsogeny__x_var)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^3 - x, x^3*y + x*y, 6, 0, 1, 2)

            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: E = EllipticCurve(F, [1,1,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, x); phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + x over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 1 over Finite Field in alpha of size 2^4

            sage: psig = two_torsion_part(E,R,x,None)(phi._EllipticCurveIsogeny__x_var)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^3 + x, x^3*y + x^2 + x*y, 1, 0, 1, 2)

            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: R.<x> = GF(7)[]
            sage: f = x^3 + 6*x^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 2*x + 5 over Finite Field of size 7
            sage: psig = two_torsion_part(E,R,f,None)
            sage: psig = two_torsion_part(E,R,f,None)(phi._EllipticCurveIsogeny__x_var)
            sage: phi._EllipticCurveIsogeny__init_even_kernel_polynomial(E,psig)
            (x^7 - 2*x^6 + 2*x^5 - x^4 + 3*x^3 - 2*x^2 - x + 3,
            x^9*y - 3*x^8*y + 2*x^7*y - 3*x^3*y + 2*x^2*y + x*y - y,
            1,
            6,
            3,
            4)


        """


        #check if the polynomial really divides the two_torsion_polynomial
        if  self.__check and E.division_polynomial(2)(self.__x_var) % psi_G  != 0 :
            raise ValueError, "The polynomial does not define a finite subgroup of the elliptic curve."

        n = psi_G.degree()
        d = n+1

        base_field = self.__base_field
        char = base_field.characteristic()

        x = self.__x_var
        y = self.__y_var

        a1,a2,a3,a4,a6 = E.ainvs()

        b2 = E.b2()
        b4 = E.b4()

        if (1 == n):
            x0 = -psi_G.constant_coefficient()

            # determine y0
            if (2 == char):
                y0 = (x0**3 + a2*x0**2 + a4*x0 + a6).sqrt()
            else:
                y0 = -(a1*x0 + a3)/2

            (v,w) = compute_vw_kohel_even_deg1(x0,y0,a1,a2,a4)

            phi = (x*psi_G + v)*psi_G
            omega = (y*psi_G**2 - v*(a1*psi_G + (y - y0)))*psi_G

        elif (3 == n):
            s = psi_G.univariate_polynomial().list()
            s1 = -s[n-1]
            s2 = s[n-2]
            s3 = -s[n-3]

            psi_G_pr = psi_G.derivative(x)
            psi_G_prpr = psi_G_pr.derivative(x)

            phi = (psi_G_pr**2) + (-2*psi_G_prpr + (4*x - s1))*psi_G

            phi_pr = phi.derivative(x)

            psi_2 = 2*y + a1*x + a3

            omega = (psi_2*(phi_pr*psi_G - phi*psi_G_pr) - (a1*phi + a3*psi_G)*psi_G)/2

            phi = phi*psi_G
            omega = omega*psi_G

            (v,w) = compute_vw_kohel_even_deg3(b2,b4,s1,s2,s3)

        else:
            raise ValueError, "input polynomial must be of degree 1 or 3, not %d" % n

        return (phi, omega, v, w, n, d)


    def __init_odd_kernel_polynomial(self, E, psi):
        r"""
        Private function that initializes the isogeny from a kernel polynomial.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__init_odd_kernel_polynomial(E, x+6)
            (x^3 - 2*x^2 + 3*x + 2, x^3*y - 3*x^2*y + x*y, 2, 6, 1, 3)

            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: alpha = F.gen()
            sage: E = EllipticCurve(F, [1,1,F.gen(),F.gen()^2+1,1])
            sage: f = x + alpha^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + (alpha^2+1)*x + 1 over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + alpha*x + alpha^3 over Finite Field in alpha of size 2^4

            sage: R.<x,y> = F[]
            sage: f = x + alpha^2 + 1
            sage: phi._EllipticCurveIsogeny__init_odd_kernel_polynomial(E, f)
            (x^3 + (alpha^2 + 1)*x + (alpha^3 + alpha^2 + alpha),
             x^3*y + (alpha^2 + 1)*x^2*y + (alpha^2 + alpha + 1)*x^2 + (alpha^2 + 1)*x*y + (alpha^2 + alpha)*x + (alpha)*y + (alpha),
             alpha^2 + alpha + 1,
             alpha^3 + alpha^2 + alpha,
             1,
             3)

        """

        n = psi.degree()
        d = 2*n + 1

        # check if the polynomial really divides the torsion polynomial :
        if self.__check and E.division_polynomial(d)(self.__x_var) % psi != 0:
            raise ValueError, "The polynomial does not define a finite subgroup of the elliptic curve."

        x = self.__x_var

        b2 = E.b2()
        b4 = E.b4()
        b6 = E.b6()

        psi_coeffs = psi.univariate_polynomial().list()

        s1 = 0; s2 = 0; s3 = 0

        if (1 <= n):
            s1 = -psi_coeffs[n-1]

        if (2 <= n):
            s2 = psi_coeffs[n-2]

        if (3 <= n):
            s3 = -psi_coeffs[n-3]

        # initializing these allows us to calculate E2.
        (v,w) = compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,n);

        # initialize the polynomial temporary variables

        psi_pr = psi.derivative(x)
        psi_prpr = psi_pr.derivative(x)

        phi = (4*x**3 + b2*x**2 + 2*b4*x + b6)*(psi_pr**2 - psi_prpr*psi) - \
                (6*x**2 + b2*x + b4)*psi_pr*psi + (d*x - 2*s1)*psi**2

        phi_pr = phi.derivative(x)

        omega = 0
        if (2 != self.__base_field.characteristic()):
            omega = self.__compute_omega_fast(E, psi, psi_pr, phi, phi_pr)
        else:
            omega = self.__compute_omega_general(E, psi, psi_pr, phi, phi_pr)

        return (phi, omega, v, w, n, d)


    #
    # This is the vast omega computation that works when characteristic is not 2
    #
    def __compute_omega_fast(self, E, psi, psi_pr, phi, phi_pr):
        r"""
        Private function that initializes the omega polynomial (from Kohel's formulas)
        in the case that the characteristic of the underlying field is not 2.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 1 over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

            sage: R.<x,y> = GF(7)[]
            sage: psi = phi._EllipticCurveIsogeny__psi
            sage: psi_pr = psi.derivative(x)
            sage: fi = phi._EllipticCurveIsogeny__phi
            sage: fi_pr = fi.derivative(x)
            sage: phi._EllipticCurveIsogeny__compute_omega_fast(E, psi, psi_pr, fi, fi_pr)
            x^3*y - 3*x^2*y + x*y

        """

        a1 = E.a1()
        a3 = E.a3()

        x = self.__x_var; # 'x'
        y = self.__y_var; # 'y'

        psi_2 = 2*y + a1*x + a3

        # note, the formula below is correct
        # the formula in Kohel's thesis has some typos
        # notably the first plus sign should be a minus
        # as it is here below.

        omega = phi_pr*psi*psi_2/2 - phi*psi_pr*psi_2 - \
                (a1*phi + a3*psi**2)*psi/2

        return omega


    def __compute_omega_general(self, E, psi, psi_pr, phi, phi_pr):
        r"""
        Private function that initializes the omega polynomial (from Kohel's formulas)
        in the case of general characteristic of the underlying field.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: F = GF(2^4, 'alpha'); R.<x> = F[]
            sage: alpha = F.gen()
            sage: E = EllipticCurve(F, [1,1,F.gen(),F.gen()^2+1,1])
            sage: f = x + alpha^2 + 1
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + (alpha^2+1)*x + 1 over Finite Field in alpha of size 2^4 to Elliptic Curve defined by y^2 + x*y + alpha*y = x^3 + x^2 + alpha*x + alpha^3 over Finite Field in alpha of size 2^4

            sage: R.<x,y> = F[]
            sage: psi = phi._EllipticCurveIsogeny__psi
            sage: psi_pr = psi.derivative(x)
            sage: fi = phi._EllipticCurveIsogeny__phi
            sage: fi_pr = fi.derivative(x)
            sage: phi._EllipticCurveIsogeny__compute_omega_general(E, psi, psi_pr, fi, fi_pr)
            x^3*y + (alpha^2 + 1)*x^2*y + (alpha^2 + alpha + 1)*x^2 + (alpha^2 + 1)*x*y + (alpha^2 + alpha)*x + (alpha)*y + (alpha)


        """

        a1,a2,a3,a4,a6 = E.ainvs()

        b2 = E.b2()
        b4 = E.b4()

        n = psi.degree()
        d = 2*n+1

        x = self.__x_var
        y = self.__y_var

        psi_2 = 2*y + a1*x + a3

        psi_coeffs = psi.univariate_polynomial().list()

        if (0 < n):
            s1 = -psi_coeffs[n-1]
        else:
            s1 = 0

        psi_prpr = 0
        cur_x_pow = 1

        #
        # Note: we now get the "derivatives" of psi
        # these are not actually the derivatives
        # furthermore, the formulas in Kohel's
        # thesis are wrong, the correct formulas
        # are coded below
        #

        for j  in xrange(0,n-1):
            psi_prpr = psi_prpr + \
                binomial(j+2,2)*psi_coeffs[(j+2)]*cur_x_pow
            cur_x_pow = x*cur_x_pow

        psi_prprpr = 0
        cur_x_pow = 1

        for j in xrange(0,n-2):
            psi_prprpr = psi_prprpr + \
                (3*binomial(j+3,3))*psi_coeffs[(j+3)]*cur_x_pow
            cur_x_pow = x*cur_x_pow


        omega = phi_pr*psi*y - phi*psi_pr*psi_2 + \
                ((a1*x + a3)*(psi_2**2)*(psi_prpr*psi_pr-psi_prprpr*psi) + \
                (a1*psi_2**2 - 3*(a1*x + a3)*(6*x**2 + b2*x + b4))*psi_prpr*psi + \
                (a1*x**3 + 3*a3*x**2 + (2*a2*a3 - a1*a4)*x + (a3*a4 - 2*a1*a6))*psi_pr**2 + \
                (-(3*a1*x**2 + 6*a3*x + (-a1*a4 + 2*a2*a3)) + \
                (a1*x + a3)*(d*x - 2*s1) )*psi_pr*psi + (a1*s1 + a3*n)*psi**2)*psi

        return omega


    def __compute_via_kohel_numeric(self, xP, yP):
        r"""
        Private function that computes a numeric result of this isogeny (via Kohel's formulas.)

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: P = E((0,1)); phi(P)
            (2 : 0 : 1)
            sage: P = E((1,1)); phi(P)
            (0 : 1 : 0)
            sage: phi._EllipticCurveIsogeny__compute_via_kohel_numeric(0, 1)
            (2, 0)
            sage: phi._EllipticCurveIsogeny__compute_via_kohel_numeric(1, 1)
            (0 : 1 : 0)

        """

        E2 = self.__E2

        # first check if this is a kernel point
        # to avoid a divide by 0 error later
        if(0 == self.__inner_kernel_polynomial(x=xP)):
            return E2(0)

        (xP_out, yP_out) = self.__compute_via_kohel(xP,yP)

        # for some dang reason in some cases
        # when the base_field is a number field
        # xP_out and yP_out do not get evaluated to field elements
        # but rather constant polynomials.
        # So in this case, we do some explicit casting to make sure
        # everything comes out right

        if is_NumberField(self.__base_field) and \
                (1 < self.__base_field.degree()) :
            xP_out = self.__poly_ring(xP_out).constant_coefficient()
            yP_out = self.__poly_ring(yP_out).constant_coefficient()

        return (xP_out,yP_out)


    def __compute_via_kohel(self, xP, yP):
        r"""
        Private function that applies Kohel's formulas.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: P = E((0,1)); phi(P)
            (2 : 0 : 1)
            sage: phi.rational_maps()
            ((x^3 - 2*x^2 + 3*x + 2)/(x^2 - 2*x + 1),
             (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))
            sage: phi._EllipticCurveIsogeny__compute_via_kohel(0,1)
            (2, 0)
            sage: R.<x,y> = GF(7)[]
            sage: phi._EllipticCurveIsogeny__compute_via_kohel(x,y)
            ((x^3 - 2*x^2 + 3*x + 2)/(x^2 - 2*x + 1),
             (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))

        """

        x = self.__x_var
        y = self.__y_var

        psi_out = self.__psi(x=xP,y=yP)
        phi_out = self.__phi(x=xP,y=yP)
        omega_out =self.__omega(x=xP, y=yP)

        psi_inv_out = 1/psi_out

        psi_inv_sq_out = psi_inv_out**2

        X_out = phi_out*(psi_inv_sq_out)
        Y_out = omega_out*(psi_inv_sq_out*psi_inv_out)

        return (X_out, Y_out)


    def __initialize_rational_maps_via_kohel(self):
        r"""
        Private function that computes and initializes the rational maps of this isogeny.

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: phi.rational_maps()
            ((x^3 - 2*x^2 + 3*x + 2)/(x^2 - 2*x + 1),
             (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))
            sage: phi._EllipticCurveIsogeny__initialize_rational_maps_via_kohel()
            ((x^3 - 2*x^2 + 3*x + 2)/(x^2 - 2*x + 1),
             (x^3*y - 3*x^2*y + x*y)/(x^3 - 3*x^2 + 3*x - 1))


        """
        x = self.__x_var
        y = self.__y_var

        (X,Y) = self.__compute_via_kohel(x,y)

        return (X,Y)


    #
    # Kohel's formula computing the codomain curve
    #
    def __compute_E2_via_kohel(self):
        r"""
        Private function that computes and initializes the codomain of the isogeny (via Kohel's.)

        EXAMPLES:

        These examples inherently exercise this private function::

            sage: R.<x> = GF(7)[]
            sage: E = EllipticCurve(GF(7), [0,-1,0,0,1])
            sage: phi = EllipticCurveIsogeny(E, x+6, degree=3)
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7
            sage: phi._EllipticCurveIsogeny__compute_E2_via_kohel()
            Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 4*x + 2 over Finite Field of size 7

        """

        v = self.__v
        w = self.__w

        return compute_codomain_formula(self.__E1, v,w)



    #
    # public isogeny methods
    #

    def domain(self):
        r"""
        Returns the domain curve of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E(0,0))
            sage: phi.domain() == E
            True

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: phi = EllipticCurveIsogeny(E, [17, 1])
            sage: phi.domain()
            Elliptic Curve defined by y^2 + x*y = x^3 + x + 2 over Finite Field of size 31

        """
        return self.__E1


    def codomain(self):
        r"""
        Returns the codomain (range) curve of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.codomain()
            Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: phi = EllipticCurveIsogeny(E, [17, 1])
            sage: phi.codomain()
            Elliptic Curve defined by y^2 + x*y = x^3 + 24*x + 6 over Finite Field of size 31

        """
        return self.__E2


    def degree(self):
        r"""
        Returns the degree of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.degree()
            2
            sage: phi = EllipticCurveIsogeny(E, [0,1,0,1])
            sage: phi.degree()
            4

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: phi = EllipticCurveIsogeny(E, [17, 1])
            sage: phi.degree()
            3

        """
        return self.__degree


    def rational_maps(self):
        r"""
        This function returns this isogeny as a pair of rational maps.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,2,0,1,-1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: phi.rational_maps()
            (x, y)

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.rational_maps()
            ((x^2 + 3)/x, (x^2*y - 3*y)/x^2)


        """
        if (not self.__rational_maps_initialized):
            self.__initialize_rational_maps()
        return (self.__X_coord_rational_map, self.__Y_coord_rational_map)


    def is_separable(self):
        r"""
        This function returns a bool indicating whether or not this isogeny is separable.

        This function always returns true.  This class only implements separable isogenies.

        EXAMPLES::

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.is_separable()
            True

            sage: E = EllipticCurve('11a1')
            sage: phi = EllipticCurveIsogeny(E, E.torsion_points())
            sage: phi.is_separable()
            True


        """
        return self.__separable


    def kernel_polynomial(self):
        r"""
        Returns the kernel polynomial of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,2,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.kernel_polynomial()
            x

            sage: E = EllipticCurve('11a1')
            sage: phi = EllipticCurveIsogeny(E, E.torsion_points())
            sage: phi.kernel_polynomial()
            x^2 - 21*x + 80

            sage: E = EllipticCurve(GF(17), [1,-1,1,-1,1])
            sage: phi = EllipticCurveIsogeny(E, [1])
            sage: phi.kernel_polynomial()
            1

            sage: E = EllipticCurve(GF(31), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E, [0,3,0,1])
            sage: phi.kernel_polynomial()
            x^3 + 3*x


        """
        if(None == self.__kernel_polynomial):
            self.__init_kernel_polynomial()

        return self.__kernel_polynomial.univariate_polynomial()


    def set_pre_isomorphism(self, preWI):
        r"""
        Modifies this isogeny object to pre compose with the given Weierstrass isomorphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(31), [1,1,0,1,-1])
            sage: R.<x> = GF(31)[]
            sage: f = x^3 + 9*x^2 + x + 30
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom)
            sage: phi.rational_maps()
            ((-6*x^4 - 3*x^3 + 12*x^2 + 10*x - 1)/(x^3 + x - 12),
             (3*x^7 + x^6*y - 14*x^6 - 3*x^5 + 5*x^4*y + 7*x^4 + 8*x^3*y - 8*x^3 - 5*x^2*y + 5*x^2 - 14*x*y + 14*x - 6*y - 6)/(x^6 + 2*x^4 + 7*x^3 + x^2 + 7*x - 11))
            sage: phi(Epr((0,22)))
            (13 : 21 : 1)
            sage: phi(Epr((3,7)))
            (14 : 17 : 1)

            sage: E = EllipticCurve(GF(29), [0,0,0,1,0])
            sage: R.<x> = GF(29)[]
            sage: f = x^2 + 5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: inv_isom = WeierstrassIsomorphism(E, (1,-2,5,10))
            sage: Epr = inv_isom.codomain().codomain()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom); phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + 10*x*y + 20*y = x^3 + 27*x^2 + 6 over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 20*x over Finite Field of size 29
            sage: phi(Epr((12,1)))
            (26 : 0 : 1)
            sage: phi(Epr((2,9)))
            (0 : 0 : 1)
            sage: phi(Epr((21,12)))
            (3 : 0 : 1)
            sage: phi.rational_maps()[0]
            (x^5 - 10*x^4 - 6*x^3 - 7*x^2 - x + 3)/(x^4 - 8*x^3 + 5*x^2 - 14*x - 6)

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f); phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom)
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 - 13392*x - 1080432 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: phi(Epr((168,1188)))
            (0 : 1 : 0)

        """

        WIdom = preWI.domain().codomain()
        WIcod = preWI.codomain().codomain()

        if (type(preWI) != WeierstrassIsomorphism):
            raise ValueError, "Invalid parameter: isomorphism must be of type Weierstrass isomorphism."

        if (self.__E1 != WIcod):
            raise ValueError, "Invalid parameter: isomorphism must have codomain curve equal to this isogenies' domain."

        if (None == self.__pre_isomorphism):
            isom = preWI
            domain = WIdom
        else:
            isom = self.__pre_isomorphism*preWI
            domain = WIdom

        self.__clear_cached_values()

        self.__set_pre_isomorphism(domain, isom)

        return


    def set_post_isomorphism(self, postWI):
        r"""
        Modifies this isogeny object to post compose with the given Weierstrass isomorphism.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(31)(0))
            sage: R.<x> = GF(31)[]
            sage: phi = EllipticCurveIsogeny(E, x+18)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(), (6,8,10,12)))
            sage: phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 31 to Elliptic Curve defined by y^2 + 24*x*y + 7*y = x^3 + 22*x^2 + 16*x + 20 over Finite Field of size 31

            sage: E = EllipticCurve(j=GF(47)(0))
            sage: f = E.torsion_polynomial(3)/3
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: post_isom = E2.isomorphism_to(E)
            sage: phi.set_post_isomorphism(post_isom)
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            False
            sage: phi.switch_sign()
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            True

        Example over a number field::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2 + 2)
            sage: E = EllipticCurve(j=K(1728))
            sage: ker_list = E.torsion_points()
            sage: phi = EllipticCurveIsogeny(E, ker_list)
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: post_isom = WeierstrassIsomorphism(phi.codomain(), (a,2,3,5))
            sage: phi
            Isogeny of degree 4 from Elliptic Curve defined by y^2 = x^3 + x over Number Field in a with defining polynomial x^2 + 2 to Elliptic Curve defined by y^2 = x^3 + (-44)*x + 112 over Number Field in a with defining polynomial x^2 + 2

        """

        WIdom = postWI.domain().codomain()
        WIcod = postWI.codomain().codomain()

        if (type(postWI) != WeierstrassIsomorphism):
            raise ValueError, "Invalid parameter: isomorphism must be of type Weierstrass isomorphism."

        if (self.__E2 != WIdom):
            raise ValueError, "Invalid parameter: isomorphism must have domain curve equal to this isogenies' codomain."

        if (None == self.__post_isomorphism):
            isom = postWI
            codomain = WIcod
        else:
            isom = postWI*self.__post_isomorphism
            codomain = WIcod

        self.__clear_cached_values()

        self.__set_post_isomorphism(codomain, isom)

        return


    def get_pre_isomorphism(self):
        r"""
        Returns the pre-isomorphism of this isogeny.
        If there has been no pre-isomorphism set, this returns None.

        EXAMPLES::

            sage: E = EllipticCurve(GF(31), [1,1,0,1,-1])
            sage: R.<x> = GF(31)[]
            sage: f = x^3 + 9*x^2 + x + 30
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi.get_post_isomorphism()
            sage: Epr = E.short_weierstrass_model()
            sage: isom = Epr.isomorphism_to(E)
            sage: phi.set_pre_isomorphism(isom)
            sage: isom == phi.get_pre_isomorphism()
            True

            sage: E = EllipticCurve(GF(83), [1,0,1,1,0])
            sage: R.<x> = GF(83)[]; f = x+24
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: phi2 = EllipticCurveIsogeny(E, None, E2, 2)
            sage: phi2.get_pre_isomorphism()
            Generic morphism:
              From: Abelian group of points on Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 83
              To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 62*x + 74 over Finite Field of size 83
              Via:  (u,r,s,t) = (1, 76, 41, 3)



        """
        return self.__pre_isomorphism


    def get_post_isomorphism(self):
        r"""
        Returns the post-isomorphism of this isogeny.
        If there has been no post-isomorphism set, this returns None.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(31)(0))
            sage: R.<x> = GF(31)[]
            sage: phi = EllipticCurveIsogeny(E, x+18)
            sage: phi.get_post_isomorphism()
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (6,8,10,12))
            sage: phi.set_post_isomorphism(isom)
            sage: isom == phi.get_post_isomorphism()
            True

            sage: E = EllipticCurve(GF(83), [1,0,1,1,0])
            sage: R.<x> = GF(83)[]; f = x+24
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: E2 = phi.codomain()
            sage: phi2 = EllipticCurveIsogeny(E, None, E2, 2)
            sage: phi2.get_post_isomorphism()
            Generic morphism:
              From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 65*x + 69 over Finite Field of size 83
              To:   Abelian group of points on Elliptic Curve defined by y^2 + x*y + 77*y = x^3 + 49*x + 28 over Finite Field of size 83
              Via:  (u,r,s,t) = (82, 7, 41, 3)

        """
        return self.__post_isomorphism


    def switch_sign(self):
        r"""
        This function composes the isogeny with `[-1]` (flipping the
        coefficient between +/-1 on the `y` coordinate rational map).

        EXAMPLES::

            sage: E = EllipticCurve(GF(23), [0,0,0,1,0])
            sage: f = E.torsion_polynomial(3)/3
            sage: phi = EllipticCurveIsogeny(E, f, E)
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            False
            sage: phi.switch_sign()
            sage: phi.rational_maps() == E.multiplication_by_m(3)
            True

            sage: E = EllipticCurve(GF(17), [-2, 3, -5, 7, -11])
            sage: R.<x> = GF(17)[]
            sage: f = x+6
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 7*x + 6 over Finite Field of size 17 to Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 4*x + 8 over Finite Field of size 17
            sage: phi.rational_maps()
            ((x^2 + 6*x + 4)/(x + 6), (x^2*y - 5*x*y + 8*x - 2*y)/(x^2 - 5*x + 2))
            sage: phi.switch_sign()
            sage: phi
            Isogeny of degree 2 from Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 7*x + 6 over Finite Field of size 17 to Elliptic Curve defined by y^2 + 15*x*y + 12*y = x^3 + 3*x^2 + 4*x + 8 over Finite Field of size 17
            sage: phi.rational_maps()
            ((x^2 + 6*x + 4)/(x + 6),
             (2*x^3 - x^2*y - 5*x^2 + 5*x*y - 4*x + 2*y + 7)/(x^2 - 5*x + 2))

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: (xmap1, ymap1) = phi.rational_maps()
            sage: phi.switch_sign()
            sage: (xmap2, ymap2) = phi.rational_maps()
            sage: xmap1 == xmap2
            True
            sage: ymap1 == -ymap2 - E.a1()*xmap2 - E.a3()
            True

            sage: K.<a> = NumberField(x^2 + 1)
            sage: E = EllipticCurve(K, [0,0,0,1,0])
            sage: R.<x> = K[]
            sage: phi = EllipticCurveIsogeny(E, x-a)
            sage: phi.rational_maps()
            ((x^2 + (-a)*x - 2)/(x + (-a)), (x^2*y + (-2*a)*x*y + y)/(x^2 + (-2*a)*x - 1))
            sage: phi.switch_sign()
            sage: phi.rational_maps()
            ((x^2 + (-a)*x - 2)/(x + (-a)), (-x^2*y + (2*a)*x*y - y)/(x^2 + (-2*a)*x - 1))

        """
        self.set_post_isomorphism(WeierstrassIsomorphism(self.__E2, (-1,0,-self.__E2.a1(),-self.__E2.a3())))


    def is_normalized(self, check_by_pullback=True):
        r"""
        Returns true if this isogeny is normalized.

        This code does two things.  If the check_by_pullback flag is set, then the code checks that the coefficient on the pullback of the
        invariant differential is 1.  However, in some cases (after a translation has been applied) the underlying polynomial algebra code
        can not sufficiently simplify the pullback expression.  As such, we also cheat a little by falling back and seeing if the post isomorphism
        on this isogeny is a translation with no rescaling.

        INPUT:

        - ``check_by_pullback`` - If this flag is true then the code
          checks the coefficient on the pullback of the invariant
          differential.  (default:True)

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (3, 0, 0, 0))
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            False
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (5, 0, 0, 0))
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1, 1, 1, 1))
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            True

            sage: F = GF(2^5, 'alpha'); alpha = F.gen()
            sage: E = EllipticCurve(F, [1,0,1,1,1])
            sage: R.<x> = F[]
            sage: phi = EllipticCurveIsogeny(E, x+1)
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (alpha, 0, 0, 0))
            sage: phi.is_normalized()
            True
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            False
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1/alpha, 0, 0, 0))
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1, 1, 1, 1))
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            True

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^3 - x^2 - 10*x - 79/4
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (2, 0, 0, 0))
            sage: phi.is_normalized()
            True
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            False
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1/2, 0, 0, 0))
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1, 1, 1, 1))
            sage: phi.set_post_isomorphism(isom)
            sage: phi.is_normalized()
            True

        """

        check_prepost_isomorphism = False

        f_normalized = True

        if (check_by_pullback):

            (Xmap, Ymap) = self.rational_maps()

            E1 = self.__E1
            E2 = self.__E2

            a1 = E1.a1()
            a3 = E1.a3()

            a1pr = E2.a1()
            a3pr = E2.a3()

            x = self.__x_var
            y = self.__y_var

            Xmap_pr = Xmap.derivative(x)

            domain_inv_diff = 1/(2*y + a1*x + a3)
            codomain_inv_diff = Xmap_pr/(2*Ymap + a1pr*Xmap + a3pr)

            inv_diff_quo = domain_inv_diff/codomain_inv_diff

            if (1 == inv_diff_quo):
                f_normalized
            else:
                # For some reason, in certain cases, when the isogeny is pre or post composed with a translation
                # the resulting rational functions are too complicated for sage to simplify down to a constant
                # in this case, we do some cheating by checking if the post-composition by isogeny has
                # a non 1 scaling factor
                if ( inv_diff_quo.numerator().is_constant() and (inv_diff_quo.denominator().is_constant) ):
                    f_normalized = False
                else:
                    check_prepost_isomorphism = True
        else:
            check_prepost_isomorphism = True

        #
        # If we skip checking by the pullback of the invariant differential OR if that was inconclusive
        # We explicitly check if there is a post isomorphism and if it has a non 1 scaling factor
        # or if it is a just a translation.
        # NOTE: This only works because we are using algorithms for calculating the isogenies that calculate
        # a separable normalized isogeny, if this changes, this check will no longer be correct.
        #
        if (check_prepost_isomorphism):
            post_isom = self.__post_isomorphism
            if (None != post_isom):
                if (1 == self.__base_field(post_isom.u)):
                    f_post_normalized = True
                else:
                    f_post_normalized = False
            else:
                f_post_normalized = True

            pre_isom = self.__pre_isomorphism
            if (None != pre_isom):
                if (1 == self.__base_field(pre_isom.u)):
                    f_pre_normalized = True
                else:
                    f_pre_normalized = False
            else:
                f_pre_normalized = True

            f_normalized = f_pre_normalized and f_post_normalized

        return f_normalized


    def dual(self):
        r"""
        Computes and returns the dual isogeny of this isogeny.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 - 21*x + 80
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi_hat = phi.dual()
            sage: phi_hat.domain() == phi.codomain()
            True
            sage: phi_hat.codomain() == phi.domain()
            True
            sage: (X, Y) = phi.rational_maps()
            sage: (Xhat, Yhat) = phi_hat.rational_maps()
            sage: Xm = Xhat.subs(x=X, y=Y)
            sage: Ym = Yhat.subs(x=X, y=Y)
            sage: (Xm, Ym) == E.multiplication_by_m(5)
            True

            sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
            sage: R.<x> = GF(37)[]
            sage: f = x^3 + x^2 + 28*x + 33
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi_hat = phi.dual()
            sage: phi_hat.codomain() == phi.domain()
            True
            sage: phi_hat.domain() == phi.codomain()
            True
            sage: (X, Y) = phi.rational_maps()
            sage: (Xhat, Yhat) = phi_hat.rational_maps()
            sage: Xm = Xhat.subs(x=X, y=Y)
            sage: Ym = Yhat.subs(x=X, y=Y)
            sage: (Xm, Ym) == E.multiplication_by_m(7)
            False
            sage: (Xm, -Ym) == E.multiplication_by_m(7)
            True

            sage: E = EllipticCurve(GF(31), [0,0,0,1,8])
            sage: R.<x> = GF(31)[]
            sage: f = x^2 + 17*x + 29
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi_hat = phi.dual()
            sage: phi_hat.codomain() == phi.domain()
            True
            sage: phi_hat.domain() == phi.codomain()
            True
            sage: (X, Y) = phi.rational_maps()
            sage: (Xhat, Yhat) = phi_hat.rational_maps()
            sage: Xm = Xhat.subs(x=X, y=Y)
            sage: Ym = Yhat.subs(x=X, y=Y)
            sage: (Xm, Ym) == E.multiplication_by_m(5)
            False
            sage: (Xm, -Ym) == E.multiplication_by_m(5)
            True

        """

        if (self.__base_field.characteristic() in [2,3]):
            raise NotImplemented

        if (None != self.__dual):
            return self.__dual

        E1 = self.__intermediate_codomain
        E2pr = self.__intermediate_domain

        F = self.__base_field

        d = self.__degree

        isom = WeierstrassIsomorphism(E2pr, (1/F(d), 0, 0, 0))

        E2 = isom.codomain().codomain()

        pre_isom = self.__E2.isomorphism_to(E1)
        post_isom = E2.isomorphism_to(self.__E1)

        phi_hat = EllipticCurveIsogeny(E1, None, E2, d)

        phi_hat.set_pre_isomorphism(pre_isom)
        phi_hat.set_post_isomorphism(post_isom)

        self.__dual = phi_hat

        return phi_hat


    #
    # Overload Morphism methods that we want to
    #

    def _composition_(self, right, homset):
        r"""
        Composition operator function inherited from morphism class.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [E(0), E((0,1)), E((0,-1))])
            sage: phi*phi
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: phi._composition_(phi, phi.parent())
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError


    def is_injective(self):
        r"""
        Method inherited from the morphism class.
        Returns True if and only if this isogeny has trivial kernel.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 + x - 29/5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi.is_injective()
            False
            sage: phi = EllipticCurveIsogeny(E, R(1))
            sage: phi.is_injective()
            True

            sage: F = GF(7)
            sage: E = EllipticCurve(j=F(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,-1)), E((0,1))])
            sage: phi.is_injective()
            False
            sage: phi = EllipticCurveIsogeny(E, E(0))
            sage: phi.is_injective()
            True

        """

        if (1 < self.__degree): return False
        return True


    def is_surjective(self):
        r"""
        For elliptic curve isogenies, always returns True (as a non-constant map of algebraic curves must be surjective).

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 + x - 29/5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi.is_surjective()
            True

            sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.is_surjective()
            True

            sage: F = GF(2^5, 'omega')
            sage: E = EllipticCurve(j=F(0))
            sage: R.<x> = F[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.is_surjective()
            True

        """
        return True

    def is_zero(self):
        r"""
        Member function inherited from morphism class.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,1)), E((0,-1))])
            sage: phi.is_zero()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def post_compose(self, left):
        r"""
        Member function inherited from morphism class.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,1)), E((0,-1))])
            sage: phi.post_compose(phi)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError


    def pre_compose(self, right):
        r"""
        Member function inherited from morphism class.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,1)), E((0,-1))])
            sage: phi.pre_compose(phi)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError


    def n(self):
        r"""
        Numerical Approximation inherited from Map (through morphism), nonsensical for isogenies.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,1)), E((0,-1))])
            sage: phi.n()
            Traceback (most recent call last):
            ...
            NotImplementedError: Numerical approximations do not make sense for Elliptic Curve Isogenies

        """
        raise NotImplementedError, "Numerical approximations do not make sense for Elliptic Curve Isogenies"


def truncated_reciprocal_quadratic(f, n):
    r"""
    Computes the truncated reciprocal of `f`, to precision `n`.
    This algorithm has complexity `O(dM(n))`, where `M(n)` is the cost of a multiplication and `d` is the degree of `f`.

    INPUT:

    - ``f`` - polynomial, to compute the truncated reciprocal of
    - ``n`` - integer, precision to compute reciprocal to

    OUTPUT:

    polynomial -- a polynomial `g`, such that `gf \equiv 1 \pmod {z^n}`

    ALGORITHM:

    This function uses the algorithm described in section 2.1 of [BMSS]

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, Schost, "Fast Algorithms for Isogenies."

    EXAMPLES::

        sage: f = 1 + x + 7*x^2 + 9*x^5
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import truncated_reciprocal_quadratic
        sage: R.<x> = QQ[]
        sage: f = 1 + x + 7*x^2 + 9*x^5
        sage: truncated_reciprocal_quadratic(f, 5)
        29*x^4 + 13*x^3 - 6*x^2 - x + 1
        sage: (f*truncated_reciprocal_quadratic(f, 5)).quo_rem(x**5)[1]
        1

    """

    R = f.parent()
    d = f.degree()

    f_coef = f.coeffs()
    for j in xrange(len(f_coef), n):
        f_coef.append(0)

    g_coef = [0 for j in xrange(n)]

    g_coef[0] = g0 = 1/f.constant_coefficient()

    for j in xrange(1,n):
        gj = 0
        for k in xrange(1,min(j+1,d+1)):
            gj += f_coef[k]*g_coef[j-k]
        gj = -g0*gj
        g_coef[j] = gj

    g = R(g_coef)

    return g


def truncated_reciprocal_newton(f, n):
    r"""
    Computes the truncated reciprocal of `f`, to precision `n` by newton iteration.
    This algorithm has complexity `O(M(n))`, where `M(n)` is the cost of a multiplication.

    INPUT:

    - ``f`` - polynomial, to compute the truncated reciprocal of
    - ``n`` - integer, precision to compute reciprocal to

    OUTPUT:

    polynomial -- polynomial `g`, such that `gf \equiv 1 \pmod {z^n}`

    ALGORITHM:

    This function uses the algorithm described in section 2.1 of [BMSS]

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, SCHOST, "Fast Algorithms for Isogenies."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import truncated_reciprocal_newton
        sage: R.<x> = GF(53)[]
        sage: f = 4 + 3*x^3 + 47*x^7
        sage: truncated_reciprocal_newton(f, 6)
        23*x^3 + 40
        sage: (f*truncated_reciprocal_newton(f, 6)).quo_rem(x**6)[1]
        1

    """

    hj = 1/f.constant_coefficient()

    if (f.is_constant()):
        return hj

    z = f.variables()[0]

    loop_condition = 2

    while (loop_condition < 2*n):
        hj_next = hj*(2 - f*hj)
        (hj_quo, hj) = hj_next.quo_rem(z**loop_condition)
        loop_condition *= 2

    (hj_quo, g) = hj.quo_rem(z**n)

    return g


def truncated_reciprocal(f, n, algorithm="newtoniteration"):
    r"""
    Computes the truncated reciprocal of `f`, to precision `n`.

    INPUT:

    - ``f``         - polynomial, to compute the truncated reciprocal of
    - ``n``         - integer, precision to compute reciprocal to
    - ``algorithm`` - string (default:"newtoniteration"), string that selects the algorithm to use.  Choices are "newtoniteration" or "quadratic"

    OUTPUT:

    polynomial -- a polynomial `g`, such that `gf \equiv 1 \pmod {z^n}`

    ALGORITHM:

    "newtoniteration" algorithm has complexity `O(M(n))` and "quadratic" has algorithm complexity `O(dM(n))`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import truncated_reciprocal
        sage: R.<x> = GF(101)[]
        sage: f = 4 + x - x^2 + 50*x^3 + 91*x^5
        sage: truncated_reciprocal(f, 5)
        37*x^4 + 43*x^3 + 49*x^2 + 82*x + 76
        sage: (f*truncated_reciprocal(f, 5)).quo_rem(x**5)[1]
        1

    """

    if ("quadratic"==algorithm):
        trunc_recip = truncated_reciprocal_quadratic(f, n)
    elif ("newtoniteration"==algorithm):
        trunc_recip = truncated_reciprocal_newton(f, n)
    else:
        raise ValueError, "Unknown algorithm for computing reciprocal"

    return trunc_recip


def truncated_log(f, n):
    r"""
    Computes the truncated logarithm of polynomial `f` to precision `n`.
    The complexity of this function is `O(M(n))`.

    INPUT:

    - ``f``         - polynomial, to compute the truncated logarithm of;  `f(0)` must not equal 0.
    - ``n``         - integer, precision to compute logarithm to

    OUTPUT:

    polynomial -- a polynomial `g`, such that `\exp_n (g) \equiv f \pmod {z^n}`

    ALGORITHM:

    Uses the algorithm from section 2.2 of [BMSS].

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, Schost, "Fast Algorithms for Isogenies."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import truncated_log, truncated_exp
        sage: R.<x> = GF(31)[]
        sage: f = 1 + x + 2*x^2 + 3*x^3
        sage: truncated_log(f, 4)
        3*x^3 + 2*x^2 + x
        sage: g = truncated_log(f, 4)
        sage: truncated_exp(g, 4)
        3*x^3 + 2*x^2 + x + 1

    """

    R = f.parent()
    K = R.base_ring()

    z = R.gen()

    z_mod = z**n

    g = 1 - f
    cur_pow = g

    sum_acc = 0

    for j in xrange(1,n):
        sum_acc -= K(1/j)*cur_pow
        cur_pow *= g
        cur_pow = cur_pow.quo_rem(z_mod)[1]

    return sum_acc


def truncated_exp_quadratic(f, n):
    r"""
    Computes the truncated exponential of `f`, to precision `n`, using a straight forward power series approximation.

    This algorithm has complexity `O(nM(n))`, where `M(n)` is the
    cost of a multiplication.

    INPUT:

    - ``f`` - polynomial, to compute the truncated exponential of
    - ``n`` - integer, precision to compute reciprocal to

    OUTPUT:

    polynomial -- a polynomial `g`, such that `\log_n(g) \equiv f \pmod {z^n}`.

    ALGORITHM:

    This function uses the algorithm described in section 2.2 of [BMSS]

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, Schost, "Fast Algorithms for Isogenies."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import truncated_exp_quadratic
        sage: R.<x> = GF(127)[]
        sage: f = 1 + x^10
        sage: truncated_exp_quadratic(f, 10)
        31
        sage: truncated_exp_quadratic(f, 11)
        31*x^10 + 123

    """

    R = f.parent()
    K = R.base_ring()

    z = R.gen()

    z_mod = z**n

    g = R(1)
    cur_coef = K(1)

    cur_pow = R(1)

    for j in xrange(1,n):
        cur_coef = cur_coef/K(j)
        cur_pow = cur_pow * f
        cur_pow = cur_pow.quo_rem(z_mod)[1]
        g += cur_coef*cur_pow

    return g


def truncated_exp_fast(f, n):
    r"""
    Computes the truncated exponential of `f`, to precision `n`, using an efficient newton iteration.
    This algorithm has complexity `O(M(n))`, where `M(n)` is the cost of a multiplication.

    INPUT:

    - ``f`` - polynomial, to compute the truncated exponential of
    - ``n`` - integer, precision to compute reciprocal to

    OUTPUT:

    polynomial -- a polynomial `g`, such that `\log_n(g) \equiv f \pmod {z^n}`.

    ALGORITHM:

    This function uses the algorithm described in section 2.2 of [BMSS]

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, Schost, "Fast Algorithms for Isogenies."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import truncated_exp_fast

        sage: R.<x> = GF(27, 'a')[]
        sage: f = 1 + x^2 + 2*x^3 + x^4
        sage: truncated_exp_fast(f, 4)
        x^3 + 2
        sage: truncated_exp_fast(f, 3)
        2*x^2 + 2

    """

    R = f.parent()
    K = R.base_ring()

    z = R.gen()

    gi = R(1)

    j = 0
    m = 0

    while (m < n):

        m = 2**j + 1

        if (n < m):
            m = n

        g_nexti = gi*(1 + f - truncated_log(gi, m))
        gi = g_nexti.quo_rem(z**m)[1]

        j += 1

    return gi


def truncated_exp(f, n, algorithm="fast"):
    r"""
    Computes the truncated exponential of `f`, to precision `n`, using an efficient newton iteration.
    This algorithm has complexity `O(M(n))`, where `M(n)` is the cost of a multiplication.

    INPUT:

    - ``f``         - polynomial, to compute the truncated exponential of
    - ``n``         - integer, precision to compute reciprocal to
    - ``algorithm`` - string (default:"fast") indicates which algorithm to use, choices are "fast" or "quadratic"

    OUTPUT:

    polynomial -- a polynomial `g`, such that `\log_n(g) \equiv f \pmod {z^n}`.

    ALGORITHM:

    algorithm "fast" uses newton iteration and has complexity O(M(n)), algorithm "quadratic" has complexity O(n*M(n))

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import truncated_exp

        sage: R.<x> = GF(29)[]
        sage: f = 1+x+3*x^3
        sage: truncated_exp(f, 3)
        x + 2

    """

    if ("quadratic"==algorithm):
        trunc_exp = truncated_exp_quadratic(f, n)
    elif ("fast"==algorithm):
        trunc_exp = truncated_exp_fast(f, n)
    else:
        raise ValueError, "Unknown algorithm for computing truncated exponential."

    return trunc_exp


def solve_linear_differential_system(a, b, c, alpha, z, n):
    r"""
    Solves a system of linear differential equations:
    `af' + bf = c`, `f'(0) = \alpha`
    where `a`, `b`, `c`, `f`, `f'` are polynomials in variable `z`, and `f`, `f'` are computed to precision `n`.

    EXAMPLES:

    The following examples inherently exercises this function::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_kernel_polynomial, solve_linear_differential_system

        sage: R.<x> = GF(47)[]
        sage: E = EllipticCurve(GF(47), [0,0,0,1,0])
        sage: E2 = EllipticCurve(GF(47), [0,0,0,16,0])
        sage: compute_isogeny_kernel_polynomial(E, E2, 4, algorithm="starks")
        z^3 + z

    """

    z_mod = z**(n-1)

    a_recip = truncated_reciprocal(a, n-1)

    B = (b*a_recip).quo_rem(z_mod)[1]

    C = (c*a_recip).quo_rem(z_mod)[1]

    int_B = B.integral()

    J = truncated_exp(int_B, n)

    J_recip = truncated_reciprocal(J, n)

    CJ = (C*J).quo_rem(z_mod)[1]

    int_CJ = CJ.integral()

    f = (J_recip*(alpha + int_CJ)).quo_rem(z*z_mod)[1]

    return f


def compute_pe_quadratic(R, A, B, ell):
    r"""
    Computes the truncated Weierstrass function of an elliptic curve defined by short Weierstrass model: `y^2 = x^3 + Ax + B`.
    Uses an algorithm that is of complexity `O(\ell^2)`.

    Let `p` be the characteristic of the underlying field: Then we must have either `p=0`, or `p > 2\ell + 3`.

    INPUT:

     - ``poly_ring`` - polynomial ring, to compute the `\wp` function in (assumes that the generator is `z^2` for efficiency of storage/operations.)
     - ``A``         - field element corresponding to the `x` coefficient in the Weierstrass equation of an elliptic curve
     - ``B``         - field element corresponding to the constant coefficient in the Weierstrass equation of an elliptic curve
     - ``ell``       - degree of `z` to compute the truncated function to.  If `p` is the characteristic of the underlying field.     If `p > 0` then we must have `2\ell + 3 < p`.

    OUTPUT:

    polynomial -- the element in ``poly_ring`` that corresponds to the truncated function to precision `2\ell`.

    ALGORITHM:

    This function uses the algorithm described in section 3.2 of [BMSS].

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, Schost, "Fast Algorithms for Isogenies."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_pe_quadratic

        sage: R.<x> = GF(37)[]
        sage: compute_pe_quadratic(R, GF(37)(1), GF(37)(8), 7)
        (7*x^7 + 9*x^5 + x^4 + 20*x^3 + 22*x^2 + 1)/x

        sage: R.<x> = QQ[]
        sage: compute_pe_quadratic(R, 1, 8, 7)
        (-16/5775*x^7 + 23902/238875*x^6 + 24/385*x^5 + 1/75*x^4 - 8/7*x^3 - 1/5*x^2 + 1)/x

    """

    c = [0 for j in xrange(ell)]
    c[0] = -A/5
    c[1] = -B/7

    Z = R.gen()
    pe = Z**-1 + c[0]*Z + c[1]*Z**2

    for k in xrange(3, ell):
        t = 0
        for j in xrange(1,k-1):
            t += c[j-1]*c[k-2-j]
        ck = (3*t)/((k-2)*(2*k+3))
        pe += ck*Z**k
        c[k-1] = ck

    return pe


def compute_pe_fast(poly_ring, A, B, ell):
    r"""
    Computes the truncated Weierstrass function of an elliptic curve defined by short Weierstrass model: `y^2 = x^3 + Ax + B`.
    It does this with time complexity `O(M(n))`.

    Let `p` be the characteristic of the underlying field: Then we must have either `p=0`, or `p > 2\ell + 3`.

    INPUT:

     - ``poly_ring`` - polynomial ring, to compute the function in (assumes that the generator is `z^2` for efficiency of storage/operations.)
     - ``A``         - field element corresponding to the `x` coefficient in the Weierstrass equation of an elliptic curve
     - ``B``         - field element corresponding to the constant coefficient in the Weierstrass equation of an elliptic curve
     - ``ell``       - degree of `z` to compute the truncated function to.  If `p` is the characteristic of the underlying field and `p > 0`, then we must have `2\ell + 3 < p`.

    OUTPUT:

    polynomial -- the element in ``poly_ring`` that corresponds to the truncated `\wp` function to precision `2\ell`.

    ALGORITHM:

    This function uses the algorithm described in section 3.3 of
    [BMSS].

    .. note::

       Occasionally this function will fail to give the right answer,
       it faithfully implements the above published algorithm, so
       compute_pe_quadratic should be used as a fallback.

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, Schost, "Fast Algorithms for Isogenies."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_pe_fast
        sage: R.<x> = QQ[]
        sage: compute_pe_fast(R, 1, 8, 7)
        (-16/5775*x^7 + 23902/238875*x^6 + 24/385*x^5 + 1/75*x^4 - 8/7*x^3 - 1/5*x^2 + 1)/x

        sage: R.<x> = GF(37)[]
        sage: compute_pe_fast(R, GF(37)(1), GF(37)(8), 5)
        (9*x^5 + x^4 + 20*x^3 + 22*x^2 + 1)/x

    """

    z = poly_ring.gen()

    f1 = z

    s = 2

    # solve the nonlinear differential equation

    n = 2*ell + 4

    while (s < n):
        f1pr = f1.derivative()

        next_s = 2*s - 1

        if ( n < next_s):
            next_s = n

        a = 2*f1pr
        b = -(6*B*(f1**5) + 4*A*(f1**3))
        c = B*(f1**6) + A*f1**4 + 1 - (f1pr**2)

        z_mod = z**(next_s)

        a = a.quo_rem(z_mod)[1]
        b = b.quo_rem(z_mod)[1]
        c = c.quo_rem(z_mod)[1]

        f2 = solve_linear_differential_system(a, b, c, 0, z, next_s)

        f1 = f1 + f2

        s = next_s

    R = f1


    Q = (R**2).quo_rem(z**(2*ell+5))[1]

    pe_denom = Q.quo_rem(z**2)[0]

    pe_numer1 = truncated_reciprocal(pe_denom, 2*ell+1)

    pe_numer1_list = pe_numer1.list()

    pe_numer2 = 0

    # now we go through and make this a polynomial in the even powers
    for j in xrange(ell+1):
        pe_numer2 = pe_numer2*z + pe_numer1_list[2*(ell - j)]

    pe = pe_numer2/(z)

    return pe


def compute_pe(R, E, ell, algorithm=None):
    r"""
    Computes the truncated Weierstrass function on an elliptic curve defined by short Weierstrass model: `y^2 = x^3 + Ax + B`.
    Uses the algorithm specified by the algorithm parameter.

    INPUT:

    - ``R``         - polynomial ring, the ring to compute the truncated `\wp` function in (treats the generator as a power of 2 for ease of storage/use.)
    - ``E``         - Elliptic Curve, must be in short Weierstrass form `0 = a_1 =  a_2 = a_3`
    - ``ell``       - precision to compute the truncated function to
    - ``algorithm`` - string (default:None) an algorithm identifier indicating using the "fast" or "quadratic" algorithm. If the algorithm is None, then this function determines the best algorithm to use.

    OUTPUT:

    polynomial - a polynomial corresponding to the truncated Weierstrass `\wp` function in ``R``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_pe
        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 10) * (x + 12) * (x + 16)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: pe = compute_pe(R, E, 7, algorithm="quadratic")
        sage: pe = pe(x^2)
        sage: pe
        (7*x^14 + 9*x^10 + x^8 + 20*x^6 + 22*x^4 + 1)/x^2

    """
    Ea_inv = E.a_invariants()
    if ( (0,0,0) != (Ea_inv[0], Ea_inv[1], Ea_inv[2]) ):
        raise ValueError, "elliptic curve parameter must be in short Weierstrass form"

    p = R.base_ring().characteristic()

    # if the algorithm is not set, try to determine algorithm from input
    if (None == algorithm):
        if (0 == p) or (p < 2*ell + 5):
            algorithm = "fast"
        elif (p < 2*ell + 3):
            algorithm = "quadratic"
        else:
            raise NotImplementedError, "algorithms for computing pe function for that characteristic / precision pair not implemented."

    A = Ea_inv[3]
    B = Ea_inv[4]

    if ("quadratic"==algorithm):

        if (0 < p) and (p < 2*ell + 3):
            raise ValueError, "For computing pe via quadratic algorithm, characteristic of underlying field must be greater than or equal to 2*ell + 3"

        pe = compute_pe_quadratic(R, A, B, ell)

    elif ("fast"==algorithm):

        if (0 < p) and (p < 2*ell + 4):
            raise ValueError, "For computing pe via the fast algorithm, characteristic of underlying field must be greater than or equal to 2*ell + 4"

        pe = compute_pe_fast(R, A, B, ell)

    else:
        raise ValueError, "unknown algorithm for computing pe."

    return pe


def starks_find_r_and_t(T, Z):
    r"""
    Helper function for starks algorithm.

    INPUT:

    - ``T`` - Rational function in ``Z``.
    - ``Z`` - Variable of the rational function ``T``.

    OUTPUT:

    tuple -- `(r, t)` where `r` is the largest exponent such that the
    coefficient `t` of `1/Z^r` is nonzero, where `T` is regarded as
    a polynomial in `1/Z`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import starks_find_r_and_t
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_pe
        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 10) * (x + 12) * (x + 16)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: pe = compute_pe(R, E, 7, algorithm="quadratic")
        sage: pe = pe(x^2)
        sage: starks_find_r_and_t(pe, x)
        (2, 1)

    """
    Tdenom = T.denominator()
    Tdenom_lc = Tdenom.leading_coefficient()
    Tnumer = (1/Tdenom_lc)*T.numerator()

    r = Tdenom.degree()

    t_r = Tnumer.constant_coefficient()

    if (0 == r) and (t_r == 0):
        Tnumer_quo = Tnumer
        Tnumer_rem = 0
        while (0 == Tnumer_rem):
            (Tnumer_quo, Tnumer_rem) = Tnumer_quo.quo_rem(Z)
            r -= 1
        r += 1
        t_r = Tnumer_rem

    return (r, t_r)


def compute_isogeny_starks(E1, E2, ell, pe_algorithm="fast"):
    r"""
    Computes the degree ``ell`` isogeny between ``E1`` and ``E2`` via
    Stark's algorithm.  There must be a degree ``ell``, separable,
    normalized isogeny from ``E1`` to ``E2``.

    INPUT:

    - ``E1``  - an elliptic curve in short Weierstrass form.
    - ``E2``  - an elliptic curve in short Weierstrass form.
    - ``ell`` - the degree of the isogeny from E1 to E2.

    OUTPUT:

    polynomial -- over the field of definition of ``E1``, ``E2``, that is the kernel polynomial of the isogeny from ``E1`` to ``E2``.

    ALGORITHM:

    This function uses Starks Algorithm as presented in section 6.2 of
    [BMSS].

    .. note::

       As published there, the algorithm is incorrect, and a correct
       version (with slightly different notation) can be found in
       [M09].  The algorithm originates in [S72]

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, SCHOST, "Fast Algorithms for Isogenies."
    - [M09] Moody, "The Diffie-Hellman Problem and Generalization of Verheul's Theorem"
    - [S72] Stark, "Class-numbers of complex quadratic fields."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks, compute_sequence_of_maps

        sage: E = EllipticCurve(GF(97), [1,0,1,1,0])
        sage: R.<x> = GF(97)[]; f = x^5 + 27*x^4 + 61*x^3 + 58*x^2 + 28*x + 21
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: (isom1, isom2, E1pr, E2pr, ker_poly) = compute_sequence_of_maps(E, E2, 11)
        sage: compute_isogeny_starks(E1pr, E2pr, 11, pe_algorithm="quadratic")
        z^10 + 37*z^9 + 53*z^8 + 66*z^7 + 66*z^6 + 17*z^5 + 57*z^4 + 6*z^3 + 89*z^2 + 53*z + 8

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_starks(E, E2, 5)
        z^4 + 14*z^3 + z^2 + 34*z + 21
        sage: f**2
        x^4 + 14*x^3 + x^2 + 34*x + 21

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: R.<x> = QQ[]
        sage: f = x
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_starks(E, E2, 2, pe_algorithm="fast")
        z

    """

    K = E1.base_field()

    R = PolynomialRing(K, 'z')
    z = R.gen()

    S = PolynomialRing(K, 'Z')
    Z = S.gen()

    pe1 = compute_pe(S, E1, 2*ell, pe_algorithm)
    pe2 = compute_pe(S, E2, 2*ell, pe_algorithm)

    n = 1
    q = [R(1), R(0)]

    T = pe2

    while ( q[n].degree() < (ell-1) ):
#        print 'n=', n
#        print 'T = ', T

        n = n + 1
        a_n = 0

        (r, t_r) = starks_find_r_and_t(T, Z)

        while ( 0 <= r ):
#            print 'r=',r
#            print 't_r=',t_r
            a_n = a_n + t_r*z**r
            T = T - t_r*pe1**r
            (r, t_r) = starks_find_r_and_t(T, Z)

        q_n = a_n*q[n-1] + q[n-2]
        q.append(q_n)

        if (n == ell+1):
            break

        (r, t_r) = starks_find_r_and_t(T, Z)

        Tdenom_lc = T.denominator().leading_coefficient()
        Tnumer = (1/Tdenom_lc)*T.numerator()

#        print 'Tnumer before divide out=', Tnumer

        Tdenom_next = Z**(-r)

#        print 'r=', r

#        print 'Tdenom_next = ', Tdenom_next

        # compute the highest power of z**2 that divides the numerator
        (Tnumer, Tnumer_rem) = Tnumer.quo_rem(Tdenom_next);

#        if (0 != Tnumer_rem):
#            print 'ERROR expected remainder =0 was not 0.'

#        print 'Tnumer after divide out = ', Tnumer
        Tnumer_next = truncated_reciprocal(Tnumer, 2*ell)
#        print "recip check", Tnumer*Tnumer_next

#        print 'Tnumer_next =', Tnumer_next
#        print 'Tdenom_next =', Tdenom_next

        T = Tnumer_next/Tdenom_next

#        print 'q_n=', q_n
#        print 'a_n=', a_n

#    print q

    qn = q[n]
    qn = (1/qn.leading_coefficient())*qn

    return qn


def split_kernel_polynomial(E1, ker_poly, ell):
    r"""
    Internal helper function for ``compute_isogeny_kernel_polynomial``.

    Given a full kernel polynomial (where two torsion `x`-coordinates
    are roots of multiplicity 1, and all other roots have multiplicity
    2.)  of degree `\ell-1`, returns the maximum separable divisor.
    (i.e. the kernel polynomial with roots of multiplicity at most 1).

    EXAMPLES:

    The following example implicitly exercises this function::

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 10) * (x + 12) * (x + 16)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import split_kernel_polynomial
        sage: ker_poly = compute_isogeny_starks(E, E2, 7, pe_algorithm="quadratic"); ker_poly
        z^6 + 2*z^5 + 20*z^4 + 11*z^3 + 36*z^2 + 35*z + 16
        sage: split_kernel_polynomial(E, ker_poly, 7)
        z^3 + z^2 + 28*z + 33

    """

    poly_ring = ker_poly.parent()

    z = poly_ring.gen(0)

    ker_poly_2tor = two_torsion_part(E1, poly_ring, ker_poly, ell)
    ker_poly_quo = poly_ring(ker_poly/ker_poly_2tor)
    ker_poly_quo_sqrt = ker_poly_quo.gcd(ker_poly_quo.derivative(z))
    ker_poly = ker_poly_2tor*ker_poly_quo_sqrt
    ker_poly = (1/ker_poly.leading_coefficient())*ker_poly

    return ker_poly


def compute_isogeny_kernel_polynomial(E1, E2, ell, algorithm="starks"):
    r"""
    Computes the degree ``ell`` isogeny between ``E1`` and ``E2``.

    There must be a degree ``ell``, separable, normalized isogeny from
    ``E1`` to ``E2``.  If no algorithm is specified, this function
    determines the best algorithm to use.

    INPUT:

    - ``E1``        - an elliptic curve in short Weierstrass form.

    - ``E2``        - an elliptic curve in short Weierstrass form.

    - ``ell``       - the degree of the isogeny from ``E1`` to ``E2``.

    - ``algorithm`` - string (default:"starks") if None, this function automatically determines best algorithm to use.
                   Otherwise uses the algorithm specified by the string.  Valid values are "starks" or "fastElkies"

    OUTPUT:

    polynomial -- over the field of definition of ``E1``, ``E2``, that is the kernel polynomial of the isogeny from ``E1`` to ``E2``.

    .. note::

       When using Stark's algorithm, occasionally the fast pe
       computation fails, so we retry with the quadratic algorithm, which works in all cases (for valid inputs.)

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_kernel_polynomial

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_kernel_polynomial(E, E2, 5)
        z^2 + 7*z + 13
        sage: f
        x^2 + 7*x + 13

        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_isogeny_kernel_polynomial(E, E2, 4)
        z^3 + z

    """

    if ("starks"==algorithm):
        ker_poly = compute_isogeny_starks(E1, E2, ell)
    elif ("fastElkies"==algorithm):
        raise NotImplementedError
    else:
        raise ValueError, "algorithm parameter was for an unknown algorithm."

    #
    # Everything that follows here is how we get the kernel polynomial in the form we want
    # i.e. a separable polynomial (no repeated roots.)
    #
    ker_poly = split_kernel_polynomial(E1, ker_poly, ell)

    # in case we catastrophically failed using the fast pe algorithm fall back to quadratic algorithm
    # this is a bug and should be fixed, but this is a work around for now
    if (ker_poly.degree() < ell/2) and ("starks"==algorithm):
        ker_poly = compute_isogeny_starks(E1, E2, ell, pe_algorithm="quadratic")
        ker_poly = split_kernel_polynomial(E1, ker_poly, ell)

    return ker_poly


def compute_intermediate_curves(E1, E2):
    r"""
    Computes isomorphism from ``E1`` to an intermediate domain and an
    isomorphism from an intermediate codomain to ``E2``.

    Intermediate domain and intermediate codomain, are in short
    Weierstrass form.

    This is used so we can compute `\wp` functions from the short Weierstrass model more easily.

    The underlying field must be of characteristic not equal to 2,3.

    INPUT:

    - ``E1`` - an elliptic curve
    - ``E2`` - an elliptic curve

    OUTPUT:

    tuple -- (``pre_isomorphism``, ``post_isomorphism``, ``intermediate_domain``, ``intermediate_codomain``):

    - ``intermediate_domain``: a short Weierstrass model isomorphic to ``E1``
    - ``intermediate_codomain``: a short Weierstrass model isomorphic to ``E2``
    - ``pre_isomorphism``: normalized isomorphism from ``E1`` to intermediate_domain
    - ``post_isomorphism``: normalized isomorphism from intermediate_codomain to ``E2``

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_intermediate_curves
        sage: E = EllipticCurve(GF(83), [1,0,1,1,0])
        sage: R.<x> = GF(83)[]; f = x+24
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_intermediate_curves(E, E2)
        (Elliptic Curve defined by y^2 = x^3 + 62*x + 74 over Finite Field of size 83,
         Elliptic Curve defined by y^2 = x^3 + 65*x + 69 over Finite Field of size 83,
         Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 83
          To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 62*x + 74 over Finite Field of size 83
          Via:  (u,r,s,t) = (1, 76, 41, 3),
         Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 65*x + 69 over Finite Field of size 83
          To:   Abelian group of points on Elliptic Curve defined by y^2 + x*y + 77*y = x^3 + 49*x + 28 over Finite Field of size 83
          Via:  (u,r,s,t) = (1, 7, 42, 80))

        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_intermediate_curves(E, E2)
        (Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1,
         Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1,
         Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1
          To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1
          To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0))

    """

    if (E1.base_ring().characteristic() in [2,3]):
        raise NotImplemented

    # compute the r,s,t values that clear the denominator of E1
    a1 = E1.a1()
    a2 = E1.a2()
    a3 = E1.a3()

    s1 = -a1/2
    r1 = (s1**2 + s1*a1 - a2)/3
    t1 = (-r1*a1 - a3)/2

    # compute the isomorphism from E1 to intermediate_domain
    pre_isom = WeierstrassIsomorphism(E1, (1, r1, s1, t1))

    intermediate_domain = pre_isom.codomain().codomain()

    # compute the r,s,t values that clear the denominator of E2
    a1pr = E2.a1()
    a2pr = E2.a2()
    a3pr = E2.a3()

    s2 = -a1pr/2
    r2 = (s2**2 + s2*a1pr - a2pr)/3
    t2 = (-r2*a1pr - a3pr)/2

    post_isom_inv = WeierstrassIsomorphism(E2, (1, r2, s2, t2))
    intermediate_codomain = post_isom_inv.codomain().codomain()

    post_isom = WeierstrassIsomorphism(intermediate_codomain, (1, -r2, -s2, -t2))

    return (intermediate_domain, intermediate_codomain, pre_isom, post_isom)


def compute_sequence_of_maps(E1, E2, ell):
    r"""
    Given domain ``E1`` and codomain ``E2`` such that there is a degree ``ell`` separable normalized isogeny from ``E1`` to ``E2``,    returns pre/post isomorphism, as well as intermediate domain and codomain, and kernel polynomial.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_sequence_of_maps
        sage: E = EllipticCurve('11a1')
        sage: R.<x> = QQ[]; f = x^2 - 21*x + 80
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_sequence_of_maps(E, E2, 5)
        (Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
          To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 - 31/3*x - 2501/108 over Rational Field
          Via:  (u,r,s,t) = (1, 1/3, 0, -1/2),
         Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 - 23461/3*x - 28748141/108 over Rational Field
          To:   Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
          Via:  (u,r,s,t) = (1, -1/3, 0, 1/2),
         Elliptic Curve defined by y^2 = x^3 - 31/3*x - 2501/108 over Rational Field,
         Elliptic Curve defined by y^2 = x^3 - 23461/3*x - 28748141/108 over Rational Field,
         z^2 - 61/3*z + 658/9)

        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_sequence_of_maps(E, E2, 4)
        (Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1
          To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1
          To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1
          Via:  (u,r,s,t) = (1, 0, 0, 0),
         Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1,
         Elliptic Curve defined by y^2 = x^3 + 16*x over Number Field in i with defining polynomial x^2 + 1,
         z^3 + z)

        sage: E = EllipticCurve(GF(97), [1,0,1,1,0])
        sage: R.<x> = GF(97)[]; f = x^5 + 27*x^4 + 61*x^3 + 58*x^2 + 28*x + 21
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_sequence_of_maps(E, E2, 11)
        (Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 + x*y + y = x^3 + x over Finite Field of size 97
          To:   Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 52*x + 31 over Finite Field of size 97
          Via:  (u,r,s,t) = (1, 8, 48, 44),
         Generic morphism:
          From: Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 41*x + 66 over Finite Field of size 97
          To:   Abelian group of points on Elliptic Curve defined by y^2 + x*y + 9*y = x^3 + 83*x + 6 over Finite Field of size 97
          Via:  (u,r,s,t) = (1, 89, 49, 53),
         Elliptic Curve defined by y^2 = x^3 + 52*x + 31 over Finite Field of size 97,
         Elliptic Curve defined by y^2 = x^3 + 41*x + 66 over Finite Field of size 97,
         z^5 + 67*z^4 + 13*z^3 + 35*z^2 + 77*z + 69)

    """

    (E1pr, E2pr, pre_isom, post_isom) = compute_intermediate_curves(E1, E2)

    ker_poly = compute_isogeny_kernel_polynomial(E1pr, E2pr, ell)

    return (pre_isom, post_isom, E1pr, E2pr, ker_poly)
