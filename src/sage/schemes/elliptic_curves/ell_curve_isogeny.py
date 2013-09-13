r"""
Isogenies

An isogeny `\varphi: E_1\to E_2` between two elliptic curves `E_1` and
`E_2` is a morphism of curves that sends the origin of `E_1` to the
origin of `E_2`. Such a morphism is automatically a morphism of group
schemes and the kernel is a finite subgroup scheme of `E_1`.  Such a
subscheme can either be given by a list of generators, which have to
be torsion points, or by a polynomial in the coordinate `x` of the
Weierstrass equation of `E_1`.

The usual way to create and work with isogenies is illustrated with
the following example::

    sage: k = GF(11)
    sage: E = EllipticCurve(k,[1,1])
    sage: Q = E(6,5)
    sage: phi = E.isogeny(Q)
    sage: phi
    Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11 to Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 11
    sage: P = E(4,5)
    sage: phi(P)
    (10 : 0 : 1)
    sage: phi.codomain()
    Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 11
    sage: phi.rational_maps()
    ((x^7 + 4*x^6 - 3*x^5 - 2*x^4 - 3*x^3 + 3*x^2 + x - 2)/(x^6 + 4*x^5 - 4*x^4 - 5*x^3 + 5*x^2), (x^9*y - 5*x^8*y - x^7*y + x^5*y - x^4*y - 5*x^3*y - 5*x^2*y - 2*x*y - 5*y)/(x^9 - 5*x^8 + 4*x^6 - 3*x^4 + 2*x^3))

The functions directly accessible from an elliptic curve ``E`` over a
field are ``isogeny`` and ``isogeny_codomain``.

The most useful functions that apply to isogenies are

- ``codomain``
- ``degree``
- ``domain``
- ``dual``
- ``rational_maps``
- ``kernel_polynomial``

.. Warning::

   Only cyclic isogenies are implemented (except for [2]). Some
   algorithms may need the isogeny to be normalized.

AUTHORS:

- Daniel Shumow <shumow@gmail.com>: 2009-04-19: initial version

- Chris Wuthrich : 7/09: changes: add check of input, not the full list is needed.
  10/09: eliminating some bugs.

- John Cremona and Jenny Cooley: 2009-07..11: implement `l`-isogenies
  for `l` = 2, 3, 5, 7 13 (the genus 0 cases) and also for `l` = 11,
  17, 19, 37, 43, 67 or 163 over `\QQ` (the sporadic cases with only
  finitely many `j`-invariants each).
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

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.all import Integer, ZZ
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.all import is_Polynomial
from sage.schemes.elliptic_curves.all import EllipticCurve
from sage.schemes.elliptic_curves.all import is_EllipticCurve

from sage.rings.number_field.number_field_base import is_NumberField

from sage.rings.rational_field import is_RationalField, QQ

from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

from sage.sets.set import Set

from sage.misc.cachefunc import cached_function

#
# Private function for parsing input to determine the type of
# algorithm
#
def isogeny_determine_algorithm(E, kernel, codomain, degree, model):
    r"""
    Helper function that allows the various isogeny functions to infer
    the algorithm type from the parameters passed in.

    If ``kernel`` is a list of points on the EllipticCurve `E`, then
    we assume the algorithm to use is Velu.

    If ``kernel`` is a list of coefficients or a univariate polynomial
    we try to use the Kohel's algorithms.

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

    - ``kernel`` - Either a list of points in the kernel of the isogeny, or a
                   kernel polynomial (specified as a either a univariate
                   polynomial or a coefficient list.)

    - ``degree`` - an integer, (default:``None``)  optionally specified degree
                   of the kernel.

    OUTPUT:

    (elliptic curve) the codomain of the separable normalized isogeny
    from this kernel

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
    Given parameters `v` and `w` (as in Velu / Kohel / etc formulas)
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
    The formula for computing `v` and `w` using Kohel's formulas for
    isogenies of degree 2.

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
    The formula for computing `v` and `w` using Kohel's formulas for
    isogenies of degree 3.

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
    This function computes the codomain from the kernel polynomial as
    per Kohel's formulas.

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

    Returns the greatest common divisor of ``psi`` and the 2 torsion
    polynomial of `E`.

    EXAMPLES:

    Every function that computes the kernel polynomial via Kohel's
    formulas will call this function::

        sage: E = EllipticCurve(GF(19), [1,2,3,4,5])
        sage: R.<x> = GF(19)[]
        sage: phi = EllipticCurveIsogeny(E, x + 13)
        sage: isogeny_codomain_from_kernel(E, x + 13) == phi.codomain()
        True
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import two_torsion_part
        sage: two_torsion_part(E, R, x+13, 2)
        x + 13

    """
    if (None==degree) or (0 == degree % 2):

        x = poly_ring.gens()[0]
        psi_2 = E.two_division_polynomial(x)
        psi_G = poly_ring(psi.gcd(psi_2))

    else:

        psi_G = poly_ring(1)

    return psi_G


class EllipticCurveIsogeny(Morphism):
    r"""
    Class Implementing Isogenies of Elliptic Curves

    This class implements cyclic, separable, normalized isogenies of
    elliptic curves.

    Several different algorithms for computing isogenies are
    available.  These include:

    - Velu's Formulas: Velu's original formulas for computing
      isogenies.  This algorithm is selected by giving as the
      ``kernel`` parameter a list of points which generate a finite
      subgroup.

    - Kohel's Formulas: Kohel's original formulas for computing
      isogenies.  This algorithm is selected by giving as the
      ``kernel`` parameter a monic polynomial (or a coefficient list
      (little endian)) which will define the kernel of the isogeny.

    INPUT:

    - ``E``         - an elliptic curve, the domain of the isogeny to
                      initialize.

    - ``kernel``    - a kernel, either a point in ``E``, a list of points
                      in ``E``, a monic kernel polynomial, or ``None``.
                      If initializing from a domain/codomain, this must be
                      set to None.

    - ``codomain``  - an elliptic curve (default:``None``).  If ``kernel``
                      is ``None``, then this must be the codomain of a cyclic,
                      separable, normalized isogeny, furthermore, ``degree``
                      must be the degree of the isogeny from ``E`` to
                      ``codomain``. If ``kernel`` is not ``None``, then this
                      must be isomorphic to the codomain of the cyclic normalized
                      separable isogeny defined by ``kernel``, in this case, the
                      isogeny is post composed with an isomorphism so that this
                      parameter is the codomain.

    - ``degree``    - an integer (default:``None``).
                      If ``kernel`` is ``None``, then this is the degree of the
                      isogeny from ``E`` to ``codomain``.
                      If ``kernel`` is not ``None``, then this is used to determine
                      whether or not to skip a gcd of the kernel polynomial with the
                      two torsion polynomial of ``E``.

    - ``model``     - a string (default:``None``).  Only supported variable is
                      ``minimal``, in which case if ``E`` is a curve over the
                      rationals, then the codomain is set to be the unique global
                      minimum model.

    - ``check`` (default: ``True``) checks if the input is valid to define an isogeny

    EXAMPLES:

    A simple example of creating an isogeny of a field of small
    characteristic::

        sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
        sage: phi = EllipticCurveIsogeny(E, E((0,0)) ); phi
        Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7 to Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 7
        sage: phi.degree() == 2
        True
        sage: phi.kernel_polynomial()
        x
        sage: phi.rational_maps()
        ((x^2 + 1)/x, (x^2*y - y)/x^2)
        sage: phi == loads(dumps(phi))   # not tested - pickling http://trac.sagemath.org/sage_trac/ticket/11599
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

    We can create an isogeny that has kernel equal to the full 2
    torsion::

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

    A more complicated example over the rationals (of odd degree)::

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

    We can also do this same example over the number field defined by
    the irreducible two torsion polynomial of `E`::

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

    The following example shows how to specify an isogeny from domain
    and codomain::

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

    However only cyclic normalized isogenies can be constructed this
    way. So it won't find the isogeny [3]::

        sage: E.isogeny(None, codomain=E,degree=9)
        Traceback (most recent call last):
        ...
        ValueError: The two curves are not linked by a cyclic normalized isogeny of degree 9

    Also the presumed isogeny between the domain and codomain must be
    normalized::

        sage: E2.isogeny(None,codomain=E,degree=5)
        Traceback (most recent call last):
        ...
        ValueError: The two curves are not linked by a cyclic normalized isogeny of degree 5
        sage: phi.dual()
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: phi.dual().is_normalized()
        False

    Here an example of a construction of a endomorphisms with cyclic
    kernel on a CM-curve::

        sage: K.<i> = NumberField(x^2+1)
        sage: E = EllipticCurve(K, [1,0])
        sage: RK.<X> = K[]
        sage: f = X^2 - 2/5*i + 1/5
        sage: phi= E.isogeny(f)
        sage: isom = phi.codomain().isomorphism_to(E)
        sage: phi.set_post_isomorphism(isom)
        sage: phi.codomain() == phi.domain()
        True
        sage: phi.rational_maps()
        (((4/25*i + 3/25)*x^5 + (4/5*i - 2/5)*x^3 - x)/(x^4 + (-4/5*i + 2/5)*x^2 + (-4/25*i - 3/25)),
         ((11/125*i + 2/125)*x^6*y + (-23/125*i + 64/125)*x^4*y + (141/125*i + 162/125)*x^2*y + (3/25*i - 4/25)*y)/(x^6 + (-6/5*i + 3/5)*x^4 + (-12/25*i - 9/25)*x^2 + (2/125*i - 11/125)))
    """

    ####################
    # member variables
    ####################

    __check = None
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

        if type(kernel) != type([1,1]) and kernel in E :
            # a single point was given, we put it in a list
            # the first condition assures that [1,1] is treated as x+1
            kernel = [kernel]

        # if the kernel is None and the codomain isn't
        # calculate the kernel polynomial
        pre_isom = None
        post_isom = None

        self.__check = check

        if (kernel is None) and (codomain is not None):

            if (degree is None):
                raise ValueError, "If specifying isogeny by domain and codomain, degree parameter must be set."

            # save the domain/codomain: really used now (trac #7096)
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

        if (pre_isom is not None):
            self.set_pre_isomorphism(pre_isom)

        if (post_isom is not None):
            self.__set_post_isomorphism(old_codomain, post_isom)   #(trac #7096)

        # Inheritance house keeping

        self.__perform_inheritance_housekeeping()

        return


    def __call__(self, P, output_base_ring=None):
        r"""
        Function that implements the call-ability of elliptic curve
        isogenies.

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

            sage: E = EllipticCurve(GF(17), [0,-1,0,-3,-1])
            sage: Q = E((16,0))
            sage: tau = E.isogeny([Q],E)
            sage: tau(Q)
            (0 : 1 : 0)

        TESTS (trac 10888)::

            sage: K.<th> = NumberField(x^2+3)
            sage: E = EllipticCurve(K,[7,0])
            sage: phi = E.isogeny(E(0,0))
            sage: P = E(-3,4*th)
            sage: phi(P)
            (-16/3 : 8/9*th : 1)
            sage: Q = phi(P)
            sage: phihat = phi.dual()
            sage: phihat(Q)
            (-1/48 : 127/576*th : 1)

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
        if (self.__pre_isomorphism is not None):
            temp_xP = self.__prei_x_coord_ratl_map(xP, yP)
            temp_yP = self.__prei_y_coord_ratl_map(xP, yP)
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
        if (self.__post_isomorphism is not None):
            tempX = self.__posti_x_coord_ratl_map(outP[0], outP[1])
            tempY = self.__posti_y_coord_ratl_map(outP[0], outP[1])
            outP = [tempX, tempY]

        if change_output_ring:
            if (output_base_ring is None):
                output_base_ring = E_P.base_ring()
            outE2 = self.__E2.change_ring(output_base_ring)
        else:
            output_base_ring = self.__E2.base_ring()
            outE2 = self.__E2
            outP = self.__E2.point(outP,check=False)

        R = output_base_ring

        return outE2.point([R(outP[0]), R(outP[1]), R(1)], check=False)


    def __getitem__(self, i):
        self.__initialize_rational_maps()
        if (i < 0) or (i > 2):
            raise IndexError

        if i == 0:
            return self.__X_coord_rational_map
        else:
            return self.__Y_coord_rational_map

    def __iter__(self):
        self.__initialize_rational_maps()
        return iter((self.__X_coord_rational_map, self.__Y_coord_rational_map))

    def __hash__(self):
        r"""
        Function that implements the hash ability of Isogeny objects.

        This hashes the underlying kernel polynomial so that equal
        isogeny objects have the same hash value.  Also, this hashes
        the base field, and domain and codomain curves as well, so
        that isogenies with the same kernel polynomial (over different
        base fields / curves) hash to different values.

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

            sage: E = EllipticCurve('49a3')
            sage: R.<X> = QQ[]
            sage: EllipticCurveIsogeny(E,X^3-13*X^2-58*X+503,check=False)
            Isogeny of degree 7 from Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 107*x + 552 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 5252*x - 178837 over Rational Field

        """

        if (self.__this_hash is not None):
            return self.__this_hash

        ker_poly_list = self.__kernel_polynomial_list

        if (ker_poly_list is None):
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

        This function works by comparing the underlying kernel
        objects.

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
            sage: E = EllipticCurve('11a1')
            sage: phi = E.isogeny(E(5,5))
            sage: psi = E.isogeny(phi.kernel_polynomial())
            sage: phi == psi
            True
            sage: phi.dual() == psi.dual()
            True


        """
        if (not isinstance(other, EllipticCurveIsogeny)):
            return -1

        if (self.__kernel_polynomial is None):
            self.__init_kernel_polynomial()

        return cmp(self.__kernel_polynomial, other.kernel_polynomial())


    def __neg__(self):
        r"""
        Function to implement unary negation (-) operator on
        isogenies. Returns a copy of this isogeny that has been
        negated.

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
        Special sage specific function that implement the
        functionality to display the isogeny self as a string.

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
        Special sage specific function that implements functionality
        to display an isogeny object as a latex string.

        This function returns a latex string representing the isogeny
        self as the `x` and `y` coordinate rational functions.

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
        A private function to clear the hash if the codomain has been
        modified by a pre or post isomorphism.

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
        Internal helper function, sets values on the super classes of
        this class.

        EXAMPLES:

        The following examples will implicitly exercise this
        function::

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
        An internal function for EllipticCurveIsogeny objects that
        sets up the member variables necessary for algebra.

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
    def __initialize_rational_maps(self, precomputed_maps=None):
        r"""
        Private function that computes and initializes the rational
        maps.

        INPUT:

        - ``

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

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: Eshort = E.short_weierstrass_model()
            sage: phi = E.isogeny(E(0), Eshort)
            sage: phiX, phiY = phi.rational_maps()
            sage: phiX(1,2), phiY(1,2)
            (63, 864)
        """
        if self.__rational_maps_initialized:
            return

        if precomputed_maps is None:
            if ("velu"==self.__algorithm):
                (X_map, Y_map) = self.__initialize_rational_maps_via_velu()
            if ("kohel"==self.__algorithm):
                (X_map, Y_map) = self.__initialize_rational_maps_via_kohel()
        else:
            X_map, Y_map = precomputed_maps

        if self.__prei_x_coord_ratl_map is not None:
            prei_X_map = self.__prei_x_coord_ratl_map
            prei_Y_map = self.__prei_y_coord_ratl_map
            X_map, Y_map = X_map.subs(x=prei_X_map, y=prei_Y_map), \
                           Y_map.subs(x=prei_X_map, y=prei_Y_map)

        if self.__posti_x_coord_ratl_map is not None:
            X_map, Y_map = \
            self.__posti_x_coord_ratl_map.subs(x=X_map, y=Y_map), \
            self.__posti_y_coord_ratl_map.subs(x=X_map, y=Y_map)

        self.__X_coord_rational_map = X_map
        self.__Y_coord_rational_map = Y_map
        self.__rational_maps_initialized = True


    def __init_kernel_polynomial(self):
        r"""
        Private function that initializes the kernel polynomial (if
        the algorithm does not take it as a parameter.)

        EXAMPLES:

        The following examples inherently exercise this function::

            sage: E = EllipticCurve(j=GF(7)(1728))
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.kernel_polynomial()
            x
            sage: phi._EllipticCurveIsogeny__init_kernel_polynomial()
            [0, 1]

        """

        if (self.__kernel_polynomial_list is not None):
            return self.__kernel_polynomial_list

        if ("velu" == self.__algorithm):
            ker_poly_list = self.__init_kernel_polynomial_velu()
        else:
            raise InputError, "The kernel polynomial should already be defined!"

        return ker_poly_list


    def __set_pre_isomorphism(self, domain, isomorphism):
        r"""
        Private function to set the pre isomorphism and domain (and
        keep track of the domain of the isogeny.)

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

        if (self.__kernel_polynomial is not None):
            ker_poly = self.__kernel_polynomial
            ker_poly = ker_poly.subs(x=self.__prei_x_coord_ratl_map)
            kp_lc = ker_poly.univariate_polynomial().leading_coefficient()
            ker_poly = (1/kp_lc)*ker_poly
            self.__kernel_polynomial = ker_poly

        self.__perform_inheritance_housekeeping()

        return;


    def __set_post_isomorphism(self, codomain, isomorphism):
        r"""
        Private function to set the post isomorphism and codomain (and
        keep track of the codomain of the isogeny.)

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
        Private function to set up the post isomorphism given the
        codomain.

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

        if (model is not None):

            if (codomain is not None):
                raise ValueError, "Cannot specify a codomain and model flag simultaneously."

            if ('minimal' == model):

                if (not is_RationalField(oldE2.base_field())):
                    raise ValueError, "specifying minimal for model flag only valid with curves over the rational numbers."

                newE2 = oldE2.minimal_model()
                post_isom = oldE2.isomorphism_to(newE2)

            else:
                raise ValueError, "Unknown value of model flag."

        elif (codomain is not None):
            if (not is_EllipticCurve(codomain)):
                raise ValueError,  "Codomain parameter must be an elliptic curve."

            if (not oldE2.is_isomorphic(codomain)):
                raise ValueError, "Codomain parameter must be isomorphic to computed codomain isogeny"

            newE2 = codomain
            post_isom = oldE2.isomorphism_to(newE2)

        if (post_isom is not None):
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
        Private function that initializes the isogeny from a list of
        points which generate the kernel (For Velu's formulas.)

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
        Private function that sorts the list of points in the kernel
        (For Velu's formulas). Sorts out the 2 torsion points, and
        puts them in a dictionary.

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
        Private function that computes the codomain via Velu's
        formulas.

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
        Private function for Velu's formulas, helper function to help
        perform the summation.

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
        Private function that sorts the list of points in the kernel
        (for Velu's formulas). Sorts out the 2 torsion points, and
        puts them in a diction

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
        # first check if the point is in the kernel
        if ( self.__kernel_2tor.has_key(xP) or self.__kernel_non2tor.has_key(xP) ) :
            return self.__intermediate_codomain(0)

        outP = self.__compute_via_velu(xP,yP)

        return outP


    def __compute_via_velu(self, xP, yP):
        r"""
        Private function for Velu's formulas, to perform the
        summation.

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
        Private function for Velu's formulas, helper function to
        initialize the rational maps.

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

        if (self.__pre_isomorphism is not None):
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
        Private function that initializes the isogeny from a kernel
        polynomial.

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

        if kernel_polynomial[-1] != 1:
            raise ValueError, "The kernel polynomial must be monic."

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
        Private function that initializes the isogeny from a kernel
        polynomial, for Kohel's algorithm in the even degree case.

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
        if  self.__check and E.division_polynomial(2, x=self.__x_var) % psi_G  != 0 :
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
        Private function that initializes the isogeny from a kernel
        polynomial.

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

            sage: E = EllipticCurve(j=-262537412640768000)
            sage: f = (E.isogenies_prime_degree()[0]).kernel_polynomial()
            sage: f.degree()
            81
            sage: E.isogeny(kernel=f)  # long time (25s on sage.math, 2012)
            Isogeny of degree 163 from Elliptic Curve defined by y^2 + y = x^3 - 2174420*x + 1234136692 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - 57772164980*x - 5344733777551611 over Rational Field

        """
        n = psi.degree()
        d = 2*n + 1

        # check if the polynomial really divides the torsion polynomial :
        if self.__check:
            alpha = psi.parent().quotient(psi).gen()
            if not E.division_polynomial(d, x=alpha).is_zero():
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
    # This is the fast omega computation that works when characteristic is not 2
    #
    def __compute_omega_fast(self, E, psi, psi_pr, phi, phi_pr):
        r"""
        Private function that initializes the omega polynomial (from
        Kohel's formulas) in the case that the characteristic of the
        underlying field is not 2.

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
        Private function that initializes the omega polynomial (from
        Kohel's formulas) in the case of general characteristic of the
        underlying field.

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

        A bug fixed in ticket #7907::

            sage: F = GF(128,'a')
            sage: a = F.gen()
            sage: E = EllipticCurve([1,0,0,0,(a**6+a**4+a**2+a)])
            sage: x = polygen(F)
            sage: ker =  (x^6 + (a^6 + a^5 + a^4 + a^3 + a^2 + a)*x^5 + (a^6 + a^5 + a^2 + 1)*x^4 + (a^6 + a^5 + a^4 + a^3 + a^2 + 1)*x^3 + (a^6 + a^3 + a)*x^2 + (a^4 + a^3 + 1)*x + a^5 + a^4 + a)
            sage: E.isogeny(ker)
            Isogeny of degree 13 from Elliptic Curve defined by y^2 + x*y = x^3 + (a^6+a^4+a^2+a) over Finite Field in a of size 2^7 to Elliptic Curve defined by y^2 + x*y = x^3 + (a^6+a^5+a^4+a^3+a^2+a)*x + (a^5+a^3) over Finite Field in a of size 2^7


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
        from sage.rings.arith import binomial

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
        Private function that computes a numeric result of this
        isogeny (via Kohel's formulas.)

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

        # first check if this is a kernel point
        # to avoid a divide by 0 error later
        if(0 == self.__inner_kernel_polynomial(x=xP)):
            return self.__intermediate_codomain(0)

        (xP_out, yP_out) = self.__compute_via_kohel(xP,yP)

        # for some dang reason in some cases
        # when the base_field is a number field
        # xP_out and yP_out do not get evaluated to field elements
        # but rather constant polynomials.
        # So in this case, we do some explicit casting to make sure
        # everything comes out right

        if is_NumberField(self.__base_field) and (1 < self.__base_field.degree()) :
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

        psi_out = self.__psi(xP,yP)
        phi_out = self.__phi(xP,yP)
        omega_out =self.__omega(xP, yP)

        psi_inv_out = 1/psi_out

        psi_inv_sq_out = psi_inv_out**2

        X_out = phi_out*(psi_inv_sq_out)
        Y_out = omega_out*(psi_inv_sq_out*psi_inv_out)

        return (X_out, Y_out)


    def __initialize_rational_maps_via_kohel(self):
        r"""
        Private function that computes and initializes the rational
        maps of this isogeny.

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
        Private function that computes and initializes the codomain of
        the isogeny (via Kohel's.)

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
        This function returns a bool indicating whether or not this
        isogeny is separable.

        This function always returns ``True`` as currently this class
        only implements separable isogenies.

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
        if (self.__kernel_polynomial is None):
            self.__init_kernel_polynomial()

        return self.__kernel_polynomial.univariate_polynomial()


    def set_pre_isomorphism(self, preWI):
        r"""
        Modifies this isogeny object to pre compose with the given
        Weierstrass isomorphism.

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

        if (self.__pre_isomorphism is None):
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
        Modifies this isogeny object to post compose with the given
        Weierstrass isomorphism.

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

        if (self.__post_isomorphism is None):
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
        Returns the pre-isomorphism of this isogeny.  If there has
        been no pre-isomorphism set, this returns ``None``.

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
        Returns the post-isomorphism of this isogeny.  If there has
        been no post-isomorphism set, this returns ``None``.

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
            Via:  (u,r,s,t) = (1, 7, 42, 80)

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

    def is_normalized(self, via_formal=True, check_by_pullback=True):
        r"""
        Returns ``True`` if this isogeny is normalized. An isogeny
        `\varphi\colon E\to E_2` between two given Weierstrass
        equations is said to be normalized if the constant `c` is `1`
        in `\varphi*(\omega_2) = c\cdot\omega`, where `\omega` and
        `omega_2` are the invariant differentials on `E` and `E_2`
        corresponding to the given equation.

        INPUT:

        - ``via_formal`` - (default: ``True``) If ``True`` it simply checks if
                           the leading term of the formal series is 1. Otherwise
                           it uses a deprecated algorithm involving the second
                           optional argument.

        - ``check_by_pullback`` -  (default:``True``) Deprecated.

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
        # easy algorithm using the formal expansion.
        if via_formal:
            phi_formal = self.formal(prec=5)
            return phi_formal[1] == 1

        # this is the old algorithm. it should be deprecated.
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
                f_normalized = True
            else:
                # For some reason, in certain cases, when the isogeny
                # is pre or post composed with a translation the
                # resulting rational functions are too complicated for
                # sage to simplify down to a constant in this case, we
                # do some cheating by checking if the post-composition
                # by isogeny has a non 1 scaling factor
                if ( inv_diff_quo.numerator().is_constant() and (inv_diff_quo.denominator().is_constant) ):
                    f_normalized = False
                else:
                    check_prepost_isomorphism = True
        else:
            check_prepost_isomorphism = True

        # If we skip checking by the pullback of the invariant
        # differential OR if that was inconclusive We explicitly check
        # if there is a post isomorphism and if it has a non 1 scaling
        # factor or if it is a just a translation.  NOTE: This only
        # works because we are using algorithms for calculating the
        # isogenies that calculate a separable normalized isogeny, if
        # this changes, this check will no longer be correct.
        #
        if (check_prepost_isomorphism):
            post_isom = self.__post_isomorphism
            if (post_isom is not None):
                if (1 == self.__base_field(post_isom.u)):
                    f_post_normalized = True
                else:
                    f_post_normalized = False
            else:
                f_post_normalized = True

            pre_isom = self.__pre_isomorphism
            if (pre_isom is not None):
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
        Computes and returns the dual isogeny of this isogeny. If
        `\varphi\colon E \to E_2` is the given isogeny, then the dual
        is by definition the unique isogeny `\hat\varphi\colon E_2\to
        E` such that the compositions `\hat\varphi\circ\varphi` and
        `\varphi\circ\hat\varphi` are the multiplication `[n]` by the
        degree of `\varphi` on `E` and `E_2` respectively.

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
            True

        Test (for trac ticket 7096)::

            sage: E = EllipticCurve('11a1')
            sage: phi = E.isogeny(E(5,5))
            sage: phi.dual().dual() == phi
            True

            sage: k = GF(103)
            sage: E = EllipticCurve(k,[11,11])
            sage: phi = E.isogeny(E(4,4))
            sage: phi
            Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 11*x + 11 over Finite Field of size 103 to Elliptic Curve defined by y^2 = x^3 + 25*x + 80 over Finite Field of size 103
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: phi.set_post_isomorphism(WeierstrassIsomorphism(phi.codomain(),(5,0,1,2)))
            sage: phi.dual().dual() == phi
            True

            sage: E = EllipticCurve(GF(103),[1,0,0,1,-1])
            sage: phi = E.isogeny(E(60,85))
            sage: phi.dual()
            Isogeny of degree 7 from Elliptic Curve defined by y^2 + x*y = x^3 + 84*x + 34 over Finite Field of size 103 to Elliptic Curve defined by y^2 + x*y = x^3 + x + 102 over Finite Field of size 103

        """

        if (self.__base_field.characteristic() in [2,3]):
            raise NotImplemented

        if (self.__dual is not None):
            return self.__dual

        # trac 7096
        (E1, E2pr, pre_isom, post_isom) = compute_intermediate_curves(self.codomain(), self.domain())

        F = self.__base_field

        d = self.__degree

        # trac 7096
        if F(d) == 0:
            raise NotImplementedError, "The dual isogeny is not separable, but only separable isogenies are implemented so far"

        # trac 7096
        # this should take care of the case when the isogeny is not normalized.
        u = self.formal(prec=5)[1]
        isom = WeierstrassIsomorphism(E2pr, (u/F(d), 0, 0, 0))

        E2 = isom.codomain().codomain()

        pre_isom = self.__E2.isomorphism_to(E1)
        post_isom = E2.isomorphism_to(self.__E1)

        phi_hat = EllipticCurveIsogeny(E1, None, E2, d)

        phi_hat.set_pre_isomorphism(pre_isom)
        phi_hat.set_post_isomorphism(post_isom)
        phi_hat.__perform_inheritance_housekeeping()

        assert phi_hat.codomain() == self.domain()

        # trac 7096 : this adjust a posteriori the automorphism
        # on the codomain of the dual isogeny.
        # we used _a_ Weierstrass isomorphism to get to the original
        # curve, but we may have to change it my an automorphism.
        # we impose that the composition has the degree
        # as a leading coefficient in the formal expansion.

        phi_sc = self.formal(prec=5)[1]
        phihat_sc = phi_hat.formal(prec=5)[1]

        sc = phi_sc * phihat_sc/F(d)

        if sc != 1:
            auts = phi_hat.codomain().automorphsims()
            aut = [a for a in auts if a.u == sc]
            if len(aut) != 1:
                raise ValueError, "There is a bug in dual()."
            phi_hat.set_post_isomorphism(a[0])

        self.__dual = phi_hat

        return phi_hat

    def formal(self,prec=20):
        r"""
        Computes the formal isogeny as a power series in the variable
        `t=-x/y` on the domain curve.

        INPUT:

        - ``prec`` - (default = 20), the precision with which the computations
                     in the formal group are carried out.

        EXAMPLES::

            sage: E = EllipticCurve(GF(13),[1,7])
            sage: phi = E.isogeny(E(10,4))
            sage: phi.formal()
            t + 12*t^13 + 2*t^17 + 8*t^19 + 2*t^21 + O(t^23)

            sage: E = EllipticCurve([0,1])
            sage: phi = E.isogeny(E(2,3))
            sage: phi.formal(prec=10)
            t + 54*t^5 + 255*t^7 + 2430*t^9 + 19278*t^11 + O(t^13)

            sage: E = EllipticCurve('11a2')
            sage: R.<x> = QQ[]
            sage: phi = E.isogeny(x^2 + 101*x + 12751/5)
            sage: phi.formal(prec=7)
            t - 2724/5*t^5 + 209046/5*t^7 - 4767/5*t^8 + 29200946/5*t^9 + O(t^10)


        """
        Eh = self.__E1.formal()
        f, g = self.rational_maps()
        xh = Eh.x(prec=prec)
        yh = Eh.y(prec=prec)
        fh = f(xh,yh)
        gh = g(xh,yh)
        return -fh/gh

    #
    # Overload Morphism methods that we want to
    #

    def _composition_(self, right, homset):
        r"""
        Composition operator function inherited from morphism class.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [E(0), E((0,1)), E((0,-1))])
            sage: phi._composition_(phi, phi.parent())
            Traceback (most recent call last):
            ...
            NotImplementedError

        The following should test that :meth:`_composition_` is called
        upon a product. However phi is currently improperly
        constructed (see :trac:`12880`), which triggers an assertion
        failure before the actual call ::

            sage: phi*phi
            Traceback (most recent call last):
            ...
            TypeError: Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7 is not in Category of hom sets in Category of Schemes

        Here would be the desired output::

            sage: phi*phi            # not tested
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def is_injective(self):
        r"""
        Method inherited from the morphism class.  Returns ``True`` if
        and only if this isogeny has trivial kernel.

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
        For elliptic curve isogenies, always returns ``True`` (as a
        non-constant map of algebraic curves must be surjective).

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
        Numerical Approximation inherited from Map (through morphism),
        nonsensical for isogenies.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,1)), E((0,-1))])
            sage: phi.n()
            Traceback (most recent call last):
            ...
            NotImplementedError: Numerical approximations do not make sense for Elliptic Curve Isogenies

        """
        raise NotImplementedError, "Numerical approximations do not make sense for Elliptic Curve Isogenies"

# no longer needed (trac 7096)
# def starks_find_r_and_t(T, Z):

def compute_isogeny_starks(E1, E2, ell):
    r"""
    Computes the degree ``ell`` isogeny between ``E1`` and ``E2`` via
    Stark's algorithm.  There must be a degree ``ell``, separable,
    normalized cyclic isogeny from ``E1`` to ``E2``.

    INPUT:

    - ``E1``  - an elliptic curve in short Weierstrass form.
    - ``E2``  - an elliptic curve in short Weierstrass form.
    - ``ell`` - the degree of the isogeny from E1 to E2.

    OUTPUT:

    polynomial -- over the field of definition of ``E1``, ``E2``, that is the
                  kernel polynomial of the isogeny from ``E1`` to ``E2``.

    ALGORITHM:

    This function uses Starks Algorithm as presented in section 6.2 of
    [BMSS].

    .. note::

       As published there, the algorithm is incorrect, and a correct
       version (with slightly different notation) can be found in
       [M09].  The algorithm originates in [S72]

    REFERENCES:

    - [BMSS] Boston, Morain, Salvy, Schost, "Fast Algorithms for Isogenies."
    - [M09] Moody, "The Diffie-Hellman Problem and Generalization of Verheul's Theorem"
    - [S72] Stark, "Class-numbers of complex quadratic fields."

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks, compute_sequence_of_maps

        sage: E = EllipticCurve(GF(97), [1,0,1,1,0])
        sage: R.<x> = GF(97)[]; f = x^5 + 27*x^4 + 61*x^3 + 58*x^2 + 28*x + 21
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: (isom1, isom2, E1pr, E2pr, ker_poly) = compute_sequence_of_maps(E, E2, 11)
        sage: compute_isogeny_starks(E1pr, E2pr, 11)
        x^10 + 37*x^9 + 53*x^8 + 66*x^7 + 66*x^6 + 17*x^5 + 57*x^4 + 6*x^3 + 89*x^2 + 53*x + 8

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_starks(E, E2, 5)
        x^4 + 14*x^3 + x^2 + 34*x + 21
        sage: f**2
        x^4 + 14*x^3 + x^2 + 34*x + 21

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: R.<x> = QQ[]
        sage: f = x
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_starks(E, E2, 2)
        x

    """

    K = E1.base_field()
    R = PolynomialRing(K, 'x')
    x = R.gen()

    wp1 = E1.weierstrass_p(prec=4*ell+4)  #BMSS claim 2*ell is enough, but it is not M09
    wp2 = E2.weierstrass_p(prec=4*ell+4)

    # viewed them as power series in Z = z^2
    S = LaurentSeriesRing(K, 'Z')
    Z = S.gen()
    pe1 = 1/Z
    pe2 = 1/Z
    for i in xrange(2*ell+1):
        pe1 += wp1[2*i] * Z**i
        pe2 += wp2[2*i] * Z**i
    pe1 = pe1.add_bigoh(2*ell+2)
    pe2 = pe2.add_bigoh(2*ell+2)

    #print 'wps = ',pe1
    #print 'wps2 = ',pe2

    n = 1
    q = [R(1), R(0)]
    #p = [R(0), R(1)]
    T = pe2

    while ( q[n].degree() < (ell-1) ):
        #print 'n=', n

        n += 1
        a_n = 0
        r = -T.valuation()
        while (0 <= r):
            t_r = T[-r]
            #print '    r=',r
            #print '    t_r=',t_r
            #print '    T=',T
            a_n = a_n + t_r * x**r
            T = T - t_r*pe1**r
            r = -T.valuation()


        q_n = a_n*q[n-1] + q[n-2]
        q.append(q_n)
        #p_n = a_n*p[n-1] + q[n-2]
        #p.append(p_n)

        if (n == ell+1 or T == 0):
            if (T == 0 or T.valuation()<2):
                raise ValueError("The two curves are not linked by a cyclic normalized isogeny of degree %s" % ell)
            #print 'breaks here'
            break

        T = 1/T
        #print '  a_n=', a_n
        #print '  q_n=', q_n
        #print '  p_n=', p_n
        #print '  T = ', T

    qn = q[n]
    #pn= p[n]
    #print 'final  T = ', T
    #print '  f =', pn/qn

    qn = (1/qn.leading_coefficient())*qn
    #pn = (1/qn.leading_coefficient())*pn

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
        sage: ker_poly = compute_isogeny_starks(E, E2, 7); ker_poly
        x^6 + 2*x^5 + 20*x^4 + 11*x^3 + 36*x^2 + 35*x + 16
        sage: split_kernel_polynomial(E, ker_poly, 7)
        x^3 + x^2 + 28*x + 33

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
    Computes the kernel polynomial of the degree ``ell`` isogeny
    between ``E1`` and ``E2``.  There must be a degree ``ell``,
    cyclic, separable, normalized isogeny from ``E1`` to ``E2``.

    INPUT:

    - ``E1``        - an elliptic curve in short Weierstrass form.

    - ``E2``        - an elliptic curve in short Weierstrass form.

    - ``ell``       - the degree of the isogeny from ``E1`` to ``E2``.

    - ``algorithm`` - currently only ``starks`` (default) is implemented.

    OUTPUT:

    polynomial -- over the field of definition of ``E1``, ``E2``, that is the
                  kernel polynomial of the isogeny from ``E1`` to ``E2``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_kernel_polynomial

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: compute_isogeny_kernel_polynomial(E, E2, 5)
        x^2 + 7*x + 13
        sage: f
        x^2 + 7*x + 13

        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: compute_isogeny_kernel_polynomial(E, E2, 4)
        x^3 + x

    """

    ker_poly = compute_isogeny_starks(E1, E2, ell)
    ker_poly = split_kernel_polynomial(E1, ker_poly, ell)

    return ker_poly


def compute_intermediate_curves(E1, E2):
    r"""
    Computes isomorphism from ``E1`` to an intermediate domain and an
    isomorphism from an intermediate codomain to ``E2``.

    Intermediate domain and intermediate codomain, are in short
    Weierstrass form.

    This is used so we can compute `\wp` functions from the short
    Weierstrass model more easily.

    The underlying field must be of characteristic not equal to 2,3.

    INPUT:

    - ``E1`` - an elliptic curve
    - ``E2`` - an elliptic curve

    OUTPUT:

    tuple -- (``pre_isomorphism``, ``post_isomorphism``, ``intermediate_domain``,
              ``intermediate_codomain``):

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
    Given domain ``E1`` and codomain ``E2`` such that there is a
    degree ``ell`` separable normalized isogeny from ``E1`` to ``E2``,
    returns pre/post isomorphism, as well as intermediate domain and
    codomain, and kernel polynomial.

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
         x^2 - 61/3*x + 658/9)

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
         x^3 + x)

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
         x^5 + 67*x^4 + 13*x^3 + 35*x^2 + 77*x + 69)

    """

    (E1pr, E2pr, pre_isom, post_isom) = compute_intermediate_curves(E1, E2)

    ker_poly = compute_isogeny_kernel_polynomial(E1pr, E2pr, ell)

    return (pre_isom, post_isom, E1pr, E2pr, ker_poly)


##########################################################################
# The following section is all about computing l-isogenies, where l is
# a prime.  The genus 0 cases `l` = 2, 3, 5, 7 and 13 are
# implemented over any field of characteristic not 2, 3 or `l`; over
# `\QQ` the "sporadic" cases `l` = 11, 17, 19, 37, 43, 67 or 163 with
# only finitely many `j`-invariants each. are also implemented.
##########################################################################

@cached_function
def Fricke_polynomial(l):
    r"""
    Fricke polynomial for ``l`` =2,3,5,7,13.

    For these primes (and these only) the modular curve `X_0(l)` has
    genus zero, and its field is generated by a single modular
    function called the Fricke module (or Hauptmodul), `t`.  There is
    a classical choice of such a generator `t` in each case, and the
    `j`-function is a rational function of `t` of degree `l+1` of the
    form `P(t)/t` where `P` is a polynomial of degree `l+1`.  Up to
    scaling, `t` is determined by the condition that the ramification
    points above `j=\infty` are `t=0` (with ramification degree `1`)
    and `t=\infty` (with degree `l`).  The ramification above `j=0`
    and `j=1728` may be seen in the factorizations of `j(t)` and
    `k(t)` where `k=j-1728`.

    OUTPUT:

    The polynomial `P(t)` as an element of `\ZZ[t]`.

    TESTS::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import Fricke_polynomial
        sage: Fricke_polynomial(2)
        t^3 + 48*t^2 + 768*t + 4096
        sage: Fricke_polynomial(3)
        t^4 + 36*t^3 + 270*t^2 + 756*t + 729
        sage: Fricke_polynomial(5)
        t^6 + 30*t^5 + 315*t^4 + 1300*t^3 + 1575*t^2 + 750*t + 125
        sage: Fricke_polynomial(7)
        t^8 + 28*t^7 + 322*t^6 + 1904*t^5 + 5915*t^4 + 8624*t^3 + 4018*t^2 + 748*t + 49
        sage: Fricke_polynomial(13)
        t^14 + 26*t^13 + 325*t^12 + 2548*t^11 + 13832*t^10 + 54340*t^9 + 157118*t^8 + 333580*t^7 + 509366*t^6 + 534820*t^5 + 354536*t^4 + 124852*t^3 + 15145*t^2 + 746*t + 13
    """
    Zt = PolynomialRing(ZZ,'t')
    t = Zt.gen()
    if l==2: return (t+16)**3
    elif l==3: return (t+3)**3*(t+27)
    elif l==5: return (t**2+10*t+5)**3
    elif l==7: return (t**2+5*t+1)**3 * (t**2+13*t+49)
    elif l==13: return (t**2+5*t+13)*(t**4+7*t**3+20*t**2+19*t+1)**3
    else:
        raise ValueError, "The only genus zero primes are 2, 3, 5, 7 or 13."

@cached_function
def Fricke_module(l):
    r"""
    Fricke module for ``l`` =2,3,5,7,13.

    For these primes (and these only) the modular curve `X_0(l)` has
    genus zero, and its field is generated by a single modular
    function called the Fricke module (or Hauptmodul), `t`.  There is
    a classical choice of such a generator `t` in each case, and the
    `j`-function is a rational function of `t` of degree `l+1` of the
    form `P(t)/t` where `P` is a polynomial of degree `l+1`.  Up to
    scaling, `t` is determined by the condition that the ramification
    points above `j=\infty` are `t=0` (with ramification degree `1`)
    and `t=\infty` (with degree `l`).  The ramification above `j=0`
    and `j=1728` may be seen in the factorizations of `j(t)` and
    `k(t)` where `k=j-1728`.

    OUTPUT:

    The rational function `P(t)/t`.

    TESTS::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import Fricke_module
        sage: Fricke_module(2)
        (t^3 + 48*t^2 + 768*t + 4096)/t
        sage: Fricke_module(3)
        (t^4 + 36*t^3 + 270*t^2 + 756*t + 729)/t
        sage: Fricke_module(5)
        (t^6 + 30*t^5 + 315*t^4 + 1300*t^3 + 1575*t^2 + 750*t + 125)/t
        sage: Fricke_module(7)
        (t^8 + 28*t^7 + 322*t^6 + 1904*t^5 + 5915*t^4 + 8624*t^3 + 4018*t^2 + 748*t + 49)/t
        sage: Fricke_module(13)
        (t^14 + 26*t^13 + 325*t^12 + 2548*t^11 + 13832*t^10 + 54340*t^9 + 157118*t^8 + 333580*t^7 + 509366*t^6 + 534820*t^5 + 354536*t^4 + 124852*t^3 + 15145*t^2 + 746*t + 13)/t
    """
    try:
        t = PolynomialRing(QQ,'t').gen()
        return Fricke_polynomial(l) / t
    except ValueError:
        raise ValueError, "The only genus zero primes are 2, 3, 5, 7 or 13."

@cached_function
def Psi(l, use_stored=True):
    r"""
    Generic kernel polynomial for genus one primes.

    For each of the primes `l` for which `X_0(l)` has genus zero
    (namely `l=2,3,5,7,13`), we may define an elliptic curve `E_t`
    over `\QQ(t)`, with coefficients in `\ZZ[t]`, which has good
    reduction except at `t=0` and `t=\infty` (which lie above
    `j=\infty`) and at certain other values of `t` above `j=0` when
    `l=3` (one value) or `l\equiv1\pmod{3}` (two values) and above
    `j=1728` when `l=2` (one value) or `l\equiv1 \pmod{4}` (two
    values).  (These exceptional values correspond to endomorphisms of
    `E_t` of degree `l`.)  The `l`-division polynomial of `E_t` has a
    unique factor of degree `(l-1)/2` (or 1 when `l=2`), with
    coefficients in `\ZZ[t]`, which we call the Generic Kernel
    Polynomial for `l`.  These are used, by specialising `t`, in the
    function :meth:`isogenies_prime_degree_genus_0`, which also has to
    take into account the twisting factor between `E_t` for a specific
    value of `t` and the short Weierstrass form of an elliptic curve
    with `j`-invariant `j(t)`.  This enables the computation of the
    kernel polynomials of isogenies without having to compute and
    factor division polynomials.

    All of this data is quickly computed from the Fricke modules, except
    that for `l=13` the factorization of the Generic Division Polynomial
    takes a long time, so the value have been precomputed and cached; by
    default the cached values are used, but the code here will recompute
    them when ``use_stored`` is ``False``, as in the doctests.

    INPUT:

    - ``l`` -- either 2, 3, 5, 7, or 13.

    - ``use_stored`` (boolean, default True) -- If True, use
      precomputed values, otherwise compute them on the fly.

    .. note::

       This computation takes a negligible time for `l=2,3,5,7'
       but more than 100s for `l=13`.  The reason
       for allowing dynamic computation here instead of just using
       precomputed values is for testing.

    TESTS::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import Fricke_module, Psi
        sage: assert Psi(2, use_stored=True) == Psi(2, use_stored=False)
        sage: assert Psi(3, use_stored=True) == Psi(3, use_stored=False)
        sage: assert Psi(5, use_stored=True) == Psi(5, use_stored=False)
        sage: assert Psi(7, use_stored=True) == Psi(7, use_stored=False)
        sage: assert Psi(13, use_stored=True) == Psi(13, use_stored=False) # not tested (very long time)
    """
    if not l in [2,3,5,7,13]:
        raise ValueError, "Genus zero primes are 2, 3, 5, 7 or 13."

    R = PolynomialRing(ZZ,2,'Xt')
    X,t = R.gens()

    if use_stored:
        if l==2:
            return X + t + 64
        if l==3:
            return X + t + 27
        if l==5:
            return X**2 + 2*X*(t**2 + 22*t + 125)+ (t**2 + 22*t + 89) * (t**2 + 22*t + 125)
        if l==7:
            return (X**3 + 3*(t**2 + 13*t + 49)*X**2
                    + 3*(t**2 + 13*t + 33)*(t**2 + 13*t + 49)*X
                    + (t**2 + 13*t + 49)*(t**4 + 26*t**3 + 219*t**2 + 778*t + 881))
        if l==13:
            return (t**24 + 66*t**23 + 2091*t**22 + 6*X*t**20 + 42582*t**21 + 330*X*t**19 + 627603*t**20 + 8700*X*t**18 + 7134744*t**19 + 15*X**2*t**16 + 146886*X*t**17 + 65042724*t**18 + 660*X**2*t**15 + 1784532*X*t**16 + 487778988*t**17 + 13890*X**2*t**14 + 16594230*X*t**15 + 3061861065*t**16 + 20*X**3*t**12 + 186024*X**2*t**13 + 122552328*X*t**14 + 16280123754*t**15 + 660*X**3*t**11 + 1774887*X**2*t**12 + 735836862*X*t**13 + 73911331425*t**14 + 10380*X**3*t**10 + 12787272*X**2*t**11 + 3646188342*X*t**12 + 287938949178*t**13 + 15*X**4*t**8 + 102576*X**3*t**9 + 71909658*X**2*t**10 + 15047141292*X*t**11 + 964903805434*t**12 + 330*X**4*t**7 + 707604*X**3*t**8 + 321704316*X**2*t**9 + 51955096824*X*t**10 + 2781843718722*t**11 + 3435*X**4*t**6 + 3582876*X**3*t**7 + 1155971196*X**2*t**8 + 150205315932*X*t**9 + 6885805359741*t**10 + 6*X**5*t**4 + 21714*X**4*t**5 + 13632168*X**3*t**6 + 3343499244*X**2*t**7 + 362526695094*X*t**8 + 14569390179114*t**9 + 66*X**5*t**3 + 90660*X**4*t**4 + 39215388*X**3*t**5 + 7747596090*X**2*t**6 + 725403501318*X*t**7 + 26165223178293*t**8 + 336*X**5*t**2 + 255090*X**4*t**3 + 84525732*X**3*t**4 + 14206132008*X**2*t**5 + 1189398495432*X*t**6 + 39474479008356*t**7 + X**6 + 858*X**5*t + 472143*X**4*t**2 + 132886992*X**3*t**3 + 20157510639*X**2*t**4 + 1569568001646*X*t**5 + 49303015587132*t**6 + 1014*X**5 + 525954*X**4*t + 144222780*X**3*t**2 + 21320908440*X**2*t**3 + 1622460290100*X*t**4 + 49941619724976*t**5 + 272259*X**4 + 96482100*X**3*t + 15765293778*X**2*t**2 + 1260038295438*X*t**3 + 39836631701295*t**4 + 29641924*X**3 + 7210949460*X**2*t + 686651250012*X*t**2 + 23947528862166*t**3 + 1506392823*X**2 + 231462513906*X*t + 10114876838391*t**2 + 35655266790*X + 2644809206442*t + 317295487717)
# The coefficients for l=13 are:
# X**6: 1
# X**5: (6) * (t**2 + 5*t + 13) * (t**2 + 6*t + 13)
# X**4: (3) * (t**2 + 5*t + 13) * (t**2 + 6*t + 13) * (5*t**4 + 55*t**3 + 260*t**2 + 583*t + 537)
# X**3: (4) * (t**2 + 5*t + 13) * (t**2 + 6*t + 13)**2 * (5*t**6 + 80*t**5 + 560*t**4 + 2214*t**3 + 5128*t**2 + 6568*t + 3373)
# X**2: (3) * (t**2 + 5*t + 13)**2 * (t**2 + 6*t + 13)**2 * (5*t**8 + 110*t**7 + 1045*t**6 + 5798*t**5 + 20508*t**4 + 47134*t**3 + 67685*t**2 + 54406*t + 17581)
# X**1: (6) * (t**2 + 5*t + 13)**2 * (t**2 + 6*t + 13)**3 * (t**10 + 27*t**9 + 316*t**8 + 2225*t**7 + 10463*t**6 + 34232*t**5 + 78299*t**4 + 122305*t**3 + 122892*t**2 + 69427*t + 16005)
# X**0: (t**2 + 5*t + 13)**2 * (t**2 + 6*t + 13)**3 * (t**14 + 38*t**13 + 649*t**12 + 6844*t**11 + 50216*t**10 + 271612*t**9 + 1115174*t**8 + 3520132*t**7 + 8549270*t**6 + 15812476*t**5 + 21764840*t**4 + 21384124*t**3 + 13952929*t**2 + 5282630*t + 854569)
#

    # Here the generic kernel polynomials are actually calculated:
    j = Fricke_module(l)
    k = j-1728
    from sage.misc.all import prod
    f = prod( [p for p,e in j.factor() if e==3]
             +[p for p,e in k.factor() if e==2])
    A4 = -3*t**2*j*k // f**2
    A6 = -2*t**3*j*k**2 // f**3
    E = EllipticCurve([0,0,0,A4,A6])
    assert E.j_invariant() == j
    return E.division_polynomial(l,X).factor()[0][0]


def isogenies_prime_degree_genus_0(E, l=None):
    """
    Returns list of ``l`` -isogenies with domain ``E``.

    INPUT:

    - ``E`` -- an elliptic curve.

    - ``l`` -- either None or 2, 3, 5, 7, or 13.

    OUTPUT:

    (list) When ``l`` is None a list of all isogenies of degree 2, 3,
    5, 7 and 13, otherwise a list of isogenies of the given degree.

    .. note::

       This function would normally be invoked indirectly via
       ``E.isogenies_prime_degree(l)``, which automatically calls the
       appropriate function.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_prime_degree_genus_0
        sage: E = EllipticCurve([0,12])
        sage: isogenies_prime_degree_genus_0(E, 5)
        []

        sage: E = EllipticCurve('1450c1')
        sage: isogenies_prime_degree_genus_0(E)
        [Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 300*x - 1000 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 5950*x - 182250 over Rational Field]

        sage: E = EllipticCurve('50a1')
        sage: isogenies_prime_degree_genus_0(E)
        [Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + y = x^3 - x - 2 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 - 126*x - 552 over Rational Field,
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + x*y + y = x^3 - x - 2 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 - 76*x + 298 over Rational Field]
    """
    if not l in [2, 3, 5, 7, 13, None]:
        raise ValueError, "Isogenies must be of degree 2, 3, 5, 7 or 13."
    F = E.base_ring()
    j = E.j_invariant()
    if F.characteristic() in [2, 3, l]:
        raise NotImplementedError, "2, 3, 5, 7 and 13-isogenies are not yet implemented in characteristic 2 and 3, and when the characteristic is the same as the degree of the isogeny."
    if l==2:
        return isogenies_2(E)
    if l==3:
        return isogenies_3(E)
    if j==F(0):
        if l==5:
            return isogenies_5_0(E)
        if l==7:
            return isogenies_7_0(E)
        if l==13:
            return isogenies_13_0(E)
    if j==F(1728):
        if l==5:
            return isogenies_5_1728(E)
        if l==7:
            return isogenies_7_1728(E)
        if l==13:
            return isogenies_13_1728(E)

    if l != None:
        R = PolynomialRing(F,'t')
        t = R.gen()
        f = R(Fricke_polynomial(l))
        t_list = (f-j*t).roots(multiplicities=False)
        t_list.sort()
        # The generic kernel polynomial applies to a standard curve
        # E_t with the correct j-invariant; we must compute the
        # appropriate twising factor to scale X by:
        c4, c6 = E.c_invariants()
        T = c4/(3*c6)
        jt = Fricke_module(l)
        kt = jt-1728
        from sage.misc.all import prod
        psi = Psi(l)
        X = t
        f = R(prod( [p for p,e in jt.factor() if e==3]
                 +[p for p,e in kt.factor() if e==2]))
        kernels = [R(psi(X*T*(j-1728)*t0/f(t0),t0)) for t0 in t_list]
        kernels = [ker.monic() for ker in kernels]
        E1 = EllipticCurve([-27*c4,-54*c6])
        w = E.isomorphism_to(E1)
        model = "minimal" if F is QQ else None
        isogs = [E1.isogeny(kernel=ker, model=model) for ker in kernels]
        [isog.set_pre_isomorphism(w) for isog in isogs]
        return isogs

    if l == None:
        return sum([isogenies_prime_degree_genus_0(E, l) for l in [2,3,5,7,13]],[])


# The following is data to be used in isogenies_sporadic_Q. Over Q
# there are only finitely many j-invariants of curves with l-isogenies
# where l is not equal to 2, 3, 5, 7 or 13. In these cases l is equal
# to 11, 17, 19, 37, 43, 67 or 163. We refer to these l as "sporadic".
#
# isog_table is indexed by the possible sporadic (l,j) pairs and lists
# ((a4,a6),f) where (a4,a6) are the invariants of a short Weierstrass
# model of one curve with that j-invariant,and f is the factor of
# degree (l-1)/2 of the l-division polynomial of that curve. Then
# whenever we have a curve of that j-invariant, we can compute the
# corresponding l-isogeny by just scaling f by the right twisting
# factor and using the result as a kernel-polynomial.

# Sample code to precompute this data.  Note that working over ZZ
# instead of QQ is a lot faster (in Sage 4.2):
#
#sage: j = -297756989/2; l = 17
#sage: E = EllipticCurve(ZZ,EllipticCurve(j=j).short_weierstrass_model().ainvs())
#sage: E.division_polynomial(l,polygen(ZZ)).factor()[1][0].coefficients()
#
#sage: j = -884736000; l = 43
#sage: E = EllipticCurve(ZZ,EllipticCurve(j=j).short_weierstrass_model().ainvs())
#sage: E.division_polynomial(l,polygen(ZZ)).factor()[1][0].coefficients()


isog_table = dict()

isog_table[(11,QQ('-32768'))] = ([-9504, 365904], [1294672896, -92835072, 1463616, 7920, -264, 1])

isog_table[(11,QQ('-121'))] = ([-3267, -280962], [1480352841, -56169531, -2829222, 10890, 429, 1])

isog_table[(11,QQ('-24729001'))] = ([-38907, -2953962], [-20349931239, -424530315, -134838, 53658, 429, 1])

isog_table[(17,QQ('-297756989/2'))] = ([-3940515, 3010787550], [-6458213126940667330314375, 34699336325466068070000, -72461450055340471500, 68342601718080000, -15140380554450, -25802960400, 23981220, -8160, 1])

isog_table[(17,QQ('-882216989/131072'))] = ([-856035, -341748450], [103687510635057329105625, 961598491955315190000, 1054634146768300500, -6553122389064000, -14554350284850, -2046589200, 13185540, 8160, 1])

isog_table[(19,QQ('-884736'))] = ([-608, 5776], [-34162868224, -8540717056, 6405537792, -1123778560, 84283392, -2033152, -92416, 6992, -152, 1])

isog_table[(37,QQ('-9317'))] = ([-10395, 444150],[-38324677699334121599624973029296875, -17868327793500376961572310472656250, 2569568362004197901139023084765625, -95128267987528547588017818750000, -822168183291347061312510937500, 134395594560592096297190625000, -2881389756919344324888937500, -2503855007083401977250000, 922779077075655997443750, -11503912310262102937500, -18237870962450291250, 1457822151548910000, -10087015556047500, -13677678063000, 490243338900, -2461460400, 5198445, -4410, 1])

isog_table[(37,QQ('-162677523113838677'))] = ([-269675595, -1704553285050], [-31653754873248632711650187487655160190139073510876609346911928661154296875/37, -1469048260972089939455942042937882262144594798448952781325533511718750, -1171741935131505774747142644126089902595908234671576131857702734375, -574934780393177024547076427530739751753985644656221274606250000, -193516922725803688001809624711400287605136013195315374687500, -47085563820928456130325308223963045033502182349693125000, -8472233937388712980597845725196873697064639957437500, -1124815211213953261752081095348112305023653750000, -105684015609077608033913080859605951322531250, -5911406027236569746089675554748135312500, 22343907270397352965399097794968750, 43602171843758666292581116410000, 5054350766002463251474186500, 350135768194635636171000, 16633063574896677300, 549939627039600, 12182993865, 163170, 1])

isog_table[(43,QQ('-884736000'))] = ([-13760, 621264],[-1961864562041980324821547425314935668736, 784270445793223959453256359333693751296, -120528107728500223255333768387027271680, 10335626145581464192664472924270362624, -568426570575654606865505142156820480, 21261993723422650574629752537088000, -544630471727787626557612832587776, 8870521306520473088172555763712, -54993059067301585878494740480, -1434261324709904840432549888, 50978938193065926383894528, -845761855773797582372864, 8627493611216601088000, -48299605284169187328, -32782260293713920, 3415534989828096, -34580115625984, 199359712512, -730488128, 1658080, -2064, 1])

isog_table[(67,QQ('-147197952000'))] = ([-117920, 15585808], [426552448394636714720553816389274308035895411389805883034985546818882031845376, -55876556222880738651382959148329502876096075327084935039031884373558741172224, 3393295715290183821010552313572221545212247684503012173117764703828786020352, -125729166452196578653551230178028570067747190427221869867485520072257044480, 3121342502030777257351089270834971957072933779704445667351054593298530304, -52544031605544530265465344472543470442324636919759253720520768014516224, 532110915869155495738137756847596184665209453108323879594125221167104, -399031158106622651277981701966309467713625045637309782055519780864, -101914346170769215732007802723651742508893380955930030421292613632, 2296526155500449624398016447877283594461904009374321659789443072, -31950871094301541469458501953701002806003991982768349794795520, 329792235011603804948028315065667439678526339671142107709440, -2655636715955021784085217734679612378726442691190553837568, 16825164648840434987220620681420687654501026066872664064, -81705027839007003131400500185224450729843244954288128, 273656504606483403474090105104132405333665144373248, -320807702482945680116212224172370503903312084992, -3166683390779345463318656135338172047199043584, 27871349428383710305216046431806697565585408, -132774697798318602604125735604528772808704, 436096215568182871014215818309741314048, -964687143341252402362763535357837312, 942144169187362941776488535425024, 2794850106281773765892648206336, -17236916236678037389276086272, 50979778712911923486851072, -105035658611718440992768, 161833913559276412928, -188675698610077696, 163929317513984, -102098677888, 42387952, -10184, 1])

isog_table[(163,QQ('-262537412640768000'))] = ([-34790720, 78984748304], [11937033609701771608400008302859181021419097217821452062564326844817505228663098831008549807845273606080544999996020379100442085010011180910986311230744539492209201267851475053850145465949031314131282944262561686476431332761205504055553653209756637748416809823347727611717499380045401030656, -195130071447306655273524569473852352290010571088570084255929596426827972913675959924834214553822436039185267072655230231264389277855757579340840680919752147359277865429131554516271296434318501723193704667609983954393732283396979741770275596479169906167026623493287956522391758781645062144, 1427502734760154784443139998885173873805698542265553947708985997014264802221336771918193375050938963016560184673026895611200993837000754187703921657388546609660768979949394789838000919808563779734300559653270446650866261801622261730386746625133970663551831816512698564779980290036596736, -5599207432852992026535855102327743354503965526118767584653806869099469434869389394471452436350209911401734358692516268244363981178003440912521398279843764610451103952982463228121628686179444592844803463454880601969041617668281227800515971090274318968755339891485627349912827843313664, 7275289443663173375117187980461451723630939807433676133488193878209459535146604043302393215325736700753114238886326617501903060793779318898498537340122928915051396994852249905774785383484964157718761647679761326405167986820403416875550489904598260809039588027484880575324789669888, 53852767252251176711316610474871670511285949778037067647446607721301053062542360857644949076379793988972999080753585666268471181395581176311086709516236139859119228114940726664020344711163961450649040900179940667704933096423082541578086102752168151408016128448757608609516355584, -441024601763314150109608522830440149169771622404264936395215000306188259683627792348927267455527449622616337194882061278315837263834842622422058957336768392927319786456429412773335490683756389242358011775192147290260088212928183893496023497385031567935590156495818971603271680, 1923729652470036657640234466480856403679261373466564217770510758620007251234737653848912112257887106718711179570266106782165127910894792855741763009046440188368346201617134183692667813264354753020558715833289878752200363206184231817927249353441775148372373391040829093576704, -6174291905252822700639099970860178217904525535465573253756653994759145320917644436646210500761314537511653465341404223911052812759676375813287317189557074876672717009913962857515263325087926781282727274494599622098970326546685474959128901831166082055481492871737587531776, 15963605433939241348215386312076502855685142708851007052060264896693306401589916424547272907372369971200908960298559196028000249867821387055883338894215055448441454208145184664251082508258481820123154014612667230187589409024581605410594578299015578613903351914203971584, -34687777075075826936035422448136670540047824478867437603900549472003349053758981331083701038247602826434094223697740838012005963350495346494192086319952577715870471857062004211165110944749546489357024333327851492892063261902845789001287880518909468079710544890691584, 64935614565701660925102862834870557968917923758048315817978089985628103203452108238400473842703675144415856702414743965076964752315718880263247146806212144113325233966196430978981658754610129252819800656290540439882400228191972297000590950499022332532802841477120, -106450113192928986500166510746460290010686605565802553549898115364953461571938149159855780359758326598170956839310849727549767334328481520722457671184262149140997347199434067316119593873464084542168958564857207310255104210252802288265107074620056433232296542208, 154588468341182819097156049664780175805765082709542906853062319968823210709884467518965950031344875777123989196305422245441384128571495551803479155842444104959551491130118809781140758556182128829570214830557657320822341658937375265326185627024588353675722752, -200559814675071715715147513664537685696543665835231560123119674695431140901874837885144270078673500682009909076723885595600964688749788075855566556738418513204701882392037282897429720984011892144438563658741417367819609826485149957684720592798005782904832, 233913230287794067984901680255098952718268709897211662836962202594844321932231309396660409956864579218254273924904394011519168528453166798219739516110916032144001430976170021740673624911101610999174313953331179592523559209806130064362331035786304028672, -246353703535830799051965068435640665615241585380168769205605419959931340172791084855084136413845436489354260626697043416419953596294256639939636661074827002386688825389366314922889935515780006872313799807871034941783710991661439846743406825846079488, 234982742638283161239147075499923192184558615040750093346909554959954534765162320677404212311317141015056956519607730625264033505832417237335526289717986324622989530401755659861160278474288472483456460761582961131223878005034852988312569782468608, -203278213565303083958289848507904983111914464281937008996800600743049033135780053019308180654889993715847194782810638377233384956677984229396984592095978852141937402724743415866647478390224972086303117541277057454016741537801255285401042550784, 159421011266364165773350759598961923468646295930442671951663830439105229019095602850806206849694626149081055860876320252885592409781300283802068036284087012054926328984117466238728138017296876920468741929209733382315316048593848804000661504, -113020918346823224387675962727710849438448560054728764538889949910907956081636084818673918939857811981198128625134001275445247927733612860826192033364144709296859805908448480126055187856540001252440723696556907592809391568682658983575552, 71945777186858148520189111385851657591994785437936790663752048596268908397761761228530149720647914698227096033298596468433406211022425872184679773706985835273196359946527329295352840585601589987654250252509940697181497393634479177728, -40553506248279853197899380016510447069205875473858373944916667311741237594752103405195544436215283638784303466518029972315361472909380024525230365389883507930702477745185560842399992120410631726800835247081534888734998176190169088, 19639696188973552762819015734509385187497753279018359326954702329472344393424778261565226791972589248949076879073752940457430177635689762236418319921157145929450429275277598105096782299310906043071624518383647513981402864418816, -7557145635717782395604678688732309739123765964398036300935902895029921544449817739871127195837816246027098669576710973724603305169613035210748794125702623638061784934800621663541549758981147547220745886014319375033833619456, 1657074863449131529294393850366914004166071540911557375298352670995793916631379837217699843183161870412953109096989706473414891817476345269315000071815674837414249807862563049048940888990780356828116926372826050377809920, 594537607406738614488647166334439665406200193881374352383582637957395810395355986219181323336507140663007282377906531198753011484726171889166159296232442266402815351908351664480142997422997116094990105400004717838336, -1063381335685557601789471169771716160614658527547206527560235661204306306147122486589145738731952820149525220010993135210540125800587297277176980780550749659568216404307556497848427558893515024809067113577642459136, 866504016710397283549081797838552191017551170056714418172658679484342314216337494688411296133019156758665490378562995489461120013841310945167448587950589192882728019730184100760635093600626150305226434294054912, -544456451437710104677464982275658569122639544706468906698603670935572135268235880234260773620839923609980624261711876187191086008709319550421652238153751371247896403810296918434808016324776218527300472274944, 291391938923983436973777049837536733140215048414562274102609644110225500962890424001140034799254161420207705596799832439286323749581485058440291289292376834767289200150488260103720502388531610043446984704, -137669336373928330264003425607593379960726875901008061728935571651808147336578784033262243606118618305155146994101199276509296485578362621442124562716722538959248783675765430495753630524650420015464448, 58305039769563586190846386880024028458953483779423358562963969531061253550565868881984384427294827940260948562007946246585988242539388915757370599713569828414708265704775708180369817440224330907648, -22263674613297409152475382098514675408629524088165884591926938276387287183647211206514841583522457779552325379611699721363009213935995055122374833256328128159685068070573439567341889102047870976, 7659525572873653643686667682470156611730287715779213959462404612107764289631014622079262434509991525269350939684751076168111178931547830737285816290499346814316593688339716060537705476390912, -2356369509601030092354164626593105564281508081035914666378162646456624099153984023946457193322861690246151591627985138670141283491665875837349056368982626874320059053912723460560567402496, 636306064623878566247238117586667093960745761448319133509756604897082306674457051834692637631208787669480579961449146859307914219525683813169045631019945325543755784485150162136596480, -144354698667964959131302121388611435649331689922451077357726958728526160763985021377919777057073652734054824995380136026354322489788621439866623621032883911971992887755746217820160, 24121955586756417060840806312875152536841537469465641066665477185642959713399102667004561835703554642486033432650174169972215286732914841507194929969333501923213681889898397696, -1062982554332400020436047718307967606435230424325029173087078414691100152298424407470637592479633445244263771003475924500470457802708388802132021409732119617604143984672768, -1289959173517591900879403877417447464323361120695839905879179383576917611515563011791770601288187258710836960962551744314331281738572022125094897455138775637822772609024, 713094034093241133139975542207369737377005379966568397923649557581354608722973890056780743426440880397949894405054289149004564424971813787111450018625163532324306944, -254627379967901645020979419858451704766383325190491161773871149760504423262065255127251316663797380198173202687716259997416292667070368594859769612820620832145408, 73141763506085097918144459824978126154049420184880822723900482544151820781804249152002290712777719695036970656528265835764223704583688284269989955443519127552, -17838537225679516940483486695156689114251551141417084228356667456567223074797422210814066353881812352204911178711536306766419580650031647883778868599848960, 3732949472404613028683276858400771022649020140479952475050562287824228475488013407594597423014344735572239887680502683815016952233794807477380489150464, -658090831088275932021697572468169724648153295038942023393017242911882324988827102897240587742441233130939688197251981440661637431663180512231424000, 91044773494274344573228085560493548892450695165631132364575161840790371925870648417274049055693061297980952723608086693629684733411639066361856, -7206875978668204173369106020317806614830036116031258178074085756310402308591937698714169573782276745187407867676712962154996289706219536384, -819751497058465540849522914081577930969395072622202448379227754265195552229316711754131018363430793506957620878700584807196119659446272, 524492219141524596329431820120551446539246304466694008535436881375050446937946483408452355745347329681867559263476448413099237572608, -144831995016562693238485974645571530330971707229081912571899541756458417275452663670185243929333997861822834344607487043313860608, 29914869637979316056234691641068533246353548074478964560456384907200155525359490145390698951538414992573246875654471802683392, -5051703199007242317980870828294675665554938059871421101373455640671888804716466166733214891733545898164243959713215021056, 704266224450794229753813389400889235502895246267920454892105495898823500254340026026030392861415807391779452916596736, -77090914220067457257437168730743598007998214624866646138957596697882331770536104168830311427039146990575867133952, 5185529715527985800186976473360862373518426289807725150611458977619779754244941008386724289180811582624497664, 253997773097193639890870939442154039645526260966991421810387625142025920490445698037959549218977356972032, -167663051891371188948811201282888303918118872228509878481620659513761819092806840088343319460416847872, 37111139278046773379558801627698801653586358926939257719148812112681156722739468064234356797865984, -6094794232647488071603134809336879629321484740792083961518604911291776122012806633249006682112, 837291566508665017874331015912241564512130987193765470835615200761060354034220082158108672, -100296887606046826616864046733353580911852480578758327196473644505408067503452896362496, 10671531126308719182671886696772003009468978048122050876385337250105289881860177920, -1017660211057522997139542974021461951590281854561469513981031668621068390629376, 87341651749945752980300268787472991750253192293865191008561736741703122944, -6755055811837604655292473395831392832731215030180150346004905802596352, 470435697396489049890276884918878813522261509046808135502735081472, -29431936259973770189795237369002643735145492537132570738950144, 1647955729708147681608319291051493831298410618099083509760, -82151504501901440318090519089275042334656836617109504, 3621394962905672351516373383975280471835697741824, -139946424433346723333934826799658073276284928, 4689146702972667351768128807176661303296, -134326059738838559898371882636836864, 3230090486862617839209120497664, -63628315178835325152092160, 992573003659549935360, -11668560578458368, 95539651616, -472048, 1])



def isogenies_sporadic_Q(E, l=None):
    """
    Returns list of ``l`` -isogenies with domain ``E`` (defined over `\QQ`).

    Returns a list of sporadic l-isogenies from E (l = 11, 17, 19, 37,
    43, 67 or 163). Only for elliptic curves over `\QQ`.

    INPUT:

    - ``E`` -- an elliptic curve defined over `\QQ`.

    - ``l`` -- either None or a prime number.

    OUTPUT:

    (list) If ``l`` is None, a list of all isogenies with domain ``E``
    and of degree 11, 17, 19, 37, 43, 67 or 163; otherwise a list of
    isogenies of the given degree.

    .. note::

       This function would normally be invoked indirectly via
       ``E.isogenies_prime_degree(l)``, which automatically calls the appropriate
       function.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_sporadic_Q
        sage: E = EllipticCurve('121a1')
        sage: isogenies_sporadic_Q(E, 11)
        [Isogeny of degree 11 from Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 30*x - 76 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 305*x + 7888 over Rational Field]
        sage: isogenies_sporadic_Q(E, 13)
        []
        sage: isogenies_sporadic_Q(E, 17)
        []
        sage: isogenies_sporadic_Q(E)
        [Isogeny of degree 11 from Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 30*x - 76 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 305*x + 7888 over Rational Field]

        sage: E = EllipticCurve([1, 1, 0, -660, -7600])
        sage: isogenies_sporadic_Q(E, 17)
        [Isogeny of degree 17 from Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 660*x - 7600 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 878710*x + 316677750 over Rational Field]
        sage: isogenies_sporadic_Q(E)
        [Isogeny of degree 17 from Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 660*x - 7600 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 878710*x + 316677750 over Rational Field]
        sage: isogenies_sporadic_Q(E, 11)
        []

        sage: E = EllipticCurve([0, 0, 1, -1862, -30956])
        sage: isogenies_sporadic_Q(E, 11)
        []
        sage: isogenies_sporadic_Q(E, 19)
        [Isogeny of degree 19 from Elliptic Curve defined by y^2 + y = x^3 - 1862*x - 30956 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - 672182*x + 212325489 over Rational Field]
        sage: isogenies_sporadic_Q(E)
        [Isogeny of degree 19 from Elliptic Curve defined by y^2 + y = x^3 - 1862*x - 30956 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - 672182*x + 212325489 over Rational Field]

        sage: E = EllipticCurve([0, -1, 0, -6288, 211072])
        sage: E.conductor()
        19600
        sage: isogenies_sporadic_Q(E,37)
        [Isogeny of degree 37 from Elliptic Curve defined by y^2 = x^3 - x^2 - 6288*x + 211072 over Rational Field to Elliptic Curve defined by y^2 = x^3 - x^2 - 163137088*x - 801950801728 over Rational Field]

        sage: E = EllipticCurve([1, 1, 0, -25178045, 48616918750])
        sage: E.conductor()
        148225
        sage: isogenies_sporadic_Q(E,37)
        [Isogeny of degree 37 from Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 25178045*x + 48616918750 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 970*x - 13075 over Rational Field]

        sage: E = EllipticCurve([-3440, 77658])
        sage: E.conductor()
        118336
        sage: isogenies_sporadic_Q(E,43)
        [Isogeny of degree 43 from Elliptic Curve defined by y^2 = x^3 - 3440*x + 77658 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 6360560*x - 6174354606 over Rational Field]

        sage: E = EllipticCurve([-29480, -1948226])
        sage: E.conductor()
        287296
        sage: isogenies_sporadic_Q(E,67)
        [Isogeny of degree 67 from Elliptic Curve defined by y^2 = x^3 - 29480*x - 1948226 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 132335720*x + 585954296438 over Rational Field]

        sage: E = EllipticCurve([-34790720, -78984748304])
        sage: E.conductor()
        425104
        sage: isogenies_sporadic_Q(E,163)
        [Isogeny of degree 163 from Elliptic Curve defined by y^2 = x^3 - 34790720*x - 78984748304 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 924354639680*x + 342062961763303088 over Rational Field]
    """
    if E.base_ring() != QQ:
        raise ValueError, "The elliptic curve must be defined over QQ."
    j = E.j_invariant()
    j = QQ(j)
    if l != None:
        if not l.is_prime():
            raise ValueError("%s is not prime."%l)
        if not isog_table.has_key((l,j)):
            return []
        Ew = E.short_weierstrass_model()
        E_to_Ew = E.isomorphism_to(Ew)
        c4, c6 = Ew.c_invariants()
        (a4,a6), f = isog_table[(l,j)]
        d = (c6*a4)/(18*c4*a6) # twisting factor
        R = PolynomialRing(E.base_field(),'X')
        n = len(f)
        ker = R([d**(n-i-1) * f[i] for i in range(n)])
        isog = Ew.isogeny(kernel=ker, degree=l, model="minimal", check=False)
        isog.set_pre_isomorphism(E_to_Ew)
        return [isog]
    if l is None:
        if isog_table.has_key((11,j)):
            return isogenies_sporadic_Q(E, Integer(11))
        if isog_table.has_key((17,j)):
            return isogenies_sporadic_Q(E, Integer(17))
        if isog_table.has_key((19,j)):
            return isogenies_sporadic_Q(E, Integer(19))
        if isog_table.has_key((37,j)):
            return isogenies_sporadic_Q(E, Integer(37))
        if isog_table.has_key((43,j)):
            return isogenies_sporadic_Q(E, Integer(43))
        if isog_table.has_key((67,j)):
            return isogenies_sporadic_Q(E, Integer(67))
        if isog_table.has_key((163,j)):
            return isogenies_sporadic_Q(E, Integer(163))
        else:
            return []

def isogenies_2(E):
    """
    Returns a list of all 2-isogenies with domain ``E``.

    INPUT:

    - ``E`` -- an elliptic curve.

    OUTPUT:

    (list) 2-isogenies with domain ``E``.  In general these are
    normalised, but over `\QQ` the codomain is a minimal model.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_2
        sage: E = EllipticCurve('14a1'); E
        Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field
        sage: [phi.codomain().ainvs() for phi in isogenies_2(E)]
        [(1, 0, 1, -36, -70)]

        sage: E = EllipticCurve([1,2,3,4,5]); E
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field
        sage: [phi.codomain().ainvs() for phi in isogenies_2(E)]
        []
        sage: E = EllipticCurve(QQbar, [9,8]); E
        Elliptic Curve defined by y^2 = x^3 + 9*x + 8 over Algebraic Field
        sage: isogenies_2(E) # not implemented
    """
    f2 = E.division_polynomial(2)
    x2 = f2.roots(multiplicities=False)
    x2.sort()
    x = f2.parent().gen()
    ff = [x-x2i for x2i in x2]
    model = "minimal" if E.base_field() is QQ else None
    isogs = [E.isogeny(f, model=model) for f in ff]
    return isogs

def isogenies_3(E):
    """
    Returns a list of all 3-isogenies with domain ``E``.

    INPUT:

    - ``E`` -- an elliptic curve.

    OUTPUT:

    (list) 3-isogenies with domain ``E``.  In general these are
    normalised, but over `\QQ` the codomain is a minimal model.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_3
        sage: E = EllipticCurve(GF(17), [1,1])
        sage: [phi.codomain().ainvs() for phi in isogenies_3(E)]
        [(0, 0, 0, 9, 7), (0, 0, 0, 0, 1)]

        sage: E = EllipticCurve(GF(17^2,'a'), [1,1])
        sage: [phi.codomain().ainvs() for phi in isogenies_3(E)]
        [(0, 0, 0, 9, 7), (0, 0, 0, 0, 1), (0, 0, 0, 5*a + 1, a + 13), (0, 0, 0, 12*a + 6, 16*a + 14)]

        sage: E = EllipticCurve('19a1')
        sage: [phi.codomain().ainvs() for phi in isogenies_3(E)]
        [(0, 1, 1, 1, 0), (0, 1, 1, -769, -8470)]

        sage: E = EllipticCurve([1,1])
        sage: [phi.codomain().ainvs() for phi in isogenies_3(E)]
        []
    """
    f3 = E.division_polynomial(3)
    x3 = f3.roots(multiplicities=False)
    x3.sort()
    x = f3.parent().gen()
    ff = [x-x3i for x3i in x3]
    model = "minimal" if E.base_field() is QQ else None
    isogs = [E.isogeny(f, model=model) for f in ff]
    return isogs

# 6 special cases: `l` = 5, 7, 13 and `j` = 0, 1728.

def isogenies_5_0(E):
    """
    Returns a list of all the 5-isogenies  with domain ``E`` when the
    j-invariant is 0.

    OUTPUT:

    (list) 5-isogenies with codomain E.  In general these are
    normalised, but over `\QQ` the codomain is a minimal model.

    .. note::

       This implementation requires that the characteristic is not 2,
       3 or 5.

    .. note::

       This function would normally be invoked indirectly via ``E.isogenies_prime_degree(5)``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_5_0
        sage: E = EllipticCurve([0,12])
        sage: isogenies_5_0(E)
        []

        sage: E = EllipticCurve(GF(13^2,'a'),[0,-3])
        sage: isogenies_5_0(E)
        [Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + (4*a+6)*x + (2*a+10) over Finite Field in a of size 13^2, Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + (12*a+5)*x + (2*a+10) over Finite Field in a of size 13^2, Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + (10*a+2)*x + (2*a+10) over Finite Field in a of size 13^2, Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + (3*a+12)*x + (11*a+12) over Finite Field in a of size 13^2, Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + (a+4)*x + (11*a+12) over Finite Field in a of size 13^2, Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + (9*a+10)*x + (11*a+12) over Finite Field in a of size 13^2]

        sage: K.<a> = NumberField(x**6-320*x**3-320)
        sage: E = EllipticCurve(K,[0,0,1,0,0])
        sage: isogenies_5_0(E)
        [Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 over Number Field in a with defining polynomial x^6 - 320*x^3 - 320 to Elliptic Curve defined by y^2 = x^3 + (a^5-400*a^2)*x + (280*a^3-3120) over Number Field in a with defining polynomial x^6 - 320*x^3 - 320,
        Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 over Number Field in a with defining polynomial x^6 - 320*x^3 - 320 to Elliptic Curve defined by y^2 = x^3 + (23/2*a^5-3700*a^2)*x + (-280*a^3+86480) over Number Field in a with defining polynomial x^6 - 320*x^3 - 320]

    """
    F = E.base_field()
    if E.j_invariant() != 0:
        raise ValueError, "j-invariant must be 0."
    if F.characteristic() in [2,3,5]:
        raise NotImplementedError, "Not implemented in characteristic 2, 3 or 5."
    if not F(5).is_square():
        return []
    Ew = E.short_weierstrass_model()
    a = Ew.a6()
    x = polygen(F)
    betas = (x**6-160*a*x**3-80*a**2).roots(multiplicities=False)
    betas.sort()
    if len(betas)==0:
        return []
    gammas = [(beta**2 *(beta**3-140*a))/(120*a) for beta in betas]
    model = "minimal" if F is QQ else None
    isogs = [Ew.isogeny(x**2+beta*x+gamma, model=model) for beta,gamma in zip(betas,gammas)]
    iso = E.isomorphism_to(Ew)
    [isog.set_pre_isomorphism(iso) for isog in isogs]
    return isogs

def isogenies_5_1728(E):
    """
    Returns a list of 5-isogenies with domain ``E`` when the j-invariant is
    1728.

    OUTPUT:

    (list) 5-isogenies with codomain E.  In general these are
    normalised; but if `-1` is a square then there are two
    endomorphisms of degree `5`, for which the codomain is the same as
    the domain curve; and over `\QQ`, the codomain is a minimal model.

    .. note::

       This implementation requires that the characteristic is not 2,
       3 or 5.

    .. note::

       This function would normally be invoked indirectly via ``E.isogenies_prime_degree(5)``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_5_1728
        sage: E = EllipticCurve([7,0])
        sage: isogenies_5_1728(E)
        []

        sage: E = EllipticCurve(GF(13),[11,0])
        sage: isogenies_5_1728(E)
        [Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 11*x over Finite Field of size 13 to Elliptic Curve defined by y^2 = x^3 + 11*x over Finite Field of size 13,
        Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + 11*x over Finite Field of size 13 to Elliptic Curve defined by y^2 = x^3 + 11*x over Finite Field of size 13]

    An example of endomorphisms of degree 5::

        sage: K.<i> = QuadraticField(-1)
        sage: E = EllipticCurve(K,[0,0,0,1,0])
        sage: isogenies_5_1728(E)
        [Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1 to Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1,
        Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1 to Elliptic Curve defined by y^2 = x^3 + x over Number Field in i with defining polynomial x^2 + 1]
        sage: _[0].rational_maps()
        (((4/25*i + 3/25)*x^5 + (4/5*i - 2/5)*x^3 - x)/(x^4 + (-4/5*i + 2/5)*x^2 + (-4/25*i - 3/25)),
         ((11/125*i + 2/125)*x^6*y + (-23/125*i + 64/125)*x^4*y + (141/125*i + 162/125)*x^2*y + (3/25*i - 4/25)*y)/(x^6 + (-6/5*i + 3/5)*x^4 + (-12/25*i - 9/25)*x^2 + (2/125*i - 11/125)))

    An example of 5-isogenies over a number field::

        sage: K.<a> = NumberField(x**4+20*x**2-80)
        sage: K(5).is_square() #necessary but not sufficient!
        True
        sage: E = EllipticCurve(K,[0,0,0,1,0])
        sage: isogenies_5_1728(E)
        [Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x over Number Field in a with defining polynomial x^4 + 20*x^2 - 80 to Elliptic Curve defined by y^2 = x^3 + (-20*a^2-39)*x + (35*a^3+112*a) over Number Field in a with defining polynomial x^4 + 20*x^2 - 80,
        Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x over Number Field in a with defining polynomial x^4 + 20*x^2 - 80 to Elliptic Curve defined by y^2 = x^3 + (-20*a^2-39)*x + (-35*a^3-112*a) over Number Field in a with defining polynomial x^4 + 20*x^2 - 80]

    """
    F = E.base_field()
    if E.j_invariant() != 1728:
        raise ValueError, "j-invariant must be 1728."
    if F.characteristic() in [2,3,5]:
        raise NotImplementedError, "Not implemented in characteristic 2, 3 or 5."
    model = "minimal" if F is QQ else None
    # quick test for a negative answer (from Fricke module)
    square5 = F(5).is_square()
    square1 = F(-1).is_square()
    if not square5 and not square1:
        return []
    Ew = E.short_weierstrass_model()
    iso = E.isomorphism_to(Ew)
    a = Ew.a4()
    x = polygen(F)
    isogs = []
    # 2 cases
    # Type 1: if -1 is a square we have 2 endomorphisms
    if square1:
        i = F(-1).sqrt()
        isogs = [Ew.isogeny(f) for f in [x**2+a/(1+2*i), x**2+a/(1-2*i)]]
        [isog.set_post_isomorphism(isog.codomain().isomorphism_to(E)) for isog in isogs]
    # Type 2: if 5 is a square we have up to 4 (non-endomorphism) isogenies
    if square5:
        betas = (x**4+20*a*x**2-80*a**2).roots(multiplicities=False)
        betas.sort()
        gammas = [a*(beta**2-2)/6 for beta in betas]
        isogs += [Ew.isogeny(x**2+beta*x+gamma, model=model) for beta,gamma in zip(betas,gammas)]
    [isog.set_pre_isomorphism(iso) for isog in isogs]
    return isogs

def isogenies_7_0(E):
    """
    Returns list of all 7-isogenies from E when the j-invariant is 0.

    OUTPUT:

    (list) 7-isogenies with codomain E.  In general these are
    normalised; but if `-3` is a square then there are two
    endomorphisms of degree `7`, for which the codomain is the same as
    the domain; and over `\QQ`, the codomain is a minimal model.

    .. note::

       This implementation requires that the characteristic is not 2,
       3 or 7.

    .. note::

       This function would normally be invoked indirectly via ``E.isogenies_prime_degree(7)``.

    EXAMPLES:

    First some examples of endomorphisms::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_7_0
        sage: K.<r> = QuadraticField(-3)
        sage: E = EllipticCurve(K, [0,1])
        sage: isogenies_7_0(E)
        [Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 1 over Number Field in r with defining polynomial x^2 + 3 to Elliptic Curve defined by y^2 = x^3 + 1 over Number Field in r with defining polynomial x^2 + 3,
        Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 1 over Number Field in r with defining polynomial x^2 + 3 to Elliptic Curve defined by y^2 = x^3 + 1 over Number Field in r with defining polynomial x^2 + 3]

        sage: E = EllipticCurve(GF(13^2,'a'),[0,-3])
        sage: isogenies_7_0(E)
        [Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2, Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2 to Elliptic Curve defined by y^2 = x^3 + 10 over Finite Field in a of size 13^2]

    Now some examples of 7-isogenies which are not endomorphisms::

        sage: K = GF(101)
        sage: E = EllipticCurve(K, [0,1])
        sage: isogenies_7_0(E)
        [Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 55*x + 100 over Finite Field of size 101, Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 83*x + 26 over Finite Field of size 101]

    Examples over a number field::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_7_0
        sage: E = EllipticCurve('27a1').change_ring(QuadraticField(-3,'r'))
        sage: isogenies_7_0(E)
        [Isogeny of degree 7 from Elliptic Curve defined by y^2 + y = x^3 + (-7) over Number Field in r with defining polynomial x^2 + 3 to Elliptic Curve defined by y^2 + y = x^3 + (-7) over Number Field in r with defining polynomial x^2 + 3,
        Isogeny of degree 7 from Elliptic Curve defined by y^2 + y = x^3 + (-7) over Number Field in r with defining polynomial x^2 + 3 to Elliptic Curve defined by y^2 + y = x^3 + (-7) over Number Field in r with defining polynomial x^2 + 3]

        sage: K.<a> = NumberField(x^6 + 1512*x^3 - 21168)
        sage: E = EllipticCurve(K, [0,1])
        sage: isogs = isogenies_7_0(E)
        sage: [phi.codomain().a_invariants() for phi in isogs]
        [(0, 0, 0, -5/294*a^5 - 300/7*a^2, -55/2*a^3 - 1133),
        (0, 0, 0, -295/1176*a^5 - 5385/14*a^2, 55/2*a^3 + 40447)]
        sage: [phi.codomain().j_invariant() for phi in isogs]
        [158428486656000/7*a^3 - 313976217600000,
        -158428486656000/7*a^3 - 34534529335296000]
    """
    if E.j_invariant()!=0:
        raise ValueError, "j-invariant must be 0."
    F = E.base_field()
    if F.characteristic() in [2,3,7]:
        raise NotImplementedError, "Not implemented when the characteristic of the base field is 2, 3 or 7."
    x = polygen(F)
    Ew = E.short_weierstrass_model()
    iso = E.isomorphism_to(Ew)
    a = Ew.a6()
    model = "minimal" if F is QQ else None

    # there will be 2 endomorphisms if -3 is a square:

    ts = (x**2+3).roots(multiplicities=False)
    ts.sort()
    kers = [7*x-(2+6*t) for t in ts]
    kers = [k(x**3/a).monic() for k in kers]
    isogs = [Ew.isogeny(k,model=model) for k in kers]
    if len(isogs)>0:
        [endo.set_post_isomorphism(endo.codomain().isomorphism_to(E)) for endo in isogs]

    # we may have up to 6 other isogenies:
    ts = (x**2-21).roots(multiplicities=False)
    for t0 in ts:
        s3 = a/(28+6*t0)
        ss = (x**3-s3).roots(multiplicities=False)
        ss.sort()
        ker = x**3 - 2*t0*x**2 - 4*t0*x + 4*t0 + 28
        kers = [ker(x/s).monic() for s in ss]
        isogs += [Ew.isogeny(k, model=model) for k in kers]

    [isog.set_pre_isomorphism(iso) for isog in isogs]
    return isogs

def isogenies_7_1728(E):
    """
    Returns list of all 7-isogenies from E when the j-invariant is 1728.

    OUTPUT:

    (list) 7-isogenies with codomain E.  In general these are
    normalised; but over `\QQ` the codomain is a minimal model.

    .. note::

       This implementation requires that the characteristic is not 2,
       3, or 7.

    .. note::

       This function would normally be invoked indirectly via ``E.isogenies_prime_degree(7)``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_7_1728
        sage: E = EllipticCurve(GF(47), [1, 0])
        sage: isogenies_7_1728(E)
        [Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 47 to Elliptic Curve defined by y^2 = x^3 + 26 over Finite Field of size 47,
        Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 47 to Elliptic Curve defined by y^2 = x^3 + 21 over Finite Field of size 47]

    An example in characteristic 53 (for which an earlier implementation did not work)::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_7_1728
        sage: E = EllipticCurve(GF(53), [1, 0])
        sage: isogenies_7_1728(E)
        []
        sage: E = EllipticCurve(GF(53^2,'a'), [1, 0])
        sage: [iso.codomain().ainvs() for iso in isogenies_7_1728(E)]
        [(0, 0, 0, 36, 19*a + 15), (0, 0, 0, 36, 34*a + 38), (0, 0, 0, 33, 39*a + 28), (0, 0, 0, 33, 14*a + 25), (0, 0, 0, 19, 45*a + 16), (0, 0, 0, 19, 8*a + 37), (0, 0, 0, 3, 45*a + 16), (0, 0, 0, 3, 8*a + 37)]

    ::

        sage: K.<a> = NumberField(x^8 + 84*x^6 - 1890*x^4 + 644*x^2 - 567)
        sage: E = EllipticCurve(K, [1, 0])
        sage: isogs = isogenies_7_1728(E)
        sage: [phi.codomain().a_invariants() for phi in isogs]
        [(0,
        0,
        0,
        35/636*a^6 + 55/12*a^4 - 79135/636*a^2 + 1127/212,
        155/636*a^7 + 245/12*a^5 - 313355/636*a^3 - 3577/636*a),
        (0,
        0,
        0,
        35/636*a^6 + 55/12*a^4 - 79135/636*a^2 + 1127/212,
        -155/636*a^7 - 245/12*a^5 + 313355/636*a^3 + 3577/636*a)]
        sage: [phi.codomain().j_invariant() for phi in isogs]
        [-526110256146528/53*a^6 + 183649373229024*a^4 - 3333881559996576/53*a^2 + 2910267397643616/53,
        -526110256146528/53*a^6 + 183649373229024*a^4 - 3333881559996576/53*a^2 + 2910267397643616/53]
        sage: E1 = isogs[0].codomain()
        sage: E2 = isogs[1].codomain()
        sage: E1.is_isomorphic(E2)
        False
        sage: E1.is_quadratic_twist(E2)
        -1
    """
    if E.j_invariant()!=1728:
        raise ValueError, "j_invariant must be 1728 (in base field)."
    F = E.base_field()
    if F.characteristic() in [2,3,7]:
        raise NotImplementedError, "Not implemented when the characteristic of the base field is 2, 3 or 7."
    Ew = E.short_weierstrass_model()
    iso = E.isomorphism_to(Ew)
    a = Ew.a4()

    ts = (Fricke_module(7)-1728).numerator().roots(F,multiplicities=False)
    if len(ts)==0:
        return []
    ts.sort()
    isogs = []
    model = "minimal" if F is QQ else None
    x = polygen(F)
    for t0 in ts:
        s2 = a/t0
        ss = (x**2-s2).roots(multiplicities=False)
        ss.sort()
        ker = 9*x**3 + (-3*t0**3 - 36*t0**2 - 123*t0)*x**2 + (-8*t0**3 - 101*t0**2 - 346*t0 + 35)*x - 7*t0**3 - 88*t0**2 - 296*t0 + 28

        kers = [ker(x/s) for s in ss]
        isogs += [Ew.isogeny(k.monic(), model=model) for k in kers]
    [isog.set_pre_isomorphism(iso) for isog in isogs]
    return isogs

def isogenies_13_0(E):
    """
    Returns list of all 13-isogenies from E when the j-invariant is 0.

    OUTPUT:

    (list) 13-isogenies with codomain E.  In general these are
    normalised; but if `-3` is a square then there are two
    endomorphisms of degree `13`, for which the codomain is the same
    as the domain.

    .. note::

       This implementation requires that the characteristic is not 2,
       3 or 13.

    .. note::

       This function would normally be invoked indirectly via ``E.isogenies_prime_degree(13)``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_13_0

    Endomorphisms of degree 13 will exist when -3 is a square::

        sage: K.<r> = QuadraticField(-3)
        sage: E = EllipticCurve(K, [0, r]); E
        Elliptic Curve defined by y^2 = x^3 + r over Number Field in r with defining polynomial x^2 + 3
        sage: isogenies_13_0(E)
        [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + r over Number Field in r with defining polynomial x^2 + 3 to Elliptic Curve defined by y^2 = x^3 + r over Number Field in r with defining polynomial x^2 + 3,
        Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + r over Number Field in r with defining polynomial x^2 + 3 to Elliptic Curve defined by y^2 = x^3 + r over Number Field in r with defining polynomial x^2 + 3]
        sage: isogenies_13_0(E)[0].rational_maps()
        (((7/338*r + 23/338)*x^13 + (-164/13*r - 420/13)*x^10 + (720/13*r + 3168/13)*x^7 + (3840/13*r - 576/13)*x^4 + (4608/13*r + 2304/13)*x)/(x^12 + (4*r + 36)*x^9 + (1080/13*r + 3816/13)*x^6 + (2112/13*r - 5184/13)*x^3 + (-17280/169*r - 1152/169)), ((18/2197*r + 35/2197)*x^18*y + (23142/2197*r + 35478/2197)*x^15*y + (-1127520/2197*r - 1559664/2197)*x^12*y + (-87744/2197*r + 5992704/2197)*x^9*y + (-6625152/2197*r - 9085824/2197)*x^6*y + (-28919808/2197*r - 2239488/2197)*x^3*y + (-1990656/2197*r - 3870720/2197)*y)/(x^18 + (6*r + 54)*x^15 + (3024/13*r + 11808/13)*x^12 + (31296/13*r + 51840/13)*x^9 + (487296/169*r - 2070144/169)*x^6 + (-940032/169*r + 248832/169)*x^3 + (1990656/2197*r + 3870720/2197)))

    An example of endomorphisms over a finite field::

        sage: K = GF(19^2,'a')
        sage: E = EllipticCurve(j=K(0)); E
        Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in a of size 19^2
        sage: isogenies_13_0(E)
        [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in a of size 19^2 to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in a of size 19^2,
        Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in a of size 19^2 to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in a of size 19^2]
        sage: isogenies_13_0(E)[0].rational_maps()
        ((6*x^13 - 6*x^10 - 3*x^7 + 6*x^4 + x)/(x^12 - 5*x^9 - 9*x^6 - 7*x^3 + 5), (-8*x^18*y - 9*x^15*y + 9*x^12*y - 5*x^9*y + 5*x^6*y - 7*x^3*y + 7*y)/(x^18 + 2*x^15 + 3*x^12 - x^9 + 8*x^6 - 9*x^3 + 7))

    A previous implementation did not work in some characteristics::

        sage: K = GF(29)
        sage: E = EllipticCurve(j=K(0))
        sage: isogenies_13_0(E)
        [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 26*x + 12 over Finite Field of size 29, Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 29 to Elliptic Curve defined by y^2 = x^3 + 16*x + 28 over Finite Field of size 29]

    ::

        sage: K = GF(101)
        sage: E = EllipticCurve(j=K(0)); E.ainvs()
        (0, 0, 0, 0, 1)
        sage: [phi.codomain().ainvs() for phi in isogenies_13_0(E)]
        [(0, 0, 0, 64, 36), (0, 0, 0, 42, 66)]

    ::

        sage: x = polygen(QQ)
        sage: f = x^12 + 78624*x^9 - 130308048*x^6 + 2270840832*x^3 - 54500179968
        sage: K.<a> = NumberField(f)
        sage: E = EllipticCurve(j=K(0)); E.ainvs()
        (0, 0, 0, 0, 1)
        sage: [phi.codomain().ainvs() for phi in isogenies_13_0(E)]
        [(0, 0, 0, -739946459/23857162861049856*a^11 - 2591641747/1062017577504*a^8 + 16583647773233/4248070310016*a^5 - 14310911337/378211388*a^2, 26146225/4248070310016*a^9 + 7327668845/14750244132*a^6 + 174618431365/756422776*a^3 - 378332499709/94552847), (0, 0, 0, 3501275/5964290715262464*a^11 + 24721025/531008788752*a^8 - 47974903745/1062017577504*a^5 - 6773483100/94552847*a^2, 6699581/4248070310016*a^9 + 1826193509/14750244132*a^6 - 182763866047/756422776*a^3 - 321460597/94552847)]
    """
    if E.j_invariant()!=0:
        raise ValueError, "j-invariant must be 0."
    F = E.base_field()
    if F.characteristic() in [2,3,13]:
        raise NotImplementedError, "Not implemented when the characteristic of the base field is 2, 3 or 13."
    Ew = E.short_weierstrass_model()
    iso = E.isomorphism_to(Ew)
    a = Ew.a6()
    model = "minimal" if F is QQ else None
    x = polygen(F)

    # there will be 2 endomorphisms if -3 is a square:
    ts = (x**2+3).roots(multiplicities=False)
    ts.sort()
    kers = [13*x**2 + (78*t + 26)*x + 24*t + 40 for t in ts]
    kers = [k(x**3/a).monic() for k in kers]
    isogs = [Ew.isogeny(k,model=model) for k in kers]
    if len(isogs)>0:
        [endo.set_post_isomorphism(endo.codomain().isomorphism_to(E)) for endo in isogs]

    # we may have up to 12 other isogenies:
    ts = (x**4 + 7*x**3 + 20*x**2 + 19*x + 1).roots(multiplicities=False)
    ts.sort()
    for t0 in ts:
        s3 = a / (6*t0**3 + 32*t0**2 + 68*t0 + 4)
        ss = (x**3-s3).roots(multiplicities=False)
        ss.sort()
        ker = (x**6 + (20*t0**3 + 106*t0**2 + 218*t0 + 4)*x**5
            + (-826*t0**3 - 4424*t0**2 - 9244*t0 - 494)*x**4
            + (13514*t0**3 + 72416*t0**2 + 151416*t0 + 8238)*x**3
            + (-101948*t0**3 - 546304*t0**2 - 1142288*t0 - 62116)*x**2
            + (354472*t0**3 + 1899488*t0**2 + 3971680*t0 + 215960)*x
            - 459424*t0**3 - 2461888*t0**2 - 5147648*t0 - 279904)
        kers = [ker(x/s).monic() for s in ss]
        isogs += [Ew.isogeny(k, model=model) for k in kers]

    [isog.set_pre_isomorphism(iso) for isog in isogs]

    return isogs

def isogenies_13_1728(E):
    """
    Returns list of all 13-isogenies from E when the j-invariant is 1728.

    OUTPUT:

    (list) 13-isogenies with codomain E.  In general these are
    normalised; but if `-1` is a square then there are two
    endomorphisms of degree `13`, for which the codomain is the same
    as the domain; and over `\QQ`, the codomain is a minimal model.

    .. note::

       This implementation requires that the characteristic is not
       2, 3 or 13.

    .. note::

       This function would normally be invoked indirectly via ``E.isogenies_prime_degree(13)``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import isogenies_13_1728

        sage: K.<i> = QuadraticField(-1)
        sage: E = EllipticCurve([0,0,0,i,0]); E.ainvs()
        (0, 0, 0, i, 0)
        sage: isogenies_13_1728(E)
        [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + i*x over Number Field in i with defining polynomial x^2 + 1 to Elliptic Curve defined by y^2 = x^3 + i*x over Number Field in i with defining polynomial x^2 + 1,
        Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + i*x over Number Field in i with defining polynomial x^2 + 1 to Elliptic Curve defined by y^2 = x^3 + i*x over Number Field in i with defining polynomial x^2 + 1]

    ::

        sage: K = GF(83)
        sage: E = EllipticCurve(K, [0,0,0,5,0]); E.ainvs()
        (0, 0, 0, 5, 0)
        sage: isogenies_13_1728(E)
        []
        sage: K = GF(89)
        sage: E = EllipticCurve(K, [0,0,0,5,0]); E.ainvs()
        (0, 0, 0, 5, 0)
        sage: isogenies_13_1728(E)
        [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 89 to Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 89,
        Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 89 to Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 89]

    ::

        sage: K = GF(23)
        sage: E = EllipticCurve(K, [1,0])
        sage: isogenies_13_1728(E)
        [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 23 to Elliptic Curve defined by y^2 = x^3 + 16 over Finite Field of size 23, Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 23 to Elliptic Curve defined by y^2 = x^3 + 7 over Finite Field of size 23]

    ::

        sage: x = polygen(QQ)
        sage: f = x^12 + 1092*x^10 - 432432*x^8 + 6641024*x^6 - 282896640*x^4 - 149879808*x^2 - 349360128
        sage: K.<a> = NumberField(f)
        sage: E = EllipticCurve(K, [1,0])
        sage: [phi.codomain().ainvs() for phi in isogenies_13_1728(E)]
        [(0,
        0,
        0,
        11090413835/20943727039698624*a^10 + 32280103535965/55849938772529664*a^8 - 355655987835845/1551387188125824*a^6 + 19216954517530195/5235931759924656*a^4 - 1079766118721735/5936430566808*a^2 + 156413528482727/8080141604822,
        214217013065/82065216155553792*a^11 + 1217882637605/427423000810176*a^9 - 214645003230565/189965778137856*a^7 + 22973355421236025/1282269002430528*a^5 - 2059145797340695/2544184528632*a^3 - 23198483147321/989405094468*a),
        (0,
        0,
        0,
        11090413835/20943727039698624*a^10 + 32280103535965/55849938772529664*a^8 - 355655987835845/1551387188125824*a^6 + 19216954517530195/5235931759924656*a^4 - 1079766118721735/5936430566808*a^2 + 156413528482727/8080141604822,
        -214217013065/82065216155553792*a^11 - 1217882637605/427423000810176*a^9 + 214645003230565/189965778137856*a^7 - 22973355421236025/1282269002430528*a^5 + 2059145797340695/2544184528632*a^3 + 23198483147321/989405094468*a)]

    """
    if E.j_invariant()!=1728:
        raise ValueError, "j-invariant must be 1728."
    F = E.base_field()
    if F.characteristic() in [2, 3, 7, 13, 167, 233, 271, 1117]:
        raise NotImplementedError, "Not implemented when the characteristic of the base field is 2, 3 or 13."
    Ew = E.short_weierstrass_model()
    iso = E.isomorphism_to(Ew)
    a = Ew.a4()
    model = "minimal" if F is QQ else None
    x = polygen(F)

    # we will have two endomorphisms if -1 is a square:
    ts = (x**2+1).roots(multiplicities=False)
    ts.sort()
    kers = [13*x**3 + (-26*i - 13)*x**2 + (-52*i - 13)*x - 2*i - 3 for i in ts]
    kers = [k(x**2/a).monic() for k in kers]
    isogs = [Ew.isogeny(k,model=model) for k in kers]
    if len(isogs)>0:
        [endo.set_post_isomorphism(endo.codomain().isomorphism_to(E)) for endo in isogs]

    # we may have up to 12 other isogenies:

    ts = (x**6 + 10*x**5 + 46*x**4 + 108*x**3 + 122*x**2 + 38*x - 1).roots(multiplicities=False)
    ts.sort()
    for t0 in ts:
        s2 = a/(66*t0**5 + 630*t0**4 + 2750*t0**3 + 5882*t0**2 + 5414*t0 + 162)
        ss = (x**2-s2).roots(multiplicities=False)
        ss.sort()
        ker = (x**6 + (-66*t0**5 - 630*t0**4 - 2750*t0**3 - 5882*t0**2
              - 5414*t0 - 162)*x**5 + (-21722*t0**5 - 205718*t0**4 -
              890146*t0**3 - 1873338*t0**2 - 1652478*t0 + 61610)*x**4
              + (-3391376*t0**5 - 32162416*t0**4 - 139397232*t0**3 -
              294310576*t0**2 - 261885968*t0 + 6105552)*x**3 +
              (-241695080*t0**5 - 2291695976*t0**4 - 9930313256*t0**3
              - 20956609720*t0**2 - 18625380856*t0 + 469971320)*x**2 +
              (-8085170432*t0**5 - 76663232384*t0**4 -
              332202985024*t0**3 - 701103233152*t0**2 -
              623190845440*t0 + 15598973056)*x - 101980510208*t0**5 -
              966973468160*t0**4 - 4190156868352*t0**3 -
              8843158270336*t0**2 - 7860368751232*t0 + 196854655936)

        kers = [ker(x/s).monic() for s in ss]
        isogs += [Ew.isogeny(k, model=model) for k in kers]

    [isog.set_pre_isomorphism(iso) for isog in isogs]

    return isogs

# Utility function for manipulating isogeny degree matrices

def fill_isogeny_matrix(M):
    """
    Returns a filled isogeny matrix giving all degrees from one giving only prime degrees.

    INPUT:

    - ``M`` -- a square symmetric matrix whose off-diagonal `i`, `j`
      entry is either a prime `l` (if the `i`'th and `j`'th curves
      have an `l`-isogeny between them), otherwise is 0.

    OUTPUT:

    (matrix) a square matrix with entries `1` on the diagonal, and in
    general the `i`, `j` entry is `d>0` if `d` is the minimal degree
    of an isogeny from the `i`'th to the `j`'th curve,

    EXAMPLES::

        sage: M = Matrix([[0, 2, 3, 3, 0, 0], [2, 0, 0, 0, 3, 3], [3, 0, 0, 0, 2, 0], [3, 0, 0, 0, 0, 2], [0, 3, 2, 0, 0, 0], [0, 3, 0, 2, 0, 0]]); M
        [0 2 3 3 0 0]
        [2 0 0 0 3 3]
        [3 0 0 0 2 0]
        [3 0 0 0 0 2]
        [0 3 2 0 0 0]
        [0 3 0 2 0 0]
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix
        sage: fill_isogeny_matrix(M)
        [ 1  2  3  3  6  6]
        [ 2  1  6  6  3  3]
        [ 3  6  1  9  2 18]
        [ 3  6  9  1 18  2]
        [ 6  3  2 18  1  9]
        [ 6  3 18  2  9  1]
    """
    from sage.matrix.all import Matrix
    from sage.rings.infinity import Infinity

    n = M.nrows()
    M0 = copy(M)
    for i in range(n):
        M0[i,i]=1

    def fix(d):
        if d==0: return Infinity
        return d

    def fix2(d):
        if d==Infinity: return 0
        return d

    def pr(M1,M2):
        return Matrix([[fix2(min([fix(M1[i,k]*M2[k,j]) for k in range(n)])) for i in range(n)] for j in range(n)])

    M1 = M0
    M2 = pr(M0,M1)
    while M1!=M2:
        M1 = M2
        M2 = pr(M0,M1)

    return M1

def unfill_isogeny_matrix(M):
    """
    Reverses the action of ``fill_isogeny_matrix``.

    INPUT:

    - ``M`` -- a square symmetric matrix of integers.

    OUTPUT:

    (matrix) a square symmetric matrix obtained from ``M`` by
    replacing non-prime entries with `0`.

    EXAMPLES::

        sage: M = Matrix([[0, 2, 3, 3, 0, 0], [2, 0, 0, 0, 3, 3], [3, 0, 0, 0, 2, 0], [3, 0, 0, 0, 0, 2], [0, 3, 2, 0, 0, 0], [0, 3, 0, 2, 0, 0]]); M
        [0 2 3 3 0 0]
        [2 0 0 0 3 3]
        [3 0 0 0 2 0]
        [3 0 0 0 0 2]
        [0 3 2 0 0 0]
        [0 3 0 2 0 0]
        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix, unfill_isogeny_matrix
        sage: M1 = fill_isogeny_matrix(M); M1
        [ 1  2  3  3  6  6]
        [ 2  1  6  6  3  3]
        [ 3  6  1  9  2 18]
        [ 3  6  9  1 18  2]
        [ 6  3  2 18  1  9]
        [ 6  3 18  2  9  1]
        sage: unfill_isogeny_matrix(M1)
        [0 2 3 3 0 0]
        [2 0 0 0 3 3]
        [3 0 0 0 2 0]
        [3 0 0 0 0 2]
        [0 3 2 0 0 0]
        [0 3 0 2 0 0]
        sage: unfill_isogeny_matrix(M1) == M
        True
    """
    from sage.matrix.all import Matrix
    from sage.rings.infinity import Infinity

    n = M.nrows()
    M1 = copy(M)
    zero = Integer(0)
    for i in range(n):
        M1[i,i] = zero
        for j in range(i):
            if not M1[i,j].is_prime():
                M1[i,j] = zero
                M1[j,i] = zero
    return M1
