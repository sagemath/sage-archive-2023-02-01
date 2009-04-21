"""
Elliptic curve constructor

AUTHORS:

- William Stein (2005): Initial version

- John Cremona (2008-01): EllipticCurve(j) fixed for all cases
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import sage.rings.all as rings

from sage.structure.sequence import Sequence
from sage.structure.element import parent
from sage.calculus.calculus import SR, SymbolicEquation


def EllipticCurve(x=None, y=None, j=None):
    r"""
    There are several ways to construct an elliptic curve:

    .. math::

       y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6.


    - EllipticCurve([a1,a2,a3,a4,a6]): Elliptic curve with given
      a-invariants. The invariants are coerced into the parent of the
      first element. If all are integers, they are coerced into the
      rational numbers.

    - EllipticCurve([a4,a6]): Same as above, but a1=a2=a3=0.

    - EllipticCurve(label): Returns the elliptic curve over Q from the
      Cremona database with the given label. The label is a string, such
      as "11a" or "37b2". The letters in the label *must* be lower case
      (Cremona's new labeling).

    - EllipticCurve(R, [a1,a2,a3,a4,a6]): Create the elliptic curve
      over R with given a-invariants. Here R can be an arbitrary ring.
      Note that addition need not be defined.


    - EllipticCurve(j): Return an elliptic curve with j-invariant
      `j`.  Warning: this is deprecated.  Use ``EllipticCurve_from_j(j)``
      or ``EllipticCurve(j=j)`` instead.


    EXAMPLES: We illustrate creating elliptic curves.

    ::

        sage: EllipticCurve([0,0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    We create a curve from a Cremona label::

        sage: EllipticCurve('37b2')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field
        sage: EllipticCurve('5077a')
        Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
        sage: EllipticCurve('389a')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    We create curves over a finite field as follows::

        sage: EllipticCurve([GF(5)(0),0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5
        sage: EllipticCurve(GF(5), [0, 0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    The following is a curve over the complex numbers::

        sage: E = EllipticCurve(CC, [0,0,1,-1,0])
        sage: E
        Elliptic Curve defined by y^2 + 1.00000000000000*y = x^3 + (-1.00000000000000)*x over Complex Field with 53 bits of precision
        sage: E.j_invariant()
        2988.97297297297

    We can also create elliptic curves by giving the Weierstrass equation::

        sage: x, y = var('x,y')
        sage: EllipticCurve(y^2 + y ==  x^3 + x - 9)
        Elliptic Curve defined by y^2 + y = x^3 + x - 9 over Rational Field

        sage: R.<x,y> = GF(5)[]
        sage: EllipticCurve(x^3 + x^2 + 2 - y^2 - y*x)
        Elliptic Curve defined by y^2 + x*y  = x^3 + x^2 + 2 over Finite Field of size 5

    We can explicitly specify the `j`-invariant::

        sage: E = EllipticCurve(j=1728); E; E.j_invariant(); E.label()
        Elliptic Curve defined by y^2 = x^3 - x over Rational Field
        1728
        '32a2'

        sage: E = EllipticCurve(j=GF(5)(2)); E; E.j_invariant()
        Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 5
        2


    TESTS::

        sage: R = ZZ['u', 'v']
        sage: EllipticCurve(R, [1,1])
        Elliptic Curve defined by y^2 = x^3 + x + 1 over Multivariate Polynomial Ring in u, v
        over Integer Ring

    We create a curve and a point over QQbar::

        sage: E = EllipticCurve(QQbar,[0,1])
        sage: E(0)
        (0 : 1 : 0)
    """
    import ell_generic, ell_finite_field, ell_number_field, ell_rational_field, ell_padic_field  # here to avoid circular includes

    if j is not None:
        return EllipticCurve_from_j(j)
    assert x is not None

    if isinstance(x, SymbolicEquation):
        x = x.left() - x.right()

    if parent(x) is SR:
        x = x._polynomial_(rings.QQ['x', 'y'])

    if rings.is_MPolynomial(x) and y is None:
        f = x
        if f.degree() != 3:
            raise ValueError, "Elliptic curves must be defined by a cubic polynomial."
        if f.degrees() == (3,2):
            x, y = f.parent().gens()
        elif f.degree() == (2,3):
            y, x = f.parent().gens()
        elif len(f.parent().gens()) == 2 or len(f.parent().gens()) == 3 and f.is_homogeneous():
            # We'd need a point too...
            raise NotImplementedError, "Construction of an elliptic curve from a generic cubic not yet implemented."
        else:
            raise ValueError, "Defining polynomial must be a cubic polynomial in two variables."

        try:
            if f.coefficient(x**3) < 0:
                f = -f
            # is there a nicer way to extract the coefficients?
            a1 = a2 = a3 = a4 = a6 = 0
            for coeff, mon in f:
                if mon == x**3:
                    assert coeff == 1
                elif mon == x**2:
                    a2 = coeff
                elif mon == x:
                    a4 = coeff
                elif mon == 1:
                    a6 = coeff
                elif mon == y**2:
                    assert coeff == -1
                elif mon == x*y:
                    a1 = -coeff
                elif mon == y:
                    a3 = -coeff
                else:
                    assert False
            return EllipticCurve([a1, a2, a3, a4, a6])
        except AssertionError:
            raise NotImplementedError, "Construction of an elliptic curve from a generic cubic not yet implemented."

    if rings.is_Ring(x):
        if rings.is_RationalField(x):
            return ell_rational_field.EllipticCurve_rational_field(x, y)
        elif rings.is_FiniteField(x):
            return ell_finite_field.EllipticCurve_finite_field(x, y)
        elif rings.is_pAdicField(x):
            return ell_padic_field.EllipticCurve_padic_field(x, y)
        elif rings.is_NumberField(x):
            return ell_number_field.EllipticCurve_number_field(x, y)
        else:
            return ell_generic.EllipticCurve_generic(x, y)

    if isinstance(x, str):
        return ell_rational_field.EllipticCurve_rational_field(x)

    if rings.is_RingElement(x) and y is None:
        from sage.misc.misc import deprecation
        deprecation("'EllipticCurve(j)' is deprecated; use 'EllipticCurve_from_j(j)' or 'EllipticCurve(j=j)' instead.")
        # Fixed for all characteristics and cases by John Cremona
        j=x
        F=j.parent().fraction_field()
        char=F.characteristic()
        if char==2:
            if j==0:
                return EllipticCurve(F, [ 0, 0, 1, 0, 0 ])
            else:
                return EllipticCurve(F, [ 1, 0, 0, 0, 1/j ])
        if char==3:
            if j==0:
                return EllipticCurve(F, [ 0, 0, 0, 1, 0 ])
            else:
                return EllipticCurve(F, [ 0, j, 0, 0, -j**2 ])
        if j == 0:
            return EllipticCurve(F, [ 0, 0, 0, 0, 1 ])
        if j == 1728:
            return EllipticCurve(F, [ 0, 0, 0, 1, 0 ])
        k=j-1728
        return EllipticCurve(F, [0,0,0,-3*j*k, -2*j*k**2])

    if not isinstance(x,list):
        raise TypeError, "invalid input to EllipticCurve constructor"

    x = Sequence(x)
    if not (len(x) in [2,5]):
        raise ValueError, "sequence of coefficients must have length at 2 or 5"
    R = x.universe()

    if isinstance(x[0], (rings.Rational, rings.Integer, int, long)):
        return ell_rational_field.EllipticCurve_rational_field(x, y)

    elif rings.is_NumberField(R):
        return ell_number_field.EllipticCurve_number_field(x, y)

    elif rings.is_pAdicField(R):
        return ell_padic_field.EllipticCurve_padic_field(x, y)

    elif isinstance(x[0], rings.FiniteFieldElement) or rings.is_IntegerMod(x[0]):
        return ell_finite_field.EllipticCurve_finite_field(x, y)

    else:

        return ell_generic.EllipticCurve_generic(x, y)

def EllipticCurve_from_c4c6(c4, c6):
    """
    Return an elliptic curve with given `c_4` and
    `c_6` invariants.

    EXAMPLES::

        sage: E = EllipticCurve_from_c4c6(17, -2005)
        sage: E
        Elliptic Curve defined by y^2  = x^3 - 17/48*x + 2005/864 over Rational Field
        sage: E.c_invariants()
        (17, -2005)
    """
    try:
        K = c4.parent()
    except AttributeError:
        K = rings.RationalField()
    if not rings.is_Field(K):
        K = K.fraction_field()
    return EllipticCurve([-K(c4)/K(48), -K(c6)/K(864)])

def EllipticCurve_from_j(j):
    """
    Return an elliptic curve with given `j`-invariant.

    EXAMPLES::

        sage: E = EllipticCurve_from_j(0); E; E.j_invariant(); E.label()
        Elliptic Curve defined by y^2 + y = x^3 over Rational Field
        0
        '27a3'

        sage: E = EllipticCurve_from_j(1728); E; E.j_invariant(); E.label()
        Elliptic Curve defined by y^2 = x^3 - x over Rational Field
        1728
        '32a2'

        sage: E = EllipticCurve_from_j(1); E; E.j_invariant()
        Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455 over Rational Field
        1

    """
    try:
        K = j.parent()
    except AttributeError:
        K = rings.RationalField()
    if not rings.is_Field(K):
        K = K.fraction_field()

    char=K.characteristic()
    if char==2:
        if j == 0:
            return EllipticCurve(K, [ 0, 0, 1, 0, 0 ])
        else:
            return EllipticCurve(K, [ 1, 0, 0, 0, 1/j ])
    if char == 3:
        if j==0:
            return EllipticCurve(K, [ 0, 0, 0, 1, 0 ])
        else:
            return EllipticCurve(K, [ 0, j, 0, 0, -j**2 ])

    if K is rings.RationalField():
        # we construct the minimal twist, i.e. the curve with minimal
        # conductor with this j_invariant:

        if j == 0:
            return EllipticCurve(K, [ 0, 0, 1, 0, 0 ]) # 27a3
        if j == 1728:
            return EllipticCurve(K, [ 0, 0, 0, -1, 0 ]) # 32a2

        n = j.numerator()
        m = n-1728*j.denominator()
        a4 = -3*n*m
        a6 = -2*n*m**2

        # Now E=[0,0,0,a4,a6] has j-invariant j=n/d

        from sage.sets.set import Set
        for p in Set(n.prime_divisors()+m.prime_divisors()):
            e = min(a4.valuation(p)//2,a6.valuation(p)//3)
            if e>0:
                p  = p**e
                a4 /= p**2
                a6 /= p**3

        # Now E=[0,0,0,a4,a6] is minimal at all p != 2,3

        tw = [-1,2,-2,3,-3,6,-6]
        E1 = EllipticCurve([0,0,0,a4,a6])
        Elist = [E1] + [E1.quadratic_twist(t) for t in tw]
        crv_cmp = lambda E,F: cmp(E.conductor(),F.conductor())
        Elist.sort(cmp=crv_cmp)
        return Elist[0]

    # defaults for all other fields:
    if j == 0:
        return EllipticCurve(K, [ 0, 0, 0, 0, 1 ])
    if j == 1728:
        return EllipticCurve(K, [ 0, 0, 0, 1, 0 ])
    k=j-1728
    return EllipticCurve(K, [0,0,0,-3*j*k, -2*j*k**2])

def EllipticCurve_from_cubic(F, P):
    r"""
    Given a nonsingular homogenous cubic polynomial F over
    `\QQ` in three variables x, y, z and a projective solution
    P=[a,b,c] to F(P)=0, find the minimal Weierstrass equation of the
    elliptic curve over `\QQ` that is isomorphic to the curve
    defined by `F=0`.

    .. note::

       USES MAGMA - This function will not work on computers that
       do not have magma installed. (HELP WANTED - somebody implement this
       independent of MAGMA.)

    EXAMPLES: First we find that the Fermat cubic is isomorphic to the
    curve with Cremona label 27a1::

        sage: E = EllipticCurve_from_cubic('x^3 + y^3 + z^3', [1,-1,0])  # optional - magma
        sage: E         # optional - magma
        Elliptic Curve defined by y^2 + y = x^3 - 7 over Rational Field
        sage: E.cremona_label()     # optional - magma
        '27a1'

    Next we find the minimal model and conductor of the Jacobian of the
    Selmer curve.

    ::

        sage: E = EllipticCurve_from_cubic('x^3 + y^3 + 60*z^3', [1,-1,0])   # optional - magma
        sage: E                # optional - magma
        Elliptic Curve defined by y^2  = x^3 - 24300 over Rational Field
        sage: E.conductor()    # optional - magma
        24300
    """
    from sage.interfaces.all import magma
    magma.eval("P<x,y,z> := ProjectivePlane(RationalField());")
    cmd = 'aInvariants(MinimalModel(EllipticCurve(Curve(Scheme(P, %s)),P!%s)));'%(F, P)
    s = magma.eval(cmd)
    return EllipticCurve(rings.RationalField(), eval(s))

def EllipticCurves_with_good_reduction_outside_S(S=[], proof=None, verbose=False):
    r"""
    Returns a sorted list of all elliptic curves defined over `Q`
    with good reduction outside the set `S` of primes.

    INPUT:

        -  ``S`` - list of primes (default: empty list).

        - ``proof`` - True/False (default True): the MW basis for
          auxiliary curves will be computed with this proof flag.

        - ``verbose`` - True/False (default False): if True, some
          details of the computation will be output.

    .. note::

        Proof flag: The algorithm used requires determining all
        S-integral points on several auxiliary curves, which in turn
        requires the computation of their generators.  This is not
        always possible (even in theory) using current knowledge.

        The value of this flag is passed to the function which
        computes generators of various auxiliary elliptic curves, in
        order to find their S-integral points.  Set to False if the
        default (True) causes warning messages, but note that you can
        then not rely on the set of curves returned being
        complete.

    EXAMPLES::

        sage: EllipticCurves_with_good_reduction_outside_S([])
        []
        sage: elist = EllipticCurves_with_good_reduction_outside_S([2])
        sage: elist
        [Elliptic Curve defined by y^2 = x^3 + 4*x over Rational Field,
        Elliptic Curve defined by y^2 = x^3 - x over Rational Field,
        ...
        Elliptic Curve defined by y^2 = x^3 - x^2 - 13*x + 21 over Rational Field]
        sage: len(elist)
        24
        sage: ', '.join([e.label() for e in elist])
        '32a1, 32a2, 32a3, 32a4, 64a1, 64a2, 64a3, 64a4, 128a1, 128a2, 128b1, 128b2, 128c1, 128c2, 128d1, 128d2, 256a1, 256a2, 256b1, 256b2, 256c1, 256c2, 256d1, 256d2'

    # Without the "Proof=False", this example gives two warnings:
        sage: elist = EllipticCurves_with_good_reduction_outside_S([11],proof=False)
        sage: len(elist)
        12
        sage: ', '.join([e.label() for e in elist])
        '11a1, 11a2, 11a3, 121a1, 121a2, 121b1, 121b2, 121c1, 121c2, 121d1, 121d2, 121d3'

        sage: elist = EllipticCurves_with_good_reduction_outside_S([2,3]) # long time (~35s)
        sage: len(elist) # long time
        752
        sage: max([e.conductor() for e in elist]) # long time
        62208
        sage: [N.factor() for N in Set([e.conductor() for e in elist])] # long time
        [2^7,
        2^8,
        2^3 * 3^4,
        2^2 * 3^3,
        2^8 * 3^4,
        2^4 * 3^4,
        2^3 * 3,
        2^7 * 3,
        2^3 * 3^5,
        3^3,
        2^8 * 3,
        2^5 * 3^4,
        2^4 * 3,
        2 * 3^4,
        2^2 * 3^2,
        2^6 * 3^4,
        2^6,
        2^7 * 3^2,
        2^4 * 3^5,
        2^4 * 3^3,
        2 * 3^3,
        2^6 * 3^3,
        2^6 * 3,
        2^5,
        2^2 * 3^4,
        2^3 * 3^2,
        2^5 * 3,
        2^7 * 3^4,
        2^2 * 3^5,
        2^8 * 3^2,
        2^5 * 3^2,
        2^7 * 3^5,
        2^8 * 3^5,
        2^3 * 3^3,
        2^8 * 3^3,
        2^5 * 3^5,
        2^4 * 3^2,
        2 * 3^5,
        2^5 * 3^3,
        2^6 * 3^5,
        2^7 * 3^3,
        3^5,
        2^6 * 3^2]
    """
    from ell_egros import (egros_from_jlist, egros_get_j)
    return egros_from_jlist(egros_get_j(S, proof=proof, verbose=verbose), S)
