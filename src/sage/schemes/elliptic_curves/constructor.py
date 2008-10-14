"""
Elliptic curve constructor

AUTHORS:
   * William Stein (2005) -- Initial version
   * John Cremona (Jan 2008) -- EllipticCurve(j) fixed for all cases
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


def EllipticCurve(x, y=None):
    r"""
    There are several ways to construct an elliptic curve:
      $$
          y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6.
      $$

        -- EllipticCurve([a1,a2,a3,a4,a6]): Elliptic curve with given
           a-invariants.  The invariants are coerced into the parent
           of the first element.  If all are integers, they are coerced
           into the rational numbers.

        -- EllipticCurve([a4,a6]): Same as above, but a1=a2=a3=0.

        -- EllipticCurve(label): Returns the elliptic curve over Q
           from the Cremona database with the given label.  The label
           is a string, such as "11a" or "37b2".  The letters in the
           label \emph{must} be lower case (Cremona's new labeling).

        -- EllipticCurve(R, [a1,a2,a3,a4,a6]): Create the elliptic
           curve over R with given a-invariants.  Here R can be an
           arbitrary ring.  Note that addition need not be defined.

        -- EllipticCurve(j): Return an elliptic curve with j-invariant
           $j$.

    EXAMPLES:
    We illustrate creating elliptic curves.

        sage: EllipticCurve([0,0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field


    We create a curve from a Cremona label:
        sage: EllipticCurve('37b2')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field
        sage: EllipticCurve('5077a')
        Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
        sage: EllipticCurve('389a')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    We create curves over a finite field as follows:
        sage: EllipticCurve([GF(5)(0),0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5
        sage: EllipticCurve(GF(5), [0, 0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    The following is a curve over the complex numbers:
        sage: E = EllipticCurve(CC, [0,0,1,-1,0])
        sage: E
        Elliptic Curve defined by y^2 + 1.00000000000000*y = x^3 + (-1.00000000000000)*x over Complex Field with 53 bits of precision
        sage: E.j_invariant()
        2988.97297297297

    TESTS:
        sage: R = ZZ['u', 'v']
        sage: EllipticCurve(R, [1,1])
        Elliptic Curve defined by y^2  = x^3 + x +1 over Multivariate Polynomial Ring in u, v
        over Integer Ring

    We create a curve and a point over QQbar:
        sage: E = EllipticCurve(QQbar,[0,1])
        sage: E(0)
        (0 : 1 : 0)
    """
    # TODO - - implement
        #sage: E = EllipticCurve(ZZ, [0, 0,1,-1,0])
        #sage: E
        #Elliptic Curve defined by y^2 + y = x^3 - x over Integer Ring

    #Of course, arithmetic on elliptic curves over Z need not be defined:
        #sage: P = E([0,0])
        #sage: P + P + P + P
        #(2, -3)
        #sage: P + P + P + P + P
        #Traceback (most recent call last):
        #...
        #ArithmeticError: Point (1/4, -5/8) is not on curve.
    #
    import ell_generic, ell_finite_field, ell_number_field, ell_rational_field, ell_padic_field  # here to avoid circular includes

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
    Return an elliptic curve with given $c_4$ and $c_6$ invariants.

    EXAMPLES:
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

def EllipticCurve_from_cubic(F, P):
    r"""
    Given a nonsingular homogenous cubic polynomial F over $\Q$ in
    three variables x, y, z and a projective solution P=[a,b,c] to
    F(P)=0, find the minimal Weierstrass equation of the elliptic
    curve over $\Q$ that is isomorphic to the curve defined by $F=0$.

    \note{USES MAGMA -- This function will not work on computers that
    do not have magma installed.  (HELP WANTED -- somebody implement
    this independent of MAGMA.)}

    EXAMPLES:
    First we find that the Fermat cubic is isomorphic to the
    curve with Cremona label 27a1:

        sage: E = EllipticCurve_from_cubic('x^3 + y^3 + z^3', [1,-1,0])  # optional -- requires magma
        sage: E         # optional
        Elliptic Curve defined by y^2 + y = x^3 - 7 over Rational Field
        sage: E.cremona_label()     # optional
        '27a1'

    Next we find the minimal model and conductor of the Jacobian
    of the Selmer curve.
        sage: E = EllipticCurve_from_cubic('x^3 + y^3 + 60*z^3', [1,-1,0])   # optional
        sage: E            # optional
        Elliptic Curve defined by y^2  = x^3 - 24300 over Rational Field
        sage: E.conductor()    # optional
        24300

    """
    from sage.interfaces.all import magma
    magma.eval("P<x,y,z> := ProjectivePlane(RationalField());")
    cmd = 'aInvariants(MinimalModel(EllipticCurve(Curve(Scheme(P, %s)),P!%s)));'%(F, P)
    s = magma.eval(cmd)
    return EllipticCurve(rings.RationalField(), eval(s))

