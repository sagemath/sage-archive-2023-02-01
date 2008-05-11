r"""
Elliptic curves over a general ring

Elliptic curves are always represented by `Weierstass Models' with
five coefficients $[a_1,a_2,a_3,a_4,a_6]$ in standard notation.  In
Magma, `Weierstrass Model' means a model with a1=a2=a3=0, which is
called `Short Weierstrass Model' in Sage; these only exist in
characteristics other than 2 and 3.

EXAMPLES:
We construct an elliptic curve over an elaborate base ring:
    sage: p = 97; a=1; b=3
    sage: R, u = PolynomialRing(GF(p), 'u').objgen()
    sage: S, v = PolynomialRing(R, 'v').objgen()
    sage: T = S.fraction_field()
    sage: E = EllipticCurve(T, [a, b]); E
    Elliptic Curve defined by y^2  = x^3 + x + 3 over Fraction Field of Univariate Polynomial Ring in v over Univariate Polynomial Ring in u over Finite Field of size 97
    sage: latex(E)
    y^2  = x^3 + x + 3

AUTHORS:
   * William Stein (2005) -- Initial version
   * Robert Bradshaw et al....
   * John Cremona (Jan 2008) -- isomorphisms, automorphisms and twists
                                in all characteristics
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

import math

from sage.rings.all import PolynomialRing

import sage.plot.all as plot

import sage.rings.arith as arith
import sage.rings.all as rings
import sage.rings.number_field as number_field
from sage.rings.all import is_Infinite
import sage.misc.misc as misc
import sage.misc.latex as latex
import sage.modular.modform as modform
import sage.functions.transcendental as transcendental

from sage.categories.morphism import IdentityMorphism
from sage.categories.homset import Hom
from sage.rings.arith import lcm

# Schemes
import sage.schemes.generic.projective_space as projective_space
import sage.schemes.generic.homset as homset

import ell_point
import constructor
import formal_group
import weierstrass_morphism as wm
from constructor import EllipticCurve


factor = arith.factor
sqrt = math.sqrt
exp = math.exp
mul = misc.mul
next_prime = arith.next_prime

oo = rings.infinity       # infinity
O = rings.O         # big oh

import sage.schemes.plane_curves.projective_curve as plane_curve

def is_EllipticCurve(x):
    """
    EXAMPLES:
        sage: E = EllipticCurve([1,2,3/4,7,19])
        sage: is_EllipticCurve(E)
        True
        sage: is_EllipticCurve(0)
        False
    """
    return isinstance(x, EllipticCurve_generic)

class EllipticCurve_generic(plane_curve.ProjectiveCurve_generic):
    """
    Elliptic curve over a generic base ring.

    EXAMPLES:
        sage: E = EllipticCurve([1,2,3/4,7,19]); E
        Elliptic Curve defined by y^2 + x*y + 3/4*y = x^3 + 2*x^2 + 7*x + 19 over Rational Field
        sage: loads(E.dumps()) == E
        True
        sage: E = EllipticCurve([1,3])
        sage: P = E([-1,1,1])
        sage: -5*P
        (179051/80089 : -91814227/22665187 : 1)
    """
    def __init__(self, ainvs, extra=None):
        """
        Constructor from [a1,a2,a3,a4,a6] or [a4,a6] (see
        constructor.py for more variants)

        EXAMPLES:
        sage: E = EllipticCurve([1,2,3,4,5]); E
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field
        sage: E = EllipticCurve(GF(7),[1,2,3,4,5]); E
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 7

        Constructor from [a4,a6] sets a1=a2=a3=0:

        sage: EllipticCurve([4,5]).ainvs()
        [0, 0, 0, 4, 5]

        Base need not be a field:

        sage: EllipticCurve(IntegerModRing(91),[1,2,3,4,5])
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Ring of integers modulo 91

        """
        if extra != None:   # possibility of two arguments
            K, ainvs = ainvs, extra
        else:
            K = ainvs[0].parent()
        assert len(ainvs) == 2 or len(ainvs) == 5
        self.__base_ring = K
        ainvs = [K(x) for x in ainvs]
        if len(ainvs) == 2:
            ainvs = [K(0),K(0),K(0)] + ainvs
        self.__ainvs = ainvs
        if self.discriminant() == 0:
            raise ArithmeticError, \
                  "Invariants %s define a singular curve."%ainvs
        PP = projective_space.ProjectiveSpace(2, K, names='xyz');
        x, y, z = PP.coordinate_ring().gens()
        a1, a2, a3, a4, a6 = ainvs
        f = y**2*z + (a1*x + a3*z)*y*z \
            - (x**3 + a2*x**2*z + a4*x*z**2 + a6*z**3)
        plane_curve.ProjectiveCurve_generic.__init__(self, PP, f)
        # TODO: cleanup, are these two point classes redundant?
        if K.is_field():
            self._point_morphism_class = self._point_class = ell_point.EllipticCurvePoint_field
        else:
            self._point_morphism_class = self._point_class = ell_point.EllipticCurvePoint

    def _defining_params_(self):
        """
        Internal function.  Returns a tuple of the base ring of this
        elliptic curve and its a-invariants, from which it can be
        reconstructed.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E._defining_params_()
            (Rational Field, [0, 0, 0, 1, 1])
            sage: EllipticCurve(*E._defining_params_()) == E
            True
        """
        return (self.__base_ring, self.__ainvs)

    def _repr_(self):
        """
        String representation of elliptic curve.

        EXAMPLES:
            sage: EllipticCurve([1,2,3,4,5])
            Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field

            sage: R.<x> = QQ['x']
            sage: K.<a> = NumberField(x^3-17)
            sage: EllipticCurve([a^2-3, -2/3*a + 3])
            Elliptic Curve defined by y^2  = x^3 + (a^2-3)*x + (-2/3*a+3) over Number Field in a
            with defining polynomial x^3 - 17
        """
        #return "Elliptic Curve with a-invariants %s over %s"%(self.ainvs(), self.base_ring())
        b = self.ainvs()
        #return "y^2 + %s*x*y + %s*y = x^3 + %s*x^2 + %s*x + %s"%\
        #       (a[0], a[2], a[1], a[3], a[4])
        a = [z._coeff_repr() for z in b]
        s = "Elliptic Curve defined by "
        s += "y^2 "
        if a[0] == "-1":
            s += "- x*y "
        elif a[0] == '1':
            s += "+ x*y "
        elif b[0]:
            s += "+ %s*x*y "%a[0]
        if a[2] == "-1":
            s += " - y"
        elif a[2] == '1':
            s += "+ y"
        elif b[2]:
            s += "+ %s*y"%a[2]
        s += " = x^3 "
        if a[1] == "-1":
            s += "- x^2 "
        elif a[1] == '1':
            s += "+ x^2 "
        elif b[1]:
            s += "+ %s*x^2 "%a[1]
        if a[3] == "-1":
            s += "- x "
        elif a[3] == '1':
            s += "+ x "
        elif b[3]:
            s += "+ %s*x "%a[3]
        if a[4] == '-1':
            s += "-1 "
        elif a[4] == '1':
            s += "+1 "
        elif b[4]:
            s += "+ %s "%a[4]
        s = s.replace("+ -","- ")
        s += "over %s"%self.base_ring()
        return s

    def _latex_(self):
        """
        Internal function.  Returns a latex string for this elliptic
        curve.  Users will normally use latex() instead.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E._latex_()
            'y^2  = x^3 + x +1 '
            sage: latex(E)
            y^2  = x^3 + x +1
        """
        b = self.ainvs()
        a = [z._latex_coeff_repr() for z in b]
        s = "y^2 "
        if a[0] == '-1':
            s += "- xy "
        elif a[0] == '1':
            s += "+ xy "
        elif b[0]:
            s += "+ %sxy "%a[0]
        if a[2] == '-1':
            s += " - y"
        elif a[2] == '1':
            s += "+ y"
        elif b[2]:
            s += "+ %sy"%a[2]
        s += " = x^3 "
        if a[1] == '-1':
            s += "- x^2 "
        elif a[1] == '1':
            s += "+ x^2 "
        elif b[1]:
            s += "+ %sx^2 "%a[1]
        if a[3] == '-1':
            s += "- x "
        elif a[3] == '1':
            s += "+ x "
        elif b[3]:
            s += "+ %sx "%a[3]
        if a[4] == '-1':
            s += "-1 "
        elif a[4] == '1':
            s += "+1 "
        elif b[4]:
            s += "+ %s "%a[4]
        s = s.replace("+ -","- ")
        return s

    def _pari_init_(self):
        """
        Internal function.  Returns a string to initialize this
        elliptic curve in the pari system.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E._pari_init_()
            'ellinit([0/1,0/1,0/1,1/1,1/1])'
        """
        return 'ellinit([%s])'%(','.join([x._pari_init_() for x in self.ainvs()]))

    def _magma_init_(self):
        """
        Internal function.  Returns a string to initialize this
        elliptic curve in the Magma subsystem.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E._magma_init_()
            'EllipticCurve([0/1,0/1,0/1,1/1,1/1])'
        """
        return 'EllipticCurve([%s])'%(','.join([x._magma_init_() for x in self.ainvs()]))

    def _symbolic_(self, SR):
        r"""
        Many elliptic curves can be converted into a symbolic expression
        using the \code{symbolic_expression} command.

        EXAMPLES:
        We find a torsion point on 11a.
            sage: E = EllipticCurve('11a')
            sage: E.torsion_subgroup().gens()
            ((5 : 5 : 1),)

        We find the corresponding symbolic equality:
            sage: eqn = symbolic_expression(E); eqn
            y^2 + y == x^3 - x^2 - 10*x - 20
            sage: print eqn
                                      2        3    2
                                     y  + y == x  - x  - 10 x - 20

        We verify that the given point is on the curve:
            sage: eqn(x=5,y=5)
            30 == 30
            sage: bool(eqn(x=5,y=5))
            True

        We create a single expression:
            sage: F = eqn.lhs() - eqn.rhs(); print F
                                      2        3    2
                                     y  + y - x  + x  + 10 x + 20
            sage: y = var('y')
            sage: print F.solve(y)
            [
                                  3      2
                        - sqrt(4 x  - 4 x  - 40 x - 79) - 1
                    y == -----------------------------------
                                         2,
                                 3      2
                         sqrt(4 x  - 4 x  - 40 x - 79) - 1
                     y == ---------------------------------
                                         2
            ]

        You can also solve for x in terms of y, but the result is horrendous.
        Continuing with the above example, we can explicitly find points
        over random fields by substituting in values for x:

            sage: v = F.solve(y)[0].rhs()
            sage: print v
                                            3      2
                                  - sqrt(4 x  - 4 x  - 40 x - 79) - 1
                                  -----------------------------------
                                                   2
            sage: v(3)
            (-sqrt(127)*I - 1)/2
            sage: v(7)
            (-sqrt(817) - 1)/2
            sage: v(-7)
            (-sqrt(1367)*I - 1)/2
            sage: v(sqrt(2))
            (-sqrt(-32*sqrt(2) - 87) - 1)/2

        We can even do arithmetic with them, as follows:
            sage: E2 = E.change_ring(SR); E2
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Symbolic Ring
            sage: P = E2.point((3, v(3), 1), check=False)
            sage: P
            (3 : (-sqrt(127)*I - 1)/2 : 1)
            sage: P + P
            (-756/127 : (sqrt(127)*I + 1)/2 + 12507*I/(127*sqrt(127)) - 1 : 1)

        We can even throw in a transcendental:
            sage: w = E2.point((pi,v(pi),1), check=False); w
            (pi : (-sqrt(4*pi^3 - 4*pi^2 - 40*pi - 79) - 1)/2 : 1)
            sage: 2*w
            ((3*pi^2 - 2*pi - 10)^2/(4*pi^3 - 4*pi^2 - 40*pi - 79) - 2*pi + 1 : (sqrt(4*pi^3 - 4*pi^2 - 40*pi - 79) + 1)/2 - (3*pi^2 - 2*pi - 10)*(-(3*pi^2 - 2*pi - 10)^2/(4*pi^3 - 4*pi^2 - 40*pi - 79) + 3*pi - 1)/sqrt(4*pi^3 - 4*pi^2 - 40*pi - 79) - 1 : 1)
        """
        a = [SR(x) for x in self.a_invariants()]
        x, y = SR.var('x, y')
        return y**2 + a[0]*x*y + a[2]*y == x**3 + a[1]*x**2 + a[3]*x + a[4]

    def __cmp__(self, other):
        """
        Standard comparison function for elliptic curves, to allow
        sorting and equality testing.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: F=EllipticCurve(QQ,[0,0,0,1,1])
            sage: E==F
            True
        """
        if not isinstance(other, EllipticCurve_generic):
            return -1
        t = cmp(self.base_ring(), other.base_ring())
        if t:
            return t
        return cmp(self.ainvs(), other.ainvs())

    def __contains__(self, P):
        """
        Returns True if and only if P defines is a point on the
        elliptic curve.  P just has to be something that can be
        coerced to a point.

        EXAMPLES:
            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: (0,0) in E
            True
            sage: (1,3) in E
            False
            sage: E = EllipticCurve([GF(7)(0), 1])
            sage: [0,0] in E
            False
            sage: [0,8] in E
            True
            sage: P = E(0,8)
            sage: P
            (0 : 1 : 1)
            sage: P in E
            True
        """
        if not isinstance(P, ell_point.EllipticCurvePoint):
            try:
                P = self(P)
            except TypeError:
                return False
        if P.curve() == self:
            return True
        x, y, a = P[0], P[1], self.ainvs()
        return y**2 + a[0]*x*y + a[2]*y == x**3 + a[1]*x**2 + a[3]*x + a[4]

    def __call__(self, *args, **kwds):
        """
        EXAMPLES:
            sage: E = EllipticCurve([0, 0, 1, -1, 0])

        The point at infinity, which is the 0 element of the group:
            sage: E(0)
            (0 : 1 : 0)

        The origin is a point on our curve:
            sage: P = E([0,0])
            sage: P
            (0 : 0 : 1)

        The curve associated to a point:
            sage: P.curve()
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        Points can be specified by given a 2-tuple or 3-tuple
            sage: E([0,0])
            (0 : 0 : 1)
            sage: E([0,1,0])
            (0 : 1 : 0)

        Over a field, points are normalized so the right-most nonzero entry is 1:
            sage: E(105, -69, 125)
            (21/25 : -69/125 : 1)

        We create points on an elliptic curve over a prime finite field.
            sage: E = EllipticCurve([GF(7)(0), 1])
            sage: E([2,3])
            (2 : 3 : 1)
            sage: E([0,0])
            Traceback (most recent call last):
            ...
            TypeError: coordinates [0, 0, 1] do not define a point on Elliptic Curve defined by y^2  = x^3 +1 over Finite Field of size 7

        We create a point on an elliptic curve over a number field.
            sage: x = polygen(RationalField())
            sage: K = NumberField(x**3 + x + 1, 'a'); a = K.gen()
            sage: E = EllipticCurve([a,a])
            sage: E
            Elliptic Curve defined by y^2  = x^3 + a*x + a over Number Field in a with defining polynomial x^3 + x + 1
            sage: E = EllipticCurve([K(1),1])
            sage: E
            Elliptic Curve defined by y^2  = x^3 + x +1 over Number Field in a with defining polynomial x^3 + x + 1
            sage: P = E([a,0,1])
            sage: P
            (a : 0 : 1)
            sage: P+P
            (0 : 1 : 0)

        Another example involving p-adics:
            sage: E = EllipticCurve('37a1')
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: R = pAdicField(3,20)
            sage: Ep = E.base_extend(R); Ep
            Elliptic Curve defined by y^2 + (1+O(3^20))*y = x^3 + (2+2*3+2*3^2+2*3^3+2*3^4+2*3^5+2*3^6+2*3^7+2*3^8+2*3^9+2*3^10+2*3^11+2*3^12+2*3^13+2*3^14+2*3^15+2*3^16+2*3^17+2*3^18+2*3^19+O(3^20))*x over 3-adic Field with capped relative precision 20
            sage: Ep(P)
            (0 : 0 : 1 + O(3^20))
        """
        if len(args) == 1 and args[0] == 0:
            R = self.base_ring()
            return self.point([R(0),R(1),R(0)], check=False)
        if isinstance(args[0],
              (ell_point.EllipticCurvePoint_field, ell_point.EllipticCurvePoint)):
            args = tuple(args[0])
        return plane_curve.ProjectiveCurve_generic.__call__(self, *args, **kwds)

    def lift_x(self, x, all=False):
        """
        Given the x-coordinate of a point on the curve, use the defining
        polynomial to find all affine points on this curve with the
        given x-coordinate.

        EXAMPLES:
            sage: E = EllipticCurve('37a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: E.lift_x(1)
            (1 : 0 : 1)
            sage: E.lift_x(2)
            (2 : 2 : 1)
            sage: E.lift_x(1/4, all=True)
            [(1/4 : -3/8 : 1), (1/4 : -5/8 : 1)]

        There are no rational points with x-cordinate 3.
            sage: E.lift_x(3)
            Traceback (most recent call last):
            ...
            ValueError: No point with x-coordinate 3 on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        However, there are two such points in $E(\R)$:
            sage: E.change_ring(RR).lift_x(3, all=True)
            [(3.00000000000000 : 4.42442890089805 : 1.00000000000000), (3.00000000000000 : -5.42442890089805 : 1.00000000000000)]

        And of course it always works in $E(\C)$:
            sage: E.change_ring(RR).lift_x(.5, all=True)
            []
            sage: E.change_ring(CC).lift_x(.5)
            (0.500000000000000 : -0.500000000000000 + 0.353553390593274*I : 1.00000000000000)

        We can perform these operations over finite fields too:
            sage: E = E.change_ring(GF(17)); E
            Elliptic Curve defined by y^2 + y = x^3 + 16*x over Finite Field of size 17
            sage: E.lift_x(7)
            (7 : 11 : 1)
            sage: E.lift_x(3)
            Traceback (most recent call last):
            ...
            ValueError: No point with x-coordinate 3 on Elliptic Curve defined by y^2 + y = x^3 + 16*x over Finite Field of size 17

        Note that there is only one lift with x-coordinate 10 in $E(\F_{17})$.
            sage: E.lift_x(10, all=True)
            [(10 : 8 : 1)]

        We can lift over more exotic rings too.
            sage: E = EllipticCurve('37a');
            sage: E.lift_x(pAdicField(17, 5)(6))
            (6 + O(17^5) : 2 + 16*17 + 16*17^2 + 16*17^3 + 16*17^4 + O(17^5) : 1 + O(17^5))
            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: E.lift_x(1+t)
            (1 + t : 2*t - t^2 + 5*t^3 - 21*t^4 + O(t^5) : 1)
            sage: K.<a> = GF(16)
            sage: E = E.change_ring(K)
            sage: E.lift_x(a^3)
            (a^3 : a^3 + a : 1)


        AUTHOR: Robert Bradshaw, 2007-04-24

        TEST:
            sage: E = EllipticCurve('37a').short_weierstrass_model().change_ring(GF(17))
            sage: E.lift_x(3, all=True)
            []
            sage: E.lift_x(7, all=True)
            [(7 : 3 : 1), (7 : 14 : 1)]
        """
        a1, a2, a3, a4, a6 = self.ainvs()
        f = ((x + a2) * x + a4) * x + a6
        K = self.base_ring()
        x += K(0)
        one = x.parent()(1)
        if a1.is_zero() and a3.is_zero():
            if f.is_square():
                if all:
                    ys = f.sqrt(all=True)
                    return [self.point([x, y, one], check=False) for y in ys]
                else:
                    return self.point([x, f.sqrt(), one], check=False)
        else:
            b = (a1*x + a3)
            D = b*b + 4*f
            if K.characteristic() == 2:
                R = PolynomialRing(K, 'y')
                F = R([-f,b,1])
                ys = F.roots()
                if all:
                    return [self.point([x, y[0], one], check=False) for y in ys]
                elif len(ys) > 0:
                    return self.point([x, ys[0][0], one], check=False)
            elif D.is_square():
                if all:
                    return [self.point([x, (-b+d)/2, one], check=False) for d in D.sqrt(all=True)]
                else:
                    return self.point([x, (-b+D.sqrt())/2, one], check=False)
        if all:
            return []
        else:
            raise ValueError, "No point with x-coordinate %s on %s"%(x, self)

    def _homset_class(self, *args, **kwds):
        """
        Internal function.  Returns the (abstract) group of points on
        this elliptic curve over a ring.

        EXAMPLES:
            sage: E=EllipticCurve(GF(5),[1,1])
            sage: E._homset_class(GF(5^10,'a'),GF(5))
            Abelian group of points on Finite Field in a of size 5^10
        """
        return homset.SchemeHomsetModule_abelian_variety_coordinates_field(*args, **kwds)

    def __getitem__(self, n):
        """
        Placeholder for standard indexing function.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E[2]
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented.
        """
        raise NotImplementedError, "not implemented."

    def __is_over_RationalField(self):
        """
        Internal function.  Returns true iff the base ring of this
        elliptic curve is the field of rational numbers.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E._EllipticCurve_generic__is_over_RationalField()
            True
            sage: E=EllipticCurve(GF(5),[1,1])
            sage: E._EllipticCurve_generic__is_over_RationalField()
            False
        """
        return isinstance(self.base_ring(), rings.RationalField)

    def change_ring(self, R):
        """
        Return the elliptic curve defined by coercing the a-invariants
        of this elliptic curve into the ring R.

        INPUT:
            R -- ring

        OUTPUT:
            an elliptic curve

        EXAMPLES:
            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: E.change_ring(GF(3))
            Elliptic Curve defined by y^2 + y = x^3 + 2*x over Finite Field of size 3
        """
        return constructor.EllipticCurve(R, [R(a) for a in self.ainvs()])

    def is_on_curve(self, x, y):
        """
        Returns True if the (x,y) is an affine point on this curve.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E.is_on_curve(0,1)
            True
            sage: E.is_on_curve(1,1)
            False
        """
        a = self.ainvs()
        return y**2 +a[0]*x*y + a[2]*y == x**3 + a[1]*x**2 + a[3]*x + a[4]

    def a_invariants(self):
        """
        The a-invariants of this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.a_invariants()
            [1, 2, 3, 4, 5]
            sage: E = EllipticCurve([0,1])
            sage: E
            Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
            sage: E.a_invariants()
            [0, 0, 0, 0, 1]
            sage: E = EllipticCurve([GF(7)(3),5])
            sage: E.a_invariants()
            [0, 0, 0, 3, 5]
        """
        return self.__ainvs

    ainvs = a_invariants

    def b_invariants(self):
        """
        The b-invariants of this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.b_invariants()
            (-4, -20, -79, -21)
            sage: E = EllipticCurve([-4,0])
            sage: E.b_invariants()
            (0, -8, 0, -16)

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.b_invariants()
            (9, 11, 29, 35)
            sage: E.b2()
            9
            sage: E.b4()
            11
            sage: E.b6()
            29
            sage: E.b8()
            35

        ALGORITHM: These are simple functions of the a invariants.

        AUTHOR: William Stein, 2005-04-25
        """
        try:
            return self.__b_invariants
        except AttributeError:
            a1,a2,a3,a4,a6 = tuple(self.ainvs())
            self.__b_invariants = a1*a1 + 4*a2, \
                                  a1*a3 + 2*a4, \
                                  a3**2 + 4*a6, \
                                  a1**2 * a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2
            return self.__b_invariants

    def b2(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.b2()
            9
        """
        try:
            return self.__b_invariants[0]
        except AttributeError:
            return self.b_invariants()[0]

    def b4(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.b4()
            11
        """
        try:
            return self.__b_invariants[1]
        except AttributeError:
            return self.b_invariants()[1]

    def b6(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.b6()
            29
        """
        try:
            return self.__b_invariants[2]
        except AttributeError:
            return self.b_invariants()[2]

    def b8(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.b8()
            35
        """
        try:
            return self.__b_invariants[3]
        except AttributeError:
            return self.b_invariants()[3]

    def c_invariants(self):
        """
        The c-invariants of this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.c_invariants()
            (496, 20008)
            sage: E = EllipticCurve([-4,0])
            sage: E.c_invariants()
            (192, 0)

        ALGORITHM: These are simple functions of the b invariants.

        AUTHOR: William Stein, 2005-04-25
        """
        try:
            return self.__c_invariants
        except AttributeError:
            b2,b4,b6,b8 = self.b_invariants()
            self.__c_invariants = b2**2 - 24*b4,\
                                  -b2**3 + 36*b2*b4 - 216*b6    # note: c6 is wrong in Silverman, but right in Cremona
            return self.__c_invariants

    def c4(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.c4()
            496
        """
        try:
            return self.__c_invariants[0]
        except AttributeError:
            pass
        return self.c_invariants()[0]

    def c6(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.c6()
            20008
        """
        try:
            return self.__c_invariants[1]
        except AttributeError:
            pass
        return self.c_invariants()[1]


    def base_extend(self, R):
        """
        Returns a new curve with the same a-invariants but defined
        over a new ring into which the original's a-invariants may be
        mapped.  R is either a ring into which they may be coerced, or
        a morphism which may be applied to them.

        EXAMPLES:
        sage: E=EllipticCurve(GF(5),[1,1]); E
        Elliptic Curve defined by y^2  = x^3 + x +1 over Finite Field of size 5
        sage: E1=E.base_extend(GF(125,'a')); E1
        Elliptic Curve defined by y^2  = x^3 + x +1 over Finite Field in a of size 5^3
        sage: F2=GF(5^2,'a'); a=F2.gen()
        sage: F4=GF(5^4,'b'); b=F4.gen()
        sage: h=F2.hom([a.charpoly().roots(ring=F4,multiplicities=False)[0]],F4)
        sage: E=EllipticCurve(F2,[1,a]); E
        Elliptic Curve defined by y^2  = x^3 + x + a over Finite Field in a of size 5^2
        sage: E.base_extend(h)
        Elliptic Curve defined by y^2  = x^3 + x + (4*b^3+4*b^2+4*b+3) over Finite Field in b of size 5^4
        """
        return constructor.EllipticCurve([R(a) for a in self.a_invariants()])

    def base_ring(self):
        """
        Returns the base ring of the elliptic curves.

        EXAMPLES:
            sage: E = EllipticCurve(GF(49, 'a'), [3,5])
            sage: E.base_ring()
            Finite Field in a of size 7^2

            sage: E = EllipticCurve([1,1])
            sage: E.base_ring()
            Rational Field

        """
       # NOT supported yet
       #     sage: E = EllipticCurve(Z, [3,5])
       #     sage: E.base_ring()
       #     Integer Ring
        return self.__base_ring

    base_field = base_ring

    def a1(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,6])
            sage: E.a1()
            1
        """
        return self.__ainvs[0]

    def a2(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,6])
            sage: E.a2()
            2
        """
        return self.__ainvs[1]

    def a3(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,6])
            sage: E.a3()
            3
        """
        return self.__ainvs[2]

    def a4(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,6])
            sage: E.a4()
            4
        """
        return self.__ainvs[3]

    def a6(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,6])
            sage: E.a6()
            6
        """
        return self.__ainvs[4]

    def gens(self):
        """
        Placeholder function to return generators of an elliptic
        curve: derived classes such as EllipticCurve_rational_field
        implement this functionality

        EXAMPLES:
            sage: R.<a1,a2,a3,a4,a6>=QQ[]
            sage: E=EllipticCurve([a1,a2,a3,a4,a6])
            sage: E.gens()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented.
            sage: E=EllipticCurve(QQ,[1,1])
            sage: E.gens()
            [(0 : 1 : 1)]
        """
        raise NotImplementedError, "not implemented."

    def gen(self, i):
        """
        Function returning the i'th generator of this elliptic curve.
        Relies on gens() being implemented.

        EXAMPLES:
            sage: R.<a1,a2,a3,a4,a6>=QQ[]
            sage: E=EllipticCurve([a1,a2,a3,a4,a6])
            sage: E.gen(0)
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented.
        """
        return self.gens()[i]


    # Twists: rewritten by John Cremona as follows:
    #
    # Quadratic twist allowed except when char=2, j=0
    # Quartic twist allowed only if j=1728!=0 (so char!=2,3)
    # Sextic  twist allowed only if j=0!=1728 (so char!=2,3)
    #
    # More complicated twists exist in theory for char=2,3 and
    # j=0=1728, but I have never worked them out or seen them used!
    #

    def quadratic_twist(self, D):
        """
        Return the quadratic twist of this curve by D, which must be nonzero except in characteristic 2.

        In characteristic!=2, D must be nonzero, and the twist is
        isomorphic to self after adjoining sqrt(D) to the base

        In characteristic==2, D is arbitrary, and the twist is
        isomorphic to self after adjoining a root of $x^2+x+D$ to the
        base

        In characteristics 2 when j==0 this is not implemented (the
        twists are more complicated than quadratic!)

        EXAMPLES:
            sage: E = EllipticCurve([GF(1103)(1), 0, 0, 107, 340]); E
            Elliptic Curve defined by y^2 + x*y  = x^3 + 107*x + 340 over Finite Field of size 1103
            sage: F=E.quadratic_twist(-1); F
            Elliptic Curve defined by y^2  = x^3 + 1102*x^2 + 609*x + 300 over Finite Field of size 1103
            sage: E.is_isomorphic(F)
            False
            sage: E.is_isomorphic(F,GF(1103^2,'a'))
            True

            A characteristic 2 example:

            sage: E=EllipticCurve(GF(2),[1,0,1,1,1])
            sage: E1=E.quadratic_twist(1)
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1,GF(4,'a'))
            True

        """
        K=self.base_ring()
        char=K.characteristic()
        D=K(D)

        if char!=2:
            if D.is_zero():
                raise ValueError, "quadratic twist requires a nonzero argument when characteristic is not 2"
            b2,b4,b6,b8=self.b_invariants()
            # E is isomorphic to  [0,b2,0,8*b4,16*b6]
            return EllipticCurve(K,[0,b2*D,0,8*b4*D**2,16*b6*D**3])

        # now char==2
        if self.j_invariant() !=0: # iff a1!=0
            a1,a2,a3,a4,a6=self.ainvs()
            E0=self.change_weierstrass_model(a1,a3/a1,0,(a1**2*a4+a3**2)/a1**3)
            # which has the form = [1,A2,0,0,A6]
            assert E0.a1()==K(1)
            assert E0.a3()==K(0)
            assert E0.a4()==K(0)
            return EllipticCurve(K,[1,E0.a2()+D,0,0,E0.a6()])
        else:
            raise ValueError, "Quadratic twist not implemented in char 2 when j=0"

    def quartic_twist(self, D):
        """
        Return the quartic twist of this curve by D, which must be nonzero.

        The characteristic must not be 2 or 3 and the j-invariant must be 1728

        EXAMPLES:
        sage: E=EllipticCurve(GF(13)(1728)); E
        Elliptic Curve defined by y^2  = x^3 + x over Finite Field of size 13
        sage: E1=E.quartic_twist(2); E1
        Elliptic Curve defined by y^2  = x^3 + 5*x over Finite Field of size 13
        sage: E.is_isomorphic(E1)
        False
        sage: E.is_isomorphic(E1,GF(13^2,'a'))
        False
        sage: E.is_isomorphic(E1,GF(13^4,'a'))
        True
        """
        K=self.base_ring()
        char=K.characteristic()
        D=K(D)

        if char==2 or char==3:
            raise ValueError, "Quartic twist not defined in chars 2,3"

        if self.j_invariant() !=K(1728):
            raise ValueError, "Quartic twist not defined when j!=1728"

        if D.is_zero():
            raise ValueError, "quartic twist requires a nonzero argument"

        c4,c6=self.c_invariants()
        # E is isomorphic to  [0,0,0,-27*c4,0]
        assert c6==0
        return EllipticCurve(K,[0,0,0,-27*c4*D,0])

    def sextic_twist(self, D):
        """
        Return the sextic twist of this curve by D, which must be nonzero.

        The characteristic must not be 2 or 3 and the j-invariant must be 0

        EXAMPLES:
        sage: E=EllipticCurve(GF(13)(0)); E
        Elliptic Curve defined by y^2  = x^3 +1 over Finite Field of size 13
        sage: E1=E.sextic_twist(2); E1
        Elliptic Curve defined by y^2  = x^3 + 11 over Finite Field of size 13
        sage: E.is_isomorphic(E1)
        False
        sage: E.is_isomorphic(E1,GF(13^2,'a'))
        False
        sage: E.is_isomorphic(E1,GF(13^4,'a'))
        False
        sage: E.is_isomorphic(E1,GF(13^6,'a'))
        True
        """
        K=self.base_ring()
        char=K.characteristic()
        D=K(D)

        if char==2 or char==3:
            raise ValueError, "Sextic twist not defined in chars 2,3"

        if self.j_invariant() !=K(0):
            raise ValueError, "Sextic twist not defined when j!=0"

        if D.is_zero():
            raise ValueError, "Sextic twist requires a nonzero argument"

        c4,c6=self.c_invariants()
        # E is isomorphic to  [0,0,0,0,-54*c6]
        assert c4==0
        return EllipticCurve(K,[0,0,0,0,-54*c6*D])

    def rst_transform(self, r, s, t):
        """
        Transforms the elliptic curve using the unimodular (u=1)
        transform with standard parameters [r,s,t].  This is just a
        special case of change_weierstrass_model().

        Returns the transformed curve.
        EXAMPLES:
        sage: R.<r,s,t>=QQ[]
        sage: E=EllipticCurve([1,2,3,4,5])
        sage: E.rst_transform(r,s,t)
        Elliptic Curve defined by y^2 + (2*s+1)*x*y + (r+2*t+3)*y = x^3 + (-s^2+3*r-s+2)*x^2 + (3*r^2-r*s-2*s*t+4*r-3*s-t+4)*x + (r^3+2*r^2-r*t-t^2+4*r-3*t+5) over Multivariate Polynomial Ring in r, s, t over Rational Field

        """
        return self.change_weierstrass_model(1,r,s,t)

    def scale_curve(self, u):
        """
        Transforms the elliptic curve using scale factor $u$,
        i.e. multiplies $c_i$ by $u^i$.  This is  another
        special case of change_weierstrass_model().

        Returns the transformed curve.

        EXAMPLES:
        sage: K=Frac(PolynomialRing(QQ,'u'))
        sage: u=K.gen()
        sage: E=EllipticCurve([1,2,3,4,5])
        sage: E.scale_curve(u)
        Elliptic Curve defined by y^2 + u*x*y + 3*u^3*y = x^3 + 2*u^2*x^2 + 4*u^4*x + 5*u^6 over Fraction Field of Univariate Polynomial Ring in u over Rational Field

        """
        if isinstance(u, (int,long)):
            u=self.base_ring()(u)       # because otherwise 1/u would round!
        return self.change_weierstrass_model(1/u,0,0,0)

    def discriminant(self):
        """
        Returns the discriminant of this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: E.discriminant()
            37
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.discriminant()
            -161051

            sage: E = EllipticCurve([GF(7)(2),1])
            sage: E.discriminant()
            1
        """
        try:
            return self.__discriminant
        except AttributeError:
            b2, b4, b6, b8 = self.b_invariants()
            self.__discriminant = -b2**2*b8 - 8*b4**3 - 27*b6**2 + 9*b2*b4*b6
            return self.__discriminant

    def j_invariant(self):
        """
        Returns the j-invariant of this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: E.j_invariant()
            110592/37
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.j_invariant()
            -122023936/161051
            sage: E = EllipticCurve([-4,0])
            sage: E.j_invariant()
            1728

            sage: E = EllipticCurve([GF(7)(2),1])
            sage: E.j_invariant()
            1


        """
        try:
            return self.__j_invariant
        except AttributeError:
            c4, _ = self.c_invariants()
            self.__j_invariant = c4**3 / self.discriminant()
            return self.__j_invariant


##     def pseudo_torsion_polynomial(self, n, x=None, cache=None):
##         r"""
##         Returns the n-th torsion polynomial (division polynomial), without
##         the 2-torsion factor if n is even, as a polynomial in $x$.

##         These are the polynomials $g_n$ defined in Mazur/Tate (``The p-adic
##         sigma function''), but with the sign flipped for even $n$, so that
##         the leading coefficient is always positive.

##         The full torsion polynomials may be recovered as follows:
##         \begin{itemize}
##         \item $\psi_n = g_n$ for odd $n$.
##         \item $\psi_n = (2y + a_1 x + a_3) g_n$ for even $n$.
##         \end{itemize}

##         Note that the $g_n$'s are always polynomials in $x$, whereas the
##         $\psi_n$'s require the appearance of a $y$.

##         SEE ALSO:
##             -- torsion_polynomial()
##             -- multiple_x_numerator()
##             -- multiple_x_denominator()

##         INPUT:
##             n -- positive integer, or the special values -1 and -2 which
##                  mean $B_6 = (2y + a_1 x + a_3)^2$ and $B_6^2$ respectively
##                  (in the notation of Mazur/Tate).
##             x -- optional ring element to use as the "x" variable. If x
##                  is None, then a new polynomial ring will be constructed over
##                  the base ring of the elliptic curve, and its generator will
##                  be used as x. Note that x does not need to be a generator of
##                  a polynomial ring; any ring element is ok. This permits fast
##                  calculation of the torsion polynomial *evaluated* on any
##                  element of a ring.
##             cache -- optional dictionary, with integer keys. If the key m
##                  is in cache, then cache[m] is assumed to be the value of
##                  pseudo_torsion_polynomial(m) for the supplied x. New entries
##                  will be added to the cache as they are computed.

##         ALGORITHM:
##             -- Recursion described in Mazur/Tate. The recursive formulae are
##             evaluated $O((log n)^2)$ times.

##         TODO:
##             -- for better unity of code, it might be good to make the regular
##             torsion_polynomial() function use this as a subroutine.

##         AUTHORS:
##             -- David Harvey (2006-09-24)

##         EXAMPLES:
##            sage: E = EllipticCurve("37a")
##            sage: E.pseudo_torsion_polynomial(1)
##            1
##            sage: E.pseudo_torsion_polynomial(2)
##            1
##            sage: E.pseudo_torsion_polynomial(3)
##            3*x^4 - 6*x^2 + 3*x - 1
##            sage: E.pseudo_torsion_polynomial(4)
##            2*x^6 - 10*x^4 + 10*x^3 - 10*x^2 + 2*x + 1
##            sage: E.pseudo_torsion_polynomial(5)
##            5*x^12 - 62*x^10 + 95*x^9 - 105*x^8 - 60*x^7 + 285*x^6 - 174*x^5 - 5*x^4 - 5*x^3 + 35*x^2 - 15*x + 2
##            sage: E.pseudo_torsion_polynomial(6)
##            3*x^16 - 72*x^14 + 168*x^13 - 364*x^12 + 1120*x^10 - 1144*x^9 + 300*x^8 - 540*x^7 + 1120*x^6 - 588*x^5 - 133*x^4 + 252*x^3 - 114*x^2 + 22*x - 1
##            sage: E.pseudo_torsion_polynomial(7)
##            7*x^24 - 308*x^22 + 986*x^21 - 2954*x^20 + 28*x^19 + 17171*x^18 - 23142*x^17 + 511*x^16 - 5012*x^15 + 43804*x^14 - 7140*x^13 - 96950*x^12 + 111356*x^11 - 19516*x^10 - 49707*x^9 + 40054*x^8 - 124*x^7 - 18382*x^6 + 13342*x^5 - 4816*x^4 + 1099*x^3 - 210*x^2 + 35*x - 3
##            sage: E.pseudo_torsion_polynomial(8)
##            4*x^30 - 292*x^28 + 1252*x^27 - 5436*x^26 + 2340*x^25 + 39834*x^24 - 79560*x^23 + 51432*x^22 - 142896*x^21 + 451596*x^20 - 212040*x^19 - 1005316*x^18 + 1726416*x^17 - 671160*x^16 - 954924*x^15 + 1119552*x^14 + 313308*x^13 - 1502818*x^12 + 1189908*x^11 - 160152*x^10 - 399176*x^9 + 386142*x^8 - 220128*x^7 + 99558*x^6 - 33528*x^5 + 6042*x^4 + 310*x^3 - 406*x^2 + 78*x - 5

##            sage: E.pseudo_torsion_polynomial(18) % E.pseudo_torsion_polynomial(6) == 0
##            True

##           An example to illustrate the relationship with torsion points.
##            sage: F = GF(11)
##            sage: E = EllipticCurve(F, [0, 2]); E
##            Elliptic Curve defined by y^2  = x^3 + 2 over Finite Field of size 11
##            sage: f = E.pseudo_torsion_polynomial(5); f
##            5*x^12 + x^9 + 8*x^6 + 4*x^3 + 7
##            sage: f.factor()
##            (5) * (x^2 + 5) * (x^2 + 2*x + 5) * (x^2 + 5*x + 7) * (x^2 + 7*x + 7) * (x^2 + 9*x + 5) * (x^2 + 10*x + 7)

##           This indicates that the x-coordinates of all the 5-torsion points
##           of $E$ are in $GF(11^2)$, and therefore the y-coordinates are in
##           $GF(11^4)$.

##            sage: K = GF(11^4, 'a')
##            sage: X = E.change_ring(K)
##            sage: f = X.pseudo_torsion_polynomial(5)
##            sage: x_coords = [root for (root, _) in f.roots()]; x_coords
##            [10*a^3 + 4*a^2 + 5*a + 6,
##             9*a^3 + 8*a^2 + 10*a + 8,
##             8*a^3 + a^2 + 4*a + 10,
##             8*a^3 + a^2 + 4*a + 8,
##             8*a^3 + a^2 + 4*a + 4,
##             6*a^3 + 9*a^2 + 3*a + 4,
##             5*a^3 + 2*a^2 + 8*a + 7,
##             3*a^3 + 10*a^2 + 7*a + 8,
##             3*a^3 + 10*a^2 + 7*a + 3,
##             3*a^3 + 10*a^2 + 7*a + 1,
##             2*a^3 + 3*a^2 + a + 7,
##             a^3 + 7*a^2 + 6*a]

##           Now we check that these are exactly the x coordinates of the
##           5-torsion points of E.
##            sage: for x in x_coords:
##            ...       y = (x**3 + 2).square_root()
##            ...       P = X([x, y])
##            ...       assert P.order(disable_warning=True) == 5

##           todo: need to show an example where the 2-torsion is missing

##         """
##         if cache is None:
##             cache = {}
##         else:
##             try:
##                 return cache[n]
##             except KeyError:
##                 pass

##         if x is None:
##             x = rings.PolynomialRing(self.base_ring(), 'x').gen()

##         b2, b4, b6, b8 = self.b_invariants()

##         n = int(n)
##         if n <= 4:
##             if n == -1:
##                 answer = 4*x**3 + b2*x**2 + 2*b4*x + b6
##             elif n == -2:
##                 answer = self.pseudo_torsion_polynomial(-1, x, cache) ** 2
##             elif n == 1 or n == 2:
##                 answer = x.parent()(1)
##             elif n == 3:
##                 answer = 3*x**4 + b2*x**3 + 3*b4*x**2 + 3*b6*x + b8
##             elif n == 4:
##                 answer = -self.pseudo_torsion_polynomial(-2, x, cache) + \
##                          (6*x**2 + b2*x + b4) * \
##                          self.pseudo_torsion_polynomial(3, x, cache)
##             else:
##                 raise ValueError, "n must be a positive integer (or -1 or -2)"
##         else:
##             if n % 2 == 0:
##                 m = (n-2) // 2
##                 g_mplus3 = self.pseudo_torsion_polynomial(m+3, x, cache)
##                 g_mplus2 = self.pseudo_torsion_polynomial(m+2, x, cache)
##                 g_mplus1 = self.pseudo_torsion_polynomial(m+1, x, cache)
##                 g_m      = self.pseudo_torsion_polynomial(m,   x, cache)
##                 g_mless1 = self.pseudo_torsion_polynomial(m-1, x, cache)
##                 answer = g_mplus1 * \
##                          (g_mplus3 * g_m**2 - g_mless1 * g_mplus2**2)
##             else:
##                 m = (n-1) // 2
##                 g_mplus2 = self.pseudo_torsion_polynomial(m+2, x, cache)
##                 g_mplus1 = self.pseudo_torsion_polynomial(m+1, x, cache)
##                 g_m      = self.pseudo_torsion_polynomial(m,   x, cache)
##                 g_mless1 = self.pseudo_torsion_polynomial(m-1, x, cache)
##                 B6_sqr   = self.pseudo_torsion_polynomial(-2, x, cache)
##                 if m % 2 == 0:
##                     answer = B6_sqr * g_mplus2 * g_m**3 - \
##                              g_mless1 * g_mplus1**3
##                 else:
##                     answer = g_mplus2 * g_m**3 - \
##                              B6_sqr * g_mless1 * g_mplus1**3

##         cache[n] = answer
##         return answer


##     def multiple_x_numerator(self, n, x=None, cache=None):
##         r"""
##         Returns the numerator of the x-coordinate of the nth multiple of
##         a point, using torsion polynomials (division polynomials).

##         The inputs n, x, cache are as described in pseudo_torsion_polynomial().

##         The result is adjusted to be correct for both even and odd n.

##         WARNING:
##           -- There may of course be cancellation between the numerator and
##           the denominator (multiple_x_denominator()). Be careful. For more
##           information on how to avoid cancellation, see Christopher Wuthrich's
##           thesis.

##         SEE ALSO:
##           -- multiple_x_denominator()

##         AUTHORS:
##            -- David Harvey (2006-09-24)

##         EXAMPLES:
##           sage: E = EllipticCurve("37a")
##           sage: P = E.gens()[0]
##           sage: x = P[0]

##           sage: (35*P)[0]
##           -804287518035141565236193151/1063198259901027900600665796
##           sage: E.multiple_x_numerator(35, x)
##           -804287518035141565236193151
##           sage: E.multiple_x_denominator(35, x)
##           1063198259901027900600665796

##           sage: (36*P)[0]
##           54202648602164057575419038802/15402543997324146892198790401
##           sage: E.multiple_x_numerator(36, x)
##           54202648602164057575419038802
##           sage: E.multiple_x_denominator(36, x)
##           15402543997324146892198790401

##         An example where cancellation occurs:
##           sage: E = EllipticCurve("88a1")
##           sage: P = E([2,2])   # fixed choice of generator
##           sage: n = E.multiple_x_numerator(11, P[0]); n
##           442446784738847563128068650529343492278651453440
##           sage: d = E.multiple_x_denominator(11, P[0]); d
##           1427247692705959881058285969449495136382746624
##           sage: n/d
##           310
##           sage: 11*P
##           (310 : -5458 : 1)

##         """
##         if cache is None:
##             cache = {}

##         if x is None:
##             x = rings.PolynomialRing(self.base_ring(), 'x').gen()

##         n = int(n)
##         if n < 2:
##             print "n must be at least 2"

##         self.pseudo_torsion_polynomial( -2, x, cache)
##         self.pseudo_torsion_polynomial(n-1, x, cache)
##         self.pseudo_torsion_polynomial(n  , x, cache)
##         self.pseudo_torsion_polynomial(n+1, x, cache)

##         if n % 2 == 0:
##             return x * cache[-1] * cache[n]**2 - cache[n-1] * cache[n+1]
##         else:
##             return x * cache[n]**2 - cache[-1] * cache[n-1] * cache[n+1]


##     def multiple_x_denominator(self, n, x=None, cache=None):
##         r"""
##         Returns the denominator of the x-coordinate of the nth multiple of
##         a point, using torsion polynomials (division polynomials).

##         The inputs n, x, cache are as described in pseudo_torsion_polynomial().

##         The result is adjusted to be correct for both even and odd n.

##         SEE ALSO:
##           -- multiple_x_numerator()

##         TODO: the numerator and denominator versions share a calculation,
##         namely squaring $\psi_n$. Maybe would be good to offer a combined
##         version to make this more efficient.

##         EXAMPLES:
##            -- see multiple_x_numerator()

##         AUTHORS:
##            -- David Harvey (2006-09-24)

##         """
##         if cache is None:
##             cache = {}

##         if x is None:
##             x = rings.PolynomialRing(self.base_ring(), 'x').gen()

##         n = int(n)
##         if n < 2:
##             print "n must be at least 2"

##         self.pseudo_torsion_polynomial(-2, x, cache)
##         self.pseudo_torsion_polynomial(n , x, cache)

##         if n % 2 == 0:
##             return cache[-1] * cache[n]**2
##         else:
##             return cache[n]**2


    def torsion_polynomial(self, n, var='x', i=0):
        """
        Returns the n-th torsion polynomial (a.k.a., division polynomial).

        INPUT:
            n -- non-negative integer
            var -- string; the indeterminate of the polynomial (default: x)
            i -- integer, either 0 (default) or 1.

        OUTPUT:
            Polynomial -- n-th torsion polynomial, which is a polynomial over
                          the base field of the elliptic curve.

        SEE ALSO: full_division_polynomial

        ALIASES: division_polynomial, torsion_polynomial

        EXAMPLES:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: E.division_polynomial(1)
            1
            sage: E.division_polynomial(2)
            4*x^3 - 4*x + 1
            sage: E.division_polynomial(3, 'z')
            3*z^4 - 6*z^2 + 3*z - 1

            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.torsion_polynomial(0)
            0
            sage: E.torsion_polynomial(1)
            1
            sage: E.torsion_polynomial(2)
            4*x^3 - 4*x^2 - 40*x - 79
            sage: E.torsion_polynomial(3)
            3*x^4 - 4*x^3 - 60*x^2 - 237*x - 21
            sage: E.torsion_polynomial(4)
            8*x^9 - 24*x^8 - 464*x^7 - 2758*x^6 + 6636*x^5 + 34356*x^4 + 53510*x^3 + 99714*x^2 + 351024*x + 459859

            sage: E = EllipticCurve([-4,0])
            sage: E.torsion_polynomial(2)
            4*x^3 - 16*x
            sage: E.torsion_polynomial(5)
            5*x^12 - 248*x^10 - 1680*x^8 + 19200*x^6 - 32000*x^4 + 51200*x^2 + 4096
            sage: E.torsion_polynomial(6)
            12*x^19 - 1200*x^17 - 18688*x^15 + 422912*x^13 - 2283520*x^11 + 9134080*x^9 - 27066368*x^7 + 19136512*x^5 + 19660800*x^3 - 3145728*x

        AUTHOR: David Kohel (kohel@maths.usyd.edu.au), 2005-04-25
        """
        n = int(n)
        try:
            return self.__torsion_polynomial[n]
        except AttributeError:
            self.__torsion_polynomial = {}
        except KeyError:
            pass
        E = self; i=int(i)
        if n < 0:
            raise ValueError, "n must be a non-negative integer."
        if i > 1 :
            raise ValueError, "i must be 0 or 1."

        R = rings.PolynomialRing(E.base_ring(), var)
        if i == 1:
            if n == 0:
                f = E.torsion_polynomial(1)
                E.__torsion_polynomial[n] = f
                return f
            else:
                x = R.gen()
                psi2 = E.torsion_polynomial(2)
                if n%2 == 0:
                    f = x * psi2 * (E.torsion_polynomial(n)//psi2)**2 - \
                        E.torsion_polynomial(n+1) * E.torsion_polynomial(n-1)
                    E.__torsion_polynomial[n] = f
                    return f
                else:
                    f = x * E.torsion_polynomial(n)**2 - \
                        (E.torsion_polynomial(n+1)//psi2) * E.torsion_polynomial(n-1)
                    E.__torsion_polynomial[n] = f
                    return f

        else:

            if n == 0: return R(0)
            if n == 1: return R(1)
            x = R.gen()
            b2, b4, b6, b8 = E.b_invariants()
            psi2 = 4*x**3 + b2*x**2 + 2*b4*x + b6
            if n == 2:
                f = psi2
                E.__torsion_polynomial[n] = f; return f
            if n == 3:
                f = 3*x**4 + b2*x**3 + 3*b4*x**2 + 3*b6*x + b8
                E.__torsion_polynomial[n] = f; return f
            if n == 4:
                f = psi2 * (2*x**6 + b2*x**5 + 5*b4*x**4 + 10*b6*x**3 \
                    + 10*b8*x**2  + (b2*b8 - b4*b6)*x + b4*b8 - b6**2)
                E.__torsion_polynomial[n] = f; return f
            if n%2 == 0:
                m = n//2
                if m%2 == 0:
                    f = E.torsion_polynomial(m) * ( \
                        (E.torsion_polynomial(m+2)//psi2) * E.torsion_polynomial(m-1)**2 - \
                        (E.torsion_polynomial(m-2)//psi2) * E.torsion_polynomial(m+1)**2)
                    E.__torsion_polynomial[n] = f; return f
                else:
                    f = psi2 * E.torsion_polynomial(m)*( \
                        E.torsion_polynomial(m+2) * (E.torsion_polynomial(m-1)//psi2)**2 - \
                        E.torsion_polynomial(m-2) * (E.torsion_polynomial(m+1)//psi2)**2)
                    E.__torsion_polynomial[n] = f; return f
            else:
                m = n//2
                if m%2 == 0:
                    f = psi2 * \
                        E.torsion_polynomial(m+2) * (E.torsion_polynomial(m)//psi2)**3 - \
                        E.torsion_polynomial(m-1) * E.torsion_polynomial(m+1)**3
                    E.__torsion_polynomial[n] = f; return f
                else:
                    f = E.torsion_polynomial(m+2) * E.torsion_polynomial(m)**3 - psi2 * \
                        E.torsion_polynomial(m-1)*(E.torsion_polynomial(m+1)//psi2)**3
                    E.__torsion_polynomial[n] = f; return f

    division_polynomial = torsion_polynomial

    def full_division_polynomial(self, m, use_divpoly=True):
        """
        Return the $m$-th bivariate division polynomial in $x$ and
        $y$.  When $m$ is odd this is exactly the same as the usual
        $m$th division polynomial.

        For the usual division polynomial only in $x$, see the
        division_polynomial function.

        INPUT:
            self -- elliptic curve in short Weierstrass form
            m    -- a positive integer
            use_divpoly -- whether to call the division_polynomial
                           function directly in case $m$ is odd.
        OUTPUT:
            a polynomial in two variables $x$, $y$.

        NOTE: The result is cached.

        REFERENCE: Exercise III.3.7 of Silverman AEC 1, 1986, page 105.

        EXAMPLES:
        We create a curve and compute the first two full division
        polynomials.
            sage: E = EllipticCurve([2,3])
            sage: E.full_division_polynomial(1)
            1
            sage: E.full_division_polynomial(2)
            2*y

        Note that for odd input the full division polynomial is
        just the usual division polynomial, but not for even input:
            sage: E.division_polynomial(2)
            4*x^3 + 8*x + 12
            sage: E.full_division_polynomial(3)
            3*x^4 + 12*x^2 + 36*x - 4
            sage: E.division_polynomial(3)
            3*x^4 + 12*x^2 + 36*x - 4
            sage: E.full_division_polynomial(4)
            4*y*x^6 + 40*y*x^4 + 240*y*x^3 - 80*y*x^2 - 96*y*x - 320*y
            sage: E.full_division_polynomial(5)
            5*x^12 + 124*x^10 + 1140*x^9 - 420*x^8 + 1440*x^7 - 4560*x^6 - 8352*x^5 - 36560*x^4 - 45120*x^3 - 10240*x^2 - 39360*x - 22976

        TESTS:
        We test that the full division polynomial as computed using
        the recurrence agrees with the norml division polynomial for
        a certain curve and all odd $n$ up to $23$:

            sage: E = EllipticCurve([23,-105])
            sage: for n in [1,3,..,23]:
            ...       assert E.full_division_polynomial(n, use_divpoly=False) == E.division_polynomial(n)
        """
        # Coerce the input m to be an integer
        m = rings.Integer(m)

        # Check whether the corresponding poly is cached already
        try:
            return self.__divpoly2[m]
        except AttributeError:
            self.__divpoly2 = {}
        except KeyError:
            pass

        # Get the invariants of the curve and make sure that the curve is
        # in short Weierstrass form.
        a0,a1,a2,A,B = self.a_invariants()
        if a0 or a1 or a2:
            raise NotImplementedError, "Full division polynomial only implemented for elliptic curves in short Weierstrass form"

        # Define the polynomial ring that will contain the full
        # division polynomial.  The one subtle point here is the
        # term order use of the variable y first.  This is done
        # so the natural reduction in the quotient ring mod y^2 - x^3 - A*x - B
        # replaces all y^2's by polys in x.  WARNING: Be careful
        # that the gens are in the order y, x!
        R, (y,x) = PolynomialRing(self.base_ring(), 2, 'y,x', order='revlex').objgens()

        # In case m is odd we can just use the usual division
        # polynomial code.  Note however that we are careful to coerce
        # into the multivariate polynomial ring, since the usual div
        # poly code returns a polynomial in one variable.  If we did
        # that here it would lead to all kinds of problems later.
        if use_divpoly and m % 2 == 1:
            return R(self.division_polynomial(m))

        # Do each case of the recurrence, exactly as in
        # Silverman's exercise in III.3.7 on page 105.
        if m <= 1:
            f = R(1)
            self.__divpoly2[m] = f
            return f
        elif m == 2:
            f = 2*y
            self.__divpoly2[m] = f
            return f
        elif m == 3:
            f = 3*x**4 + 6*A*x**2 + 12*B*x - A**2
            self.__divpoly2[m] = f
            return f
        elif m == 4:
            f = 4*y*(x**6 + 5*A*x**4 + 20*B*x**3 - 5*A**2*x**2 -
                        4*A*B*x - 8*B**2 - A**3)
            self.__divpoly2[m] = f
            return f

        # Finally we do the general part of the recurrence which
        # divides into even and odd cases.
        # We define psi to just be this full_division_polynomial function
        # evaluated at a given integer k.  This makes the code below
        # more readable.
        psi = lambda k: self.full_division_polynomial(k,use_divpoly=use_divpoly)
        if m % 2 == 1:
            n = m//2
            ans = psi(n+2) * psi(n)**3 - psi(n-1) * psi(n+1)**3
        elif m % 2 == 0:
            n = m//2
            ans = (psi(n)*(psi(n+2)*psi(n-1)**2 - psi(n-2)*psi(n+1)**2))/(2*y)

        # Create the affine quotient ring so that we can replace all y^2
        # terms by polys in x.
        Q = R.quotient(y**2 - x**3 - A*x - B)

        # Do the actual reduction and lift back to R.
        f = Q(ans).lift()

        # Cache the result and return it.
        self.__divpoly2[m] = f
        return f

    def multiplication_by_m(self, m, x_only=False):
        """
        Return the multiplication-by-m map from self to self as a rational
        function.

        INPUT:
            self -- an elliptic curve in short Weierstrass form
            m -- a positive integer
            x_only -- bool (default: False) if True, return only the x
                      coordinate of the map.

        OUTPUT:
            2-tuple -- (f(x), g(x,y)) where f and g are rational functions
                       with the degree of y in g(x,y) at most 1.

        NOTE: The result is not cached.

        EXAMPLES:
        We create an elliptic curve.
            sage: E = EllipticCurve([-1,3])

        We verify that multiplication by 1 is just the identity:
            sage: E.multiplication_by_m(1)
            (x, y)

        Multiplication by 2 is more complicated.
            sage: f = E.multiplication_by_m(2)
            sage: f
            ((x^4 + 2*x^2 - 24*x + 1)/(4*x^3 - 4*x + 12), (2*x^6 - 10*x^4 + 120*x^3 - 10*x^2 + 24*x - 142)/(16*y*x^3 - 16*y*x + 48*y))

        Grab only the x-coordinate (less work):
            sage: E.multiplication_by_m(2, x_only=True)
            (x^4 + 2*x^2 - 24*x + 1)/(4*x^3 - 4*x + 12)

        We check that it works on a point:
            sage: P = E([2,3])
            sage: f[0].subs(x=2,y=3)
            -23/36
            sage: f[1].subs(x=2,y=3)
            397/216
            sage: 2*P
            (-23/36 : 397/216 : 1)

        We do the same but with multiplication by 3:
            sage: f = E.multiplication_by_m(3)
            sage: f[0].subs(x=2,y=3)
            -10534/9025
            sage: f[1].subs(x=2,y=3)
            -1376361/857375
            sage: 3*P
            (-10534/9025 : -1376361/857375 : 1)

        And the same with multiplication by 4:
            sage: f = E.multiplication_by_m(4)
            sage: f[0].subs(x=2,y=3)
            29084737/22695696
            sage: f[1].subs(x=2,y=3)
            -211407941663/108122295744
            sage: 4*P
            (29084737/22695696 : -211407941663/108122295744 : 1)

        TESTS:
        Verify for this fairly random looking curve and point that
        multiplication by m returns the right result for the first 10
        integers.
            sage: E = EllipticCurve([23,-105])
            sage: P = E([129/4, 1479/8])
            sage: for n in [1..10]:
            ...       f = E.multiplication_by_m(n)
            ...       Q = n*P
            ...       assert f[0].subs(x=P[0],y=P[1]) == Q[0] and f[1].subs(x=P[0],y=P[1]) == Q[1]
        """
        # Define a function psi that returns the full bivariate division polynomial
        # psi(k) for this elliptic curve.   We do this for simplicity of notation
        # below.
        psi = lambda k: self.full_division_polynomial(k)
        psi_m = psi(m)

        # Grab the multivariate polynomial ring that contains the full division
        # polynomial.  NOTE the generators order y,x rather than x,y!
        R, (y,x) = psi_m.parent().objgens()

        # Special case of multiplication by 1 is easy.
        if m == 1:
            return (x, y)

        # Grab curve invariants and make sure the curve is in short Weierstrass form.
        # (It would be desirable to extend this function to work with
        # arbitrary models -- not just short ones.)
        a0,a1,a2,A,B = self.a_invariants()
        if a0 or a1 or a2:
            raise NotImplementedError, "multiplication_by_m only implemented for elliptic curves in short Weierstrass form"

        # Form the affine coordinate ring, which we'll use only for getting
        # rid of terms involving y^2.
        Q = R.quotient(y**2 - x**3 - A*x - B)
        def normalize(f):
            return Q(f.numerator()).lift() / Q(f.denominator()).lift()

        # Write down the x coordinate using the formula in
        # Silverman AEC Ex III.3.7, page 105.
        phi_m = x*psi(m)**2 - psi(m+1)*psi(m-1)
        x_coord = normalize(phi_m / psi_m**2)
        if x_only:
            # Return it if the optional parameter x_only is set.
            return x_coord

        if m == 2:
            # The formula from Silverman III.3.7 given in the other
            # case of the else is *wrong*.  I guess he made a mistake
            # and ignored a special case.  In any case, the following
            # formula from my elementary number theory book (which
            # is from Lenstra's ECM paper, actually), is right.
            lam = (3*x**2 + A)/(2*y)
            y_coord = normalize(-lam*x_coord - (y - lam*x))
        else:
            # Silverman's formula, which works for m > 2 for
            # the y coordinate.
            omega_m = (psi(m+2)*psi(m-1)**2 - psi(m-2)*psi(m+1)**2)/(4*y)
            y_coord = normalize(omega_m / psi_m**3)
        return x_coord, y_coord

    def isomorphism_to(self, other):
        """
        Given another weierstrass model \code{other} of self, return a morphism
        from self to \code{other}.

        If the curves in question are not isomorphic, raise a ValueError

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: F = E.short_weierstrass_model()
            sage: w = E.isomorphism_to(F); w
            Generic morphism:
            From: Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            To:   Abelian group of points on Elliptic Curve defined by y^2  = x^3 - 16*x + 16 over Rational Field
            Via:  (u,r,s,t) = (1/2, 0, 0, -1/2)
            sage: P = E(0,-1,1)
            sage: w(P)
            (0 : -4 : 1)
            sage: w(5*P)
            (1 : 1 : 1)
            sage: 5*w(P)
            (1 : 1 : 1)
            sage: 120*w(P) == w(120*P)
            True

          We can also handle injections to different base rings:
              sage: K.<a> = NumberField(x^3-7)
              sage: E.isomorphism_to(E.change_ring(K))
              Generic morphism:
                From: Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
                To:   Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + (-1)*x over Number Field in a with defining polynomial x^3 - 7
                Via:  (u,r,s,t) = (1, 0, 0, 0)
        """
        return wm.WeierstrassIsomorphism(self, None, other)

    def automorphisms(self, field=None):
        """
        Return the set of isomorphisms from self to itself (as a list).

        EXAMPLES:
        sage: E = EllipticCurve(QQ(0)) # a curve with j=0 over QQ
        sage: E.automorphisms();
        [Generic endomorphism of Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
        Via:  (u,r,s,t) = (-1, 0, 0, 0), Generic endomorphism of Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
        Via:  (u,r,s,t) = (1, 0, 0, 0)]

        We can also find automorphisms defined over extension fields:
        sage: K.<a> = NumberField(x^2+3) # adjoin roots of unity
        sage: E.automorphisms(K)
        [Generic endomorphism of Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Number Field in a with defining polynomial x^2 + 3
        Via:  (u,r,s,t) = (1, 0, 0, 0),
        ...
        Generic endomorphism of Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Number Field in a with defining polynomial x^2 + 3
        Via:  (u,r,s,t) = (-1/2*a - 1/2, 0, 0, 0)]

        sage: [ len(EllipticCurve(GF(q,'a')(0)).automorphisms()) for q in [2,4,3,9,5,25,7,49]]
        [2, 24, 2, 12, 2, 6, 6, 6]
        """
        if field==None:
            return [wm.WeierstrassIsomorphism(self, urst, self)
                    for urst in wm.isomorphisms(self,self)]
        E=self.change_ring(field)
        return [wm.WeierstrassIsomorphism(E, urst, E)
                for urst in wm.isomorphisms(E,E)]

    def isomorphisms(self, other, field=None):
        """
        Return the set of isomorphisms from self to other (as a list).

        EXAMPLES:
        sage: E = EllipticCurve(QQ(0)) # a curve with j=0 over QQ
        sage: F = EllipticCurve('36a1') # should be the same one
        sage: E.isomorphisms(F);
        [Generic morphism:
        From: Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
        To:   Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
        Via:  (u,r,s,t) = (-1, 0, 0, 0), Generic morphism:
        From: Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
        To:   Abelian group of points on Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
        Via:  (u,r,s,t) = (1, 0, 0, 0)]

        We can also find istomorphisms defined over extension fields:
        sage: E=EllipticCurve(GF(7),[0,0,0,1,1])
        sage: F=EllipticCurve(GF(7),[0,0,0,1,-1])
        sage: E.isomorphisms(F)
        []
        sage: E.isomorphisms(F,GF(49,'a'))
        [Generic morphism:
        From: Abelian group of points on Elliptic Curve defined by y^2  = x^3 + x +1 over Finite Field in a of size 7^2
        To:   Abelian group of points on Elliptic Curve defined by y^2  = x^3 + x + 6 over Finite Field in a of size 7^2
        Via:  (u,r,s,t) = (a + 3, 0, 0, 0), Generic morphism:
        From: Abelian group of points on Elliptic Curve defined by y^2  = x^3 + x +1 over Finite Field in a of size 7^2
        To:   Abelian group of points on Elliptic Curve defined by y^2  = x^3 + x + 6 over Finite Field in a of size 7^2
        Via:  (u,r,s,t) = (6*a + 4, 0, 0, 0)]
        """
        if field==None:
            return [wm.WeierstrassIsomorphism(self, urst, other)
                    for urst in wm.isomorphisms(self,other)]
        E=self.change_ring(field)
        F=other.change_ring(field)
        return [wm.WeierstrassIsomorphism(E, urst, F)
                for urst in wm.isomorphisms(E,F)]

    def is_isomorphic(self, other, field=None):
        """
        Returns whether or not self is isomorphic to other, i.e. they
        define the same curve over the same basering.

        If field!=None then both curves must base_extend-able to it
        and the isomorphism is then checked over that field

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: F = E.change_weierstrass_model([2,3,4,5]); F
            Elliptic Curve defined by y^2 + 4*x*y + 11/8*y = x^3 - 3/2*x^2 - 13/16*x over Rational Field
            sage: E.is_isomorphic(F)
            True
            sage: E.is_isomorphic(F.change_ring(CC))
            False
        """
        if not is_EllipticCurve(other):
            return False
        if field==None:
            if self.base_ring() != other.base_ring():
                return False
            elif self.j_invariant() != other.j_invariant(): ## easy check
                return False
            else:
                return wm.isomorphisms(self,other,True) != None
        else:
            E=self.base_extend(field)
            F=other.base_extend(field)
            if E.j_invariant() != F.j_invariant(): ## easy check
                return False
            else:
                return wm.isomorphisms(E,other,F) != None

    def change_weierstrass_model(self, *urst):
        r"""
        Return a new Weierstrass model of self under the transformation (on points)
            $$ (x,y) \mapsto (x',y') = (u^2*x+r , u^3*y + s*u^2*x' + t) $$

        EXAMPLES:
            sage: E = EllipticCurve('15a')
            sage: F1 = E.change_weierstrass_model([1/2,0,0,0]); F1
            Elliptic Curve defined by y^2 + 2*x*y + 8*y = x^3 + 4*x^2 - 160*x - 640 over Rational Field
            sage: F2 = E.change_weierstrass_model([7,2,1/3,5]); F2
            Elliptic Curve defined by y^2 + 5/21*x*y + 13/343*y = x^3 + 59/441*x^2 - 10/7203*x - 58/117649 over Rational Field
            sage: F1.is_isomorphic(F2)
            True
        """
        if isinstance(urst[0], (tuple, list)):
            urst = urst[0]
        return constructor.EllipticCurve((wm.baseWI(*urst))(self.ainvs()))

    def short_weierstrass_model(self, complete_cube=True):
        """
        Return a short Weierstrass model for self.

        INPUT:
            complete_cube -- bool (default: True); for meaning, see below.
        OUTPUT:
            an elliptic curve

        If complete_cube=True:
        Return a model of the form $y^2 = x^3 + a*x + b$ for this curve.
        The characteristic must not be 2 or 3.
        a,b = -27*c4, -54*c6

        If complete_cube=False:
        Return a model of the form $y^2 = x^3 + ax^2 + bx + c$ for this curve.
        The characteristic must not be 2.
        a,b,c = b2, 8*b4, 16*b6

        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: print E
            Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field
            sage: F = E.short_weierstrass_model()
            sage: print F
            Elliptic Curve defined by y^2  = x^3 + 4941*x + 185166 over Rational Field
            sage: E.is_isomorphic(F)
            True
            sage: F = E.short_weierstrass_model(complete_cube=False)
            sage: print F
            Elliptic Curve defined by y^2  = x^3 + 9*x^2 + 88*x + 464 over Rational Field
            sage: print E.is_isomorphic(F)
            True

            sage: E = EllipticCurve(GF(3),[1,2,3,4,5])
            sage: E.short_weierstrass_model(complete_cube=False)
            Elliptic Curve defined by y^2  = x^3 + x + 2 over Finite Field of size 3
            sage: E.short_weierstrass_model()
            Traceback (most recent call last):
            ...
            ValueError: short_weierstrass_model(): no short model for Elliptic Curve defined by y^2 + x*y  = x^3 + 2*x^2 + x + 2 over Finite Field of size 3 (characteristic is 3)

        """
        import constructor
        K = self.base_ring()
        if K.characteristic() == 2:
            raise ValueError, "short_weierstrass_model(): no short model for %s (characteristic is %s)"%(self,K.characteristic())
        if K.characteristic() == 3 and complete_cube:
            raise ValueError, "short_weierstrass_model(): no short model for %s (characteristic is %s)"%(self,K.characteristic())

        a1,a2,a3,_,_ = self.a_invariants()
        if complete_cube:
            if a1==0 and a2==0 and a3==0:
                return self
            else:
                b2,b4,b6,_ = self.b_invariants()
                if b2==0:
                    return constructor.EllipticCurve([0,0,0,8*b4,16*b6])
                else:
                    c4, c6 = self.c_invariants()
                    return constructor.EllipticCurve([0,0,0,-27*c4, -54*c6])
        else:
            if a1==0 and a3==0:
                return self
            else:
                b2,b4,b6,_ = self.b_invariants()
                return constructor.EllipticCurve([0,b2,0,8*b4,16*b6])


    ##############################################################################
    # Plotting
    ##############################################################################

    def plot(self, xmin=None, xmax=None, **args):
        """
        Draw a graph of this elliptic curve.

        INPUT:
            xmin, xmax -- points will be computed at least within this
                          rings, but possibly farther.  These may be
                          left off.
            **args -- all other options are passed to the line graphing
                      primitive.

        EXAMPLES:
            sage: E = EllipticCurve([0,-1])
            sage: plot(E, rgbcolor=hue(0.7))
        """
        RR = rings.RealField()
        K = self.base_ring()
        try:
            RR._coerce_(K(1))
        except TypeError:
            raise NotImplementedError, "Plotting of curves over %s not implemented yet"%K
        a1, a2, a3, a4, a6 = self.ainvs()
        R = rings.PolynomialRing(rings.RealField(), 'x')
        x = R.gen()
        d = 4*x**3 + (a1**2 + 4*a2)*x**2 + (2*a3*a1 + 4*a4)*x + (a3**2 + 4*a6)
        def f1(z):
            """
            Internal function for plotting first branch of the curve
            """
            return (-(a1*z + a3) + sqrt(abs(d(z))))/2
        def f2(z):
            """
            Internal function for plotting second branch of the curve
            """
            return (-(a1*z + a3) - sqrt(abs(d(z))))/2
        r = d.roots(multiplicities=False)
        r.sort()
        if xmax is None:
            xmax = r[-1] + 2
        xmax = max(xmax, r[-1]+2)
        if xmin is None:
            xmin = r[0]  - 2
        xmin = min(xmin, r[0]-2)
        if len(r) == 1:
            # one real root; 1 component
            I = [(r[0],xmax)]
        else:
            # three real roots; 2 components
            I = [(r[0],r[1]), (r[2],xmax)]
        I = [(max(a,xmin),min(b,xmax)) for a,b in I]

        g = plot.Graphics()
        try:
            plot_points = int(args['plot_points'])
            del args['plot_points']
        except KeyError:
            plot_points = 100

        for j in range(len(I)):
            a,b = I[j]
            delta = (b-a)/float(plot_points)
            v = []
            w = []
            for i in range(plot_points):
                x = a + delta*i
                v.append((x, f1(x)))
                w.append((x, f2(x)))
            v.append((b,f1(b)))
            w.append((b,f2(b)))
            if len(I) == 2 and j == 0:  # two components -- the oh.
                g += plot.line(v + list(reversed(w)) + [v[0]], **args)
            else:
                g += plot.line(list(reversed(v)) + w, **args)
        return g

    def formal_group(self):
        r"""
        The formal group associated to this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve("37a")
            sage: E.formal_group()
            Formal Group associated to the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        try:
            return self.__formal_group
        except AttributeError:
            self.__formal_group = formal_group.EllipticCurveFormalGroup(self)
            return self.__formal_group

    formal = formal_group


    def hyperelliptic_polynomials(self):
        r""" Returns a pair of polynomials g(x), h(x) such that this elliptic
        curve can be defined by the standard hyperelliptic equation
        $$y^2 + h(x)y = g(x)$$.

        EXAMPLES:
            sage: R.<a1,a2,a3,a4,a6>=QQ[]
            sage: E=EllipticCurve([a1,a2,a3,a4,a6])
            sage: E.hyperelliptic_polynomials()
            (x^3 + a2*x^2 + a4*x + a6, a1*x + a3)
        """
        K = self.base_ring()
        R = PolynomialRing(K, 'x')
        x = R.gen(0)
        a1, a2, a3, a4, a6 = self.ainvs()
        return R([a6, a4, a2, 1]), R([a3, a1])


def Hasse_bounds(q, genus=1):
    """
    Return the Hasse bounds (lb,ub) for the cardinality of a curve of
    genus g (default 1) defined over GF(q)

    EXAMPLES:
       sage: Hasse_bounds(2)
       (1, 5)
       sage: Hasse_bounds(next_prime(10^30))
       (999999999999998000000000000058, 1000000000000002000000000000058)
    """
    rq = 2*genus*q.isqrt()
    return (q+1-rq,q+1+rq)

