"""
Elliptic curves over a general ring
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

from sage.rings.all import MPolynomialRing

import sage.plot.all as plot

import sage.rings.arith as arith
import sage.rings.all as rings
import sage.rings.number_field as number_field
from sage.rings.all import is_Infinite
import sage.misc.misc as misc
import sage.misc.latex as latex
import sage.modular.modform as modform
import sage.functions.transcendental as transcendental

# Schemes
import sage.schemes.generic.projective_space as projective_space
import sage.schemes.generic.homset as homset

import ell_point
import constructor
import formal_group

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
        if K.is_field():
            self._point_class = ell_point.EllipticCurvePoint_field
        else:
            self._point_class = ell_point.EllipticCurvePoint

    def _defining_params_(self):
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
        a = self.ainvs()
        #return "y^2 + %s*x*y + %s*y = x^3 + %s*x^2 + %s*x + %s"%\
        #       (a[0], a[2], a[1], a[3], a[4])
        a = [z._coeff_repr() for z in a]
        s = "Elliptic Curve defined by "
        s += "y^2 "
        if a[0] == "-1":
            s += "- x*y "
        elif a[0] == '1':
            s += "+ x*y "
        elif a[0] != '0':
            s += "+ %s*x*y "%a[0]
        if a[2] == "-1":
            s += " - y"
        elif a[2] == '1':
            s += "+ y"
        elif a[2] != '0':
            s += "+ %s*y"%a[2]
        s += " = x^3 "
        if a[1] == "-1":
            s += "- x^2 "
        elif a[1] == '1':
            s += "+ x^2 "
        elif a[1] != '0':
            s += "+ %s*x^2 "%a[1]
        if a[3] == "-1":
            s += "- x "
        elif a[3] == '1':
            s += "+ x "
        elif a[3] != '0':
            s += "+ %s*x "%a[3]
        if a[4] == '-1':
            s += "-1 "
        elif a[4] == '1':
            s += "+1 "
        elif a[4] != '0':
            s += "+ %s "%a[4]
        s = s.replace("+ -","- ")
        s += "over %s"%self.base_ring()
        return s

    def _latex_(self):
        a = self.ainvs()
        a = [z._latex_coeff_repr() for z in a]
        s = "y^2 "
        if a[0] == '-1':
            s += "- xy "
        elif a[0] == '1':
            s += "+ xy "
        elif a[0] != '0':
            s += "+ %sxy "%a[0]
        if a[2] == '-1':
            s += " - y"
        elif a[2] == '1':
            s += "+ y"
        elif a[2] != '0':
            s += "+ %sy"%a[2]
        s += " = x^3 "
        if a[1] == '-1':
            s += "- x^2 "
        elif a[1] == '1':
            s += "+ x^2 "
        elif a[1] != '0':
            s += "+ %sx^2 "%a[1]
        if a[3] == '-1':
            s += "- x "
        elif a[3] == '1':
            s += "+ x "
        elif a[3] != '0':
            s += "+ %sx "%a[3]
        if a[4] == '-1':
            s += "-1 "
        elif a[4] == '1':
            s += "+1 "
        elif a[4] != '0':
            s += "+ %s "%a[4]
        s = s.replace("+ -","- ")
        return s

    def _pari_init_(self):
        return 'ellinit([%s])'%(','.join([x._pari_init_() for x in self.ainvs()]))

    def _magma_init_(self):
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
            (y^2 + y) == (x^3 - x^2 - 10*x - 20)
            sage: print eqn
                                      2        3    2
                                     y  + y == x  - x  - 10 x - 20

        We verify that the given point is on the curve:
            sage: eqn(x=5,y=5)
            (30) == (30)
            sage: bool(eqn(x=5,y=5))
            True

        We create a single expression:
            sage: F = eqn.lhs() - eqn.rhs(); print F
                                      2        3    2
                                     y  + y - x  + x  + 10 x + 20
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
            ((3*pi^2 - 2*pi - 10)^2/(4*pi^3 - 4*pi^2 - 40*pi - 79) - 2*pi + 1 : (sqrt(4*pi^3 - 4*pi^2 - 40*pi - 79) + 1)/2 - ((3*pi^2 - 2*pi - 10)*((-(3*pi^2 - 2*pi - 10)^2)/(4*pi^3 - 4*pi^2 - 40*pi - 79) + 3*pi - 1)/(sqrt(4*pi^3 - 4*pi^2 - 40*pi - 79))) - 1 : 1)
        """
        a = [SR(x) for x in self.a_invariants()]
        x, y = SR.var('x, y')
        return y**2 + a[0]*x*y + a[2]*y == x**3 + a[1]*x**2 + a[3]*x + a[4]

    def __cmp__(self, other):
        if not isinstance(other, EllipticCurve_generic):
            return -1
        return misc.generic_cmp(self.ainvs(), other.ainvs())

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
            TypeError: coordinates [0, 0, 1] do not define a point on Elliptic Curve
            defined by y^2  = x^3 +1 over Finite Field of size 7

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
        """
        if len(args) == 1 and args[0] == 0:
            R = self.base_ring()
            return self.point([R(0),R(1),R(0)], check=False)
        return plane_curve.ProjectiveCurve_generic.__call__(self, *args, **kwds)

    def lift_x(self, x):
        """
        Given the x-coordinate of a point on the curve, use the defining
        polynomial to find all affine points on this curve with the
        given x-coordinate.

        EXAMPLES:
            sage: E = EllipticCurve('37a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: E.lift_x(1)
            [(1 : 0 : 1), (1 : -1 : 1)]
            sage: E.lift_x(2)
            [(2 : 2 : 1), (2 : -3 : 1)]
            sage: E.lift_x(1/4)
            [(1/4 : -3/8 : 1), (1/4 : -5/8 : 1)]

        There are no rational points with x-cordinate 3.
            sage: E.lift_x(3)
            []

        However, there are two such points in $E(\R)$:
            sage: E.change_ring(RR).lift_x(3)
            [(3.00000000000000 : 4.42442890089805 : 1.00000000000000), (3.00000000000000 : -5.42442890089805 : 1.00000000000000)]

        And of course it always works in $E(\C)$:
            sage: E.change_ring(RR).lift_x(.5)
            []
            sage: E.change_ring(CC).lift_x(.5)
            [(0.500000000000000 : -0.500000000000000 + 0.353553390593274*I : 1.00000000000000), (0.500000000000000 : -0.500000000000000 - 0.353553390593274*I : 1.00000000000000)]

        We can perform these operations over finite fields too:
            sage: E = E.change_ring(GF(17)); E
            Elliptic Curve defined by y^2 + y = x^3 + 16*x over Finite Field of size 17
            sage: E.lift_x(7)
            [(7 : 11 : 1), (7 : 5 : 1)]
            sage: E.lift_x(3)
            []

        Note that there is only one lift with x-coordinate 10 in $E(\F_{17})$.
            sage: E.lift_x(10)
            [(10 : 8 : 1)]

        We can lift over more exotic rings too.
            sage: E = EllipticCurve('37a');
            sage: E.lift_x(pAdicField(17, 5)(6))
            [(6 + O(17^5) : 2 + 16*17 + 16*17^2 + 16*17^3 + 16*17^4 + O(17^5) : 1 + O(17^5)), (6 + O(17^5) : 14 + O(17^5) : 1 + O(17^5))]
            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: E.lift_x(1+t)
            [(1 + t : 2*t - t^2 + 5*t^3 - 21*t^4 + O(t^5) : 1), (1 + t : -1 - 2*t + t^2 - 5*t^3 + 21*t^4 + O(t^5) : 1)]

        TEST:
            sage: E = EllipticCurve('37a').weierstrass_model().change_ring(GF(17))
            sage: E.lift_x(3)
            []
            sage: E.lift_x(7)
            [(7 : 3 : 1), (7 : 14 : 1)]
        """
        a1, a2, a3, a4, a6 = self.ainvs()
        f = ((x + a2) * x + a4) * x + a6
        x += self.base_ring()(0)
        one = x.parent()(1)
        if a1.is_zero() and a3.is_zero():
            if f.is_zero():
                return [self.point([x,x.parent()(0), one], check=False)]
            if not f.is_square():
                return []
            else:
                y = f.sqrt(extend=False)
                return [self.point([x,  y, one], check=False),
                        self.point([x, -y, one], check=False)]
        else:
            b = (a1*x + a3)
            D = b*b + 4*f
            if D.is_zero():
                return [self.point([x, -b/2, one], check=False)]
            if not D.is_square():
                return []
            else:
                sqrtD = D.sqrt(extend=False)
                return [self.point([x, (-b+sqrtD)/2, one], check=False),
                        self.point([x, (-b-sqrtD)/2, one], check=False)]

    def _homset_class(self, *args, **kwds):
        return homset.SchemeHomsetModule_abelian_variety_coordinates_field(*args, **kwds)

    def __getitem__(self, n):
        """
        """
        raise NotImplementedError, "not implemented."

    def __is_over_RationalField(self):
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
        try:
            return self.__c_invariants[0]
        except AttributeError:
            pass
        return self.c_invariants()[0]

    def c6(self):
        try:
            return self.__c_invariants[1]
        except AttributeError:
            pass
        return self.c_invariants()[1]


    def base_extend(self, R):
        return constructor.EllipticCurve(R, [R(a) for a in self.a_invariants()])

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
        raise NotImplementedError

    def gen(self, i):
        return self.gens()[i]

    def quadratic_twist(self, D):
        """
        Return the quadratic twist of this curve by D.

        EXAMPLES:
            sage: E = EllipticCurve([GF(1103)(1), 0, 0, 107, 340]); E
            Elliptic Curve defined by y^2 + x*y  = x^3 + 107*x + 340 over Finite Field of size 1103
            sage: E.quadratic_twist(-1)
            Elliptic Curve defined by y^2  = x^3 + 770*x + 437 over Finite Field of size 1103
        """
        a = self.ainvs()
        ap = [0, 0, 0]
        ap.append(-27*D**2*a[0]**4 - 216*D**2*a[0]**2*a[1] + 648*D**2*a[0]*a[2] - 432*D**2*a[1]**2 + 1296*D**2*a[3])
        ap.append(54*D**3*a[0]**6 + 648*D**3*a[0]**4*a[1]
                  - 1944*D**3*a[0]**3*a[2] + 2592*D**3*a[0]**2*a[1]**2
                  - 3888*D**3*a[0]**2*a[3] - 7776*D**3*a[0]*a[1]*a[2]
                  + 3456*D**3*a[1]**3 - 15552*D**3*a[1]*a[3]
                  + 11664*D**3*a[2]**2 + 46656*D**3*a[4])
        import constructor
        return constructor.EllipticCurve(self.base_ring(), ap)

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

        SEE ALSO:
            pseudo_torsion_polynomial()

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

    def weierstrass_model(self):
        """
        Return a model of the form $y^2 = x^3 + a*x + b$ for this curve.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: print E
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: F = E.weierstrass_model()
            sage: print F
            Elliptic Curve defined by y^2  = x^3 - x + 1/4 over Rational Field
            sage: print F.minimal_model() == E.minimal_model()
            True

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: F = E.weierstrass_model()
            sage: print F
            Elliptic Curve defined by y^2  = x^3 + 61/16*x + 127/32 over Rational Field
            sage: print F.minimal_model() == E.minimal_model()
            True
        """
        import constructor
        c4, c6 = self.c_invariants()
        return constructor.EllipticCurve([-c4/(2**4*3), -c6/(2**5*3**3)])



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
            Graphics object consisting of 1 graphics primitive
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
            return (-(a1*z + a3) + sqrt(abs(d(z))))/2
        def f2(z):
            return (-(a1*z + a3) - sqrt(abs(d(z))))/2
        r = [t.real() for t in d.roots() if t.imag() == 0]
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
