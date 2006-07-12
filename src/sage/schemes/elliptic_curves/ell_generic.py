"""
Elliptic curves over a general ring
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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
from sage.rings.all import is_Infinity
import sage.misc.misc as misc
import sage.misc.latex as latex
import sage.modular.modform as modform
import sage.functions.transcendental as transcendental

# Schemes
import sage.schemes.generic.projective_space as projective_space

import ell_point
import constructor


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

    def __call__(self, *args):
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
        return plane_curve.ProjectiveCurve_generic.__call__(self, *args)

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
            sage: E = EllipticCurve(GF(49), [3,5])
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

    def quadratic_twist(self, D):
        """
        Return the quadratic twist of this curve by D.
        """
        a = self.ainvs()
        import constructor
        return constructor.EllipticCurve([0, 0, 0, \
                              -27*D**2*a[0]**4 - 216*D**2*a[0]**2*a[1] + 648*D**2*a[0]*a[2] - 432*D**2*a[1]**2 + 1296*D**2*a[3], \
                              54*D**3*a[0]**6 + 648*D**3*a[0]**4*a[1] - 1944*D**3*a[0]**3*a[2] + 2592*D**3*a[0]**2*a[1]**2 - 3888*D**3*a[0]**2*a[3] - 7776*D**3*a[0]*a[1]*a[2] + 3456*D**3*a[1]**3 - 15552*D**3*a[1]*a[3] + 11664*D**3*a[2]**2 + 46656*D**3*a[4]])

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


    def torsion_polynomial(self, n, i=0):
        """
        Returns the n-th torsion polynomial (a.k.a., division polynomial).

        INPUT:
            n -- non-negative integer
            i -- integer, either 0 (default) or 1.

        OUTPUT:
            Polynomial -- n-th torsion polynomial, which is a polynomial over
                          the base field of the elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: E.division_polynomial(1)
            1
            sage: E.division_polynomial(2)
            4*x^3 - 4*x + 1
            sage: E.division_polynomial(3)
            3*x^4 - 6*x^2 + 3*x - 1

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

        R = rings.PolynomialRing(E.base_ring())
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

    def _plot_(self, xmin=None, xmax=None, **args):
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
        R = rings.PolynomialRing(rings.RealField())
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


    ##############################################################################
    # Formal Groups
    ##############################################################################

    def formal_w(self, prec=20):
        """
        The formal group power series w.

        INPUT:
            prec -- integer
        OUTPUT:
            a power series

        DETAILS:
            Return the formal power series
            $$
                   w(t) = t^3 + ...
            $$
            to precision $O(t^prec)$ of Proposition IV.1.1 of [Silverman
            AEC1].  This is the formal expansion of $w = -1/y$ about the
            formal parameter $t = -x/y$ at $\\infty$.

            We compute w using the recursive procedure (4.1) on page 18
            of Bluher's ``A leisurely introduction to formal groups and
            elliptic curves'', which I downloaded from
                http://www.math.uiuc.edu/Algebraic-Number-Theory/0076/

        EXAMPLES:
            sage: e = EllipticCurve([0, 0, 1, -1, 0])
            sage: e.formal_w(10)
            t^3 + t^6 - t^7 + 2*t^9 + O(t^10)
        """
        k = self.base_ring()
        try:
            w = self.__formal_w
            pr = len(w)
        except AttributeError:
            pr = 4
            w = [k(0), k(0), k(0), k(1)]
            self.__formal_w = w

        R = rings.PowerSeriesRing(k, "t")
        if prec <= pr:
            return R(w, prec)
        a1, a2, a3, a4, a6 = self.ainvs()
        t0 = misc.cputime()

        for n in range(pr,prec):
            w.append(a1*w[n-1] + \
                     a2*w[n-2] + \
                     a3*sum([w[i]*w[n-i] for i in range(3,n-2)]) + \
                     a4*sum([w[i]*w[n-1-i] for i in range(3,n-3)])  + \
                     a6*sum([w[i]*w[j]*w[n-i-j] for i in range(3,n-2) \
                             for j in range(3,n-i-2)]))

        return R(w, prec)


##         for n in range(4,prec):
##             if n<pr:
##                 w[n] = self.__formal_w[1][n]   # already know it
##             else:
##                 w[n] = a1*w[n-1] + \
##                        a2*w[n-2] + \
##                        a3*sum([w[i]*w[n-i] for i in range(3,n-2)]) + \
##                        a4*sum([w[i]*w[n-1-i] for i in range(3,n-3)]) + \
##                        a6*sum([w[i]*w[j]*w[n-i-j] for i in range(3,n-2) \
##                                for j in range(3,n-2) if n-i-j >= 3])
##         self.__formal_w = (prec, w)
##         return w

    def formal_x(self, prec=20):
        """
        Return the formal power series x(t) in terms of the local
        parameter t = -x/y at infinity.
        """
        prec = max(prec,0)
        try:
            pr, x = self.__formal_x
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = x.parent().gen()
            return x + O(t**prec)
        w = self.formal_w(prec+5)
        t = w.parent().gen()
        x = t/w
        self.__formal_x = (prec, x)
        return x

    def formal_y(self, prec=20):
        """
        Return the formal power series y(t) in terms of the local
        parameter t = -x/y at infinity.
        """
        prec = max(prec,0)
        try:
            pr, y = self.__formal_y
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = y.parent().gen()
            return y + O(t**prec)
        w = self.formal_w(prec+6)
        t = w.parent().gen()
        y = -1/w
        self.__formal_y = (prec, y)
        return y

    def formal_log(self, prec=20):
        a = self.ainvs()
        x = self.formal_x(prec)
        y = self.formal_y(prec)
        xprime = x.derivative()
        g = xprime / (2*y + a[0]*x + a[2])
        return g.integral().power_series().add_bigoh(prec)

    def formal_inverse(self, prec=20):
        prec = max(prec,0)
        try:
            pr, inv = self.__formal_inverse
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = inv.parent().gen()
            return inv + O(t**prec)
        x = self.formal_x(prec)
        y = self.formal_y(prec)
        a1, _, a3, _, _ = self.ainvs()
        inv = x / ( y + a1*x + a3)          # page 114 of Silverman, AEC I
        self.__formal_inverse = (prec, inv)
        return inv

    def formal_n_isogeny(self, n, prec=20):
        prec = max(prec,0)
        try:
            pr, phi = self.__formal_n_isogeny[n]
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = phi.parent().gen()
            return phi + O(t**prec)
        try:
            _ = self.__formal_n_isogeny
        except AttributeError:
            self.__formal_n_isogeny = {}
        misc.todo()

    def formal_group(self, prec=10):
        prec = max(prec,0)
        try:
            pr, F = self.__formal_group
        except AttributeError:
            pr = -1
        if prec <= pr:
            return F
        R1 = rings.PowerSeriesRing(self.base_ring(),"t1")
        R2 = rings.PowerSeriesRing(R1,"t2")
        t1 = R1.gen().add_bigoh(prec)
        t2 = R2.gen().add_bigoh(prec)
        w = self.formal_w(prec)
        def tsum(n):
            return sum([t2**m * t1**(n-m-1) for m in range(n)])
        lam = sum([tsum(n)*w[n] for n in range(3,prec)])
        w1 = R1(w.coeffs()).add_bigoh(prec)
        nu = w1 - lam*t1
        a1, a2, a3, a4, a6 = self.ainvs()
        lam2 = lam*lam
        lam3 = lam2*lam
        t3 = -t1 - t2 - \
             (a1*lam + a3*lam2 + a2*nu + 2*a4*lam*nu + 3*a6*lam2*nu)/  \
             (1 + a2*lam + a4*lam2 + a6*lam3)
        inv = self.formal_inverse(prec)
        F = inv(t3)
        self.__formal_group = (prec, F)
        return F

    def formal_mult(self, n, prec=10):
        R = rings.PowerSeriesRing(self.base_ring(),"t")
        t = R._0
        if n == 1:
            return t +O(t**prec)
        if n == -1:
            return self.formal_inverse(prec)
        if n < 0:
            return self.formal_inverse(prec)(self.formal_mult(-n,prec))
        F = self.formal_group(prec)
        g = F.parent().base_ring().gen()
        for m in range(2,n+1):
            g = F(g)
        return R(g.coeffs())

    def formal_sigma(self, prec=10):
        a1,a2,a3,a4,a6 = self.ainvs()

        k = self.base_ring()
        fl = self.formal_log(prec)
        R = rings.PolynomialRing(k,'c'); c = R.gen()
        F = fl.reversion()

        S = rings.LaurentSeriesRing(R,'z')
        c = S(c)
        z = S.gen()
        F = F(z + O(z**prec))
        wp = self.formal_x()(F)
        e2 = 12*c - a1**2 - 4*a2
        g = (1/z**2 - wp + e2/12).power_series()
        h = g.integral().integral()
        sigma_of_z = z.power_series() * h.exp()

        T = rings.PowerSeriesRing(R,'t')
        fl = fl(T.gen()+O(T.gen()**prec))
        sigma_of_t = sigma_of_z(fl)
        return sigma_of_t


    ##############################################################################
    # End Formal Groups
    ##############################################################################


