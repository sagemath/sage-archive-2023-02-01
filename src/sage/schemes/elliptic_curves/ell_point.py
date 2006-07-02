"""
Points on elliptic curves
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

import sage.plot.all as plot

import ell_generic
import sage.rings.all as rings
import sage.rings.arith as arith
import sage.misc.misc as misc
from sage.schemes.generic.morphism import (SchemeMorphism_projective_coordinates_ring,
                                           SchemeMorphism_abelian_variety_coordinates_field)

oo = rings.infinity       # infinity

class EllipticCurvePoint(SchemeMorphism_projective_coordinates_ring):
    """
    A point on an elliptic curve.
    """
    def __cmp__(self, other):
        if isinstance(other, (int, long, rings.Integer)) and other == 0:
            if self.is_zero():
                return 0
            else:
                return -1
        return SchemePoint_projective_abelian_scheme.__cmp__(self, other)

class EllipticCurvePoint_field(SchemeMorphism_abelian_variety_coordinates_field):
    """
    A point on an elliptic curve over a field.  The point has coordinates
    in the base field.

    EXAMPLES:
        sage: E = EllipticCurve('37a')
        sage: E([0,0])
        (0 : 0 : 1)
        sage: E(0,0)               # brackets are optional
        (0 : 0 : 1)
        sage: E([GF(5)(0), 0])     # entries are coerced
        (0 : 0 : 1)

        sage: E(0.000, 0)
        (0 : 0 : 1)

        sage: E(1,0,0)
        Traceback (most recent call last):
        ...
        TypeError: coordinates [1, 0, 0] do not define a point on
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        sage: E = EllipticCurve([0,0,1,-1,0])
        sage: S = E(QQ); S
        Set of Rational Points of Elliptic Curve defined by y^2 + y = x^3 - x
        over Rational Field

        sage: loads(S.dumps()) == S
        True
        sage: P = E(0,0); P
        (0 : 0 : 1)
        sage: loads(P.dumps()) == P
        True
        sage: T = 100*P
        sage: loads(T.dumps()) == T
        True
    """

    def order(self):
        """
        Return the order of this point on the elliptic curve.
        If the point has infinite order, returns 0.

        EXAMPLE:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.order()
            Infinity

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.order()
            2
        """
        if self.is_zero():
            return rings.Integer(1)
        E = self.curve()
        try:
            n = int(E.pari_curve().ellorder([self[0], self[1]]))
            if n == 0:
                return oo
            return n
        except AttributeError, msg:
            raise NotImplementedError, "%s\nComputation of order of point."%msg

    def curve(self):
        return self.scheme()

    def is_zero(self):
        return self[2] == 0

    def _plot_(self, **args):
        if self.is_zero():
            return plot.text("$\\infty$", (-3,3), **args)

        else:
            return plot.point((self[0], self[1]), **args)

    def _add_(self, right):
        """
        Add self to right.
        """
        # Use Prop 7.1.7 of Cohen "A Course in Computational Algebraic Number Theory"
        if self.is_zero():
            return right
        if right.is_zero():
            return self
        E = self.curve()
        a1, a2, a3, a4, a6 = E.ainvs()
        x1, y1 = self[0], self[1]
        x2, y2 = right[0], right[1]
        if x1 == x2 and y1 == -y2 - a1*x2 - a3:
            return E(0) # point at infinity

        if x1 == x2 and y1 == y2:
            m = (3*x1*x1 + 2*a2*x1 + a4 - a1*y1) / (2*y1 + a1*x1 + a3)
        else:
            m = (y1-y2)/(x1-x2)

        x3 = -x1 - x2 - a2 + m*(m+a1)
        y3 = -y1 - a3 - a1*x3 + m*(x1-x3)
        return E.point([x3, y3, 1], check=False)

    def _sub_(self, right):
        return self + (-right)

    def __neg__(self):
        if self.is_zero():
            return self
        E, x, y = self.curve(), self[0], self[1]
        return E.point([x, -y - E.a1()*x - E.a3(), 1], check=False)

    def height(self):
        """
        The Neron-Tate canonical height of the point.

        EXAMPLES:
            sage: E = EllipticCurve('11a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: P = E([5,5]); P
            (5 : 5 : 1)
            sage: h = P.height()
            sage: RR = h.parent()
            sage: RR.scientific_notation(True)
            sage: h
            -1.4246355279999999e-248
            sage: Q = 5*P
            sage: Q.height()
            0.0000000000000000e-1

            sage: E = EllipticCurve('37a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: P = E([0,0])
            sage: P.height()
            5.1111408239968840e-2
            sage: P.order()
            Infinity
            sage: E.regulator()
            0.051111408239968840
        """
        if self.is_zero():
            return rings.RR(0)
        h = self.curve().pari_curve().ellheight([self[0], self[1]])
        return rings.RR(h)

    def xy(self):
        if self[2] == 1:
            return self[0], self[1]
        else:
            return self[0]/self[2], self[1]/self[2]

    def sigma(self, p):
        k = rings.Qp(p)
        if self.is_zero():
            return k(0)
        sigma = self.curve().sigma(p)
        return sigma(k(-self[0]/self[1]))

##########################################################################################



##     def padic_height(self, p):
##         """Returns the p-adic cyclotomic height of the point.

##         padic_height(p)

##         Input:
##            p: a prime number
##         Output:
##            p-adic cyclotomic height, which is a single p-adic number.

##         Algorithm:
##         We compute this height using the following formula, which
##         is valid for points that are in the intersection of the
##         identity component of the Neron model with the kernel of
##         reduction modulo $p$:
## $$
##      h(P) = 1/2 * sum_{ell!=p} sup(0,-ord_ell(x(P)))  + log_p(sigma_p(-x(P)/y(P)) / e),
## $$
##         where $P=(a/e^2, b/e^3)$ with $\gcd(a,e)=1$, and
##         where the first sum is over primes ell that don't equal $p$.
##         If $P$ isn't in the subgroup mentioned above, let $n$ be a positive
##         integer so that $nP$ is in that subgroup.  Then we return
##         $h(nP)/(n^2)$, which does not depend on the choice of $n$, and
##         is defined using the above formula.

##         """

##         if not arith.is_prime(p):
##             raise ArithmeticError, "p must be a prime number."
##         E = self.curve()
##         if E.conductor() % p == 0:
##             raise ArithmeticError, "p must be a prime of good reduction."

##         ap = E.ap(p)
##         if ap == 0:
##             raise ArithmeticError, "p must a prime of good ordinary reduction."

##         Np = p+1 - E.ap(p)
##         c = 1
##         for q, _ in arith.factor(E.conductor()):
##             c = arith.lcm(c, E.tamagawa_number(q))

##         n = Np*c
##         misc.verbose("Multiplying point by %s."%n)
##         R = n*self

##         x, y = R.xy()
##         sigma = E.padic_sigma(p)
##         Qp = rings.pAdicField(p)
##         print "x = ", x
##         print "-x/y = ", Qp(-x/y)
##         d = x.denominator()
##         e = d.sqrt()
##         if e*e != d:
##             raise ArithmeticError, "The denominator is not a perfect square!?"
##         t = Qp(-x/y)
##         return (sigma(t)/Qp(e)).log()



##         w = Qp(0)
##         for ell, ord in d.factor():
##             if ell != p:
##                 w += ord * Qp(ell).log()
##         # note that there is no 1/2, since we factored the square root.
##         t = Qp(-x/y)
##         print "sigm = ", sigma(t)
##         print "e = ", e
##         hP = w + (sigma(t)/e).log()
##         print "w = ", w
##         print "sigma(t) = ", sigma(t)
##         print "sigma(t)/e = ", sigma(t)/e
##         print "(sigma(t)/e).log() = ", (sigma(t)/e).log()
##         print "sum = ", w + (sigma(t)/e).log()

##         return hP/Qp(n**2)



class EllipticCurvePoint_finite_field(EllipticCurvePoint_field):
    def order(self):
        """
        Return the order of this point on the elliptic curve.
        If the point has infinite order, returns 0.

        EXAMPLE:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.order()
            Infinity

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.order()
            2
        """
        try:
            return self.__order
        except AttributeError:
            pass
        if self.is_zero():
            return rings.Integer(1)
        E = self.curve()
        K = E.base_ring()
        if K.is_prime():
            e = E._gp()
            self.__order = rings.Integer(e.ellzppointorder(list(self.xy())))
        else:
            print "WARNING -- using naive point order finding over finite field!"
            # TODO: This is very very naive!!  -- should use baby-step giant step; maybe in mwrank
            #      note that this is *not* implemented in PARI!
            P = self
            n = 1
            while not P.is_zero():
                n += 1
                P += self
            self.__order = rings.Integer(n)
        return self.__order




