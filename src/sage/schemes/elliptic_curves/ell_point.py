"""
Points on elliptic curves

EXAMPLES:
    sage: K = Qp(5)
    sage: E = EllipticCurve([K(1), K(1)])

    #sage: P = E([K(0), K(1), K(1)])
    #sage: P
    #(0 : 1 : 1)
    #sage: P + P
    #(4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + 3*5^10 + 3*5^11 + 3*5^12 + 3*5^13 + 3*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 3*5^18 + 3*5^19 + O(5^20) : 2 + 3*5^2 + 3*5^4 + 3*5^6 + 3*5^8 + 3*5^10 + 3*5^12 + 3*5^14 + 3*5^16 + 3*5^18 + O(5^20) : 1)

Arithmetic with a point over an extension of a finite field:
    sage: k.<a> = GF(5^2)
    sage: E = EllipticCurve(k,[1,0]); E
    Elliptic Curve defined by y^2  = x^3 + x over Finite Field in a of size 5^2
    sage: P = E([a,2*a+4])
    sage: 5*P
    (2*a + 3 : 2*a : 1)
    sage: P*5
    (2*a + 3 : 2*a : 1)
    sage: P + P + P + P + P
    (2*a + 3 : 2*a : 1)


    sage: F = Zmod(3)
    sage: E = EllipticCurve(F,[1,0]);
    sage: P = E([2,1])
    sage: import sys
    sage: n = sys.maxint
    sage: P*(n+1)-P*n == P
    True

AUTHORS:
   * William Stein (2005) -- Initial version
   * Robert Bradshaw et al....
   * John Cremona (Feb 2008) -- Point counting and group structure for
     non-prime fields, Frobenius endomorphism and order, elliptic logs
   * John Cremona (Aug 2008) -- Introduced EllipticCurvePoint_number_field class
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

from sage.structure.element import AdditiveGroupElement, RingElement
from sage.interfaces import gp
import sage.plot.all as plot

from sage.rings.padics.factory import Qp

import ell_generic
import sage.rings.all as rings
import sage.rings.arith as arith
import sage.misc.misc as misc
from sage.groups.all import AbelianGroup
import sage.groups.generic as generic

from sage.structure.sequence  import Sequence
from sage.schemes.generic.morphism import (SchemeMorphism_projective_coordinates_ring,
                                           SchemeMorphism_abelian_variety_coordinates_field,
                                           is_SchemeMorphism, SchemeMorphism_coordinates)

import sage.schemes.generic.scheme as scheme
from constructor import EllipticCurve

oo = rings.infinity       # infinity

class EllipticCurvePoint(SchemeMorphism_projective_coordinates_ring):
    """
    A point on an elliptic curve.
    """
    def __cmp__(self, other):
        """
        Standard comparison function for points on elliptic curves, to
        allow sorting and equality testing.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: P=E(0,1)
            sage: P.order()
            +Infinity
            sage: Q=P+P
            sage: P==Q
            False
            sage: Q+Q == 4*P
            True
        """
        if isinstance(other, (int, long, rings.Integer)) and other == 0:
            if self.is_zero():
                return 0
            else:
                return -1
        return SchemePoint_projective_abelian_scheme.__cmp__(self, other)

class EllipticCurvePoint_field(AdditiveGroupElement): # SchemeMorphism_abelian_variety_coordinates_field):
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
        Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        sage: K.<i>=NumberField(x^2+1)
        sage: E=EllipticCurve(K,[0,1,0,-160,308])
        sage: P=E(26,-120)
        sage: Q=E(2+12*i,-36+48*i)
        sage: P.order() == Q.order() == 4
        True
        sage: 2*P==2*Q
        False

        sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
        sage: E=EllipticCurve([0,0,0,0,t^2])
        sage: P=E(0,t)
        sage: P,2*P,3*P
        ((0 : t : 1), (0 : -t : 1), (0 : 1 : 0))


    TESTS:
        sage: loads(S.dumps()) == S
        True
        sage: E = EllipticCurve('37a')
        sage: P = E(0,0); P
        (0 : 0 : 1)
        sage: loads(P.dumps()) == P
        True
        sage: T = 100*P
        sage: loads(T.dumps()) == T
        True

    Test pickling an elliptic curve that has known points on it.
        sage: e = EllipticCurve([0, 0, 1, -1, 0]); g = e.gens(); loads(dumps(e)) == e
        True
    """
    def __init__(self, curve, v, check=True):
        """
        Constructor for a point on an elliptic curve

        INPUT:
            curve -- an elliptic curve
            v -- data determining a point (another point, the integer
                 0, or a tuple of coordinates)

        """
        point_homset = curve.point_homset()
        AdditiveGroupElement.__init__(self, point_homset)
        if check:
            # mostly from  SchemeMorphism_projective_coordinates_field
            d = point_homset.codomain().ambient_space().ngens()
            if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            if v == 0:
                v = (0,1,0)
            if not isinstance(v,(list,tuple)):
                raise TypeError, \
                      "Argument v (= %s) must be a scheme point, list, or tuple."%str(v)
            if len(v) != d and len(v) != d-1:
                raise TypeError, "v (=%s) must have %s components"%(v, d)
            v = Sequence(v, point_homset.value_ring())
            if len(v) == d-1:     # very common special case
                v.append(1)

            n = len(v)
            all_zero = True
            for i in range(n):
                if v[n-1-i]:
                    all_zero = False
                    c = v[n-1-i]
                    if c == 1:
                        break
                    for j in range(n-i):
                        v[j] /= c
                    break
            if all_zero:
                raise ValueError, "%s does not define a valid point since all entries are 0"%repr(v)

            point_homset.codomain()._check_satisfies_equations(v)

        self._coords = v


    def _repr_(self):
        """
        Return a string representation of this point
        """
        return self.codomain().ambient_space()._repr_generic_point(self._coords)

    def _latex_(self):
        """
        Return a latex representation of this point
        """
        return self.codomain().ambient_space()._latex_generic_point(self._coords)

    def __getitem__(self, n):
        """
        Return the n'th coordinate of this point
        """
        return self._coords[n]

    def __iter__(self):
        """
        Return the coordinates of this point as a list

        EXAMPLE:
            sage: E = EllipticCurve('37a')
            sage: list(E([0,0]))
            [0, 0, 1]
        """
        return iter(self._coords)

    def __tuple__(self):
        """
        Return the coordinates of this point as a tuple
        """
        return self._coords

    def __cmp__(self, other):
        """
        Comparison function for points to allow sorting and equality testing
        """
        if not isinstance(other, EllipticCurvePoint_field):
            try:
                other = self.codomain().ambient_space()(other)
            except TypeError:
                return -1
        return cmp(self._coords, other._coords)

    def scheme(self):
        """
        Return the scheme of this point, i.e. the curve it is on.
        This is synonymous with curve() which is perhaps more
        intuitive.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: P=E(0,1)
            sage: P.scheme()
            Elliptic Curve defined by y^2  = x^3 + x +1 over Rational Field
            sage: P.scheme() == P.curve()
            True
            sage: K.<a>=NumberField(x^2-3,'a')
            sage: P=E.base_extend(K)(1,a)
            sage: P.scheme()
            Elliptic Curve defined by y^2  = x^3 + x +1 over Number Field in a with defining polynomial x^2 - 3
        """
        #The following text is just not true: it applies to the class
        #EllipticCurvePoint, which appears to be never used, but does
        #not apply to EllipticCurvePoint_field which is simply derived
        #from AdditiveGroupElement.
        #
        #"Technically, points on curves in Sage are scheme maps from
        #  the domain Spec(F) where F is the base field of the curve to
        #  the codomain which is the curve.  See also domain() and
        #  codomain()."

        return self.codomain()

    def domain(self):
        """
        Return the domain of this point, which is Spec(F) where F is
        the field of definition.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: P=E(0,1)
            sage: P.domain()
            Spectrum of Rational Field
            sage: K.<a>=NumberField(x^2-3,'a')
            sage: P=E.base_extend(K)(1,a)
            sage: P.domain()
            Spectrum of Number Field in a with defining polynomial x^2 - 3
       """
        return self.parent().domain()

    def codomain(self):
        """
        Return the codomain of this point, which is the curve it is
        on.  Synonymous with curve() which is perhaps more intuitive.

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1])
            sage: P=E(0,1)
            sage: P.domain()
            Spectrum of Rational Field
            sage: K.<a>=NumberField(x^2-3,'a')
            sage: P=E.base_extend(K)(1,a)
            sage: P.codomain()
            Elliptic Curve defined by y^2  = x^3 + x +1 over Number Field in a with defining polynomial x^2 - 3
            sage: P.codomain() == P.curve()
            True
       """
        return self.parent().codomain()

    def order(self):
        """
        Return the order of this point on the elliptic curve.

        If the point is zero, returns 1, otherwise raise a
        NotImplementedError.

        For curves over number fields and finite fields, see below.

        EXAMPLE:
            sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
            sage: E=EllipticCurve([0,0,0,-t^2,0])
            sage: P=E(t,0)
            sage: P.order()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of order of a point not implemented over general fields.
        """
        if self.is_zero():
            return rings.Integer(1)
        raise NotImplementedError, "Computation of order of a point not implemented over general fields."

    def curve(self):
        """
        Return the curve that this point is on.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: P.curve()
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return self.scheme()

    def __nonzero__(self):
        """
        Return True if this is not the zero point on the curve.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: P = E(0); P
            (0 : 1 : 0)
            sage: P.is_zero()
            True
            sage: P = E.gens()[0]
            sage: P.is_zero()
            False
        """
        return bool(self[2])

    def has_finite_order(self):
        """
        Return True if this point has finite additive order as an element
        of the group of points on this curve.

        For fields other than number fields and finite fields, this is
        NotImplemented unless self.is_zero().

        EXAMPLES:
            sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
            sage: E=EllipticCurve([0,0,0,-t^2,0])
            sage: P = E(0)
            sage: P.has_finite_order()
            True
            sage: P=E(t,0)
            sage: P.has_finite_order()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of order of a point not implemented over general fields.
            sage: (2*P).is_zero()
            True
        """
        if self.is_zero(): return True
        return self.order() != oo

    is_finite_order = has_finite_order # for backward compatibility

    def has_infinite_order(self):
        """
        Return True if this point has infinite additive order as an element
        of the group of points on this curve.

        For fields other than number fields and finite fields, this is
        NotImplemented unless self.is_zero().

        EXAMPLES:
            sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
            sage: E=EllipticCurve([0,0,0,-t^2,0])
            sage: P = E(0)
            sage: P.has_infinite_order()
            False
            sage: P=E(t,0)
            sage: P.has_infinite_order()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of order of a point not implemented over general fields.
            sage: (2*P).is_zero()
            True
        """
        if self.is_zero(): return False
        return self.order() == oo

    def plot(self, **args):
        """
        Plot this point on an elliptic curve.

        INPUT:
            **args -- all arguments get passed directly onto the point
              plotting function.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: P.plot(pointsize=30, rgbcolor=(1,0,0))
        """
        if self.is_zero():
            return plot.text("$\\infty$", (-3,3), **args)

        else:
            return plot.point((self[0], self[1]), **args)

    def _add_(self, right):
        """
        Add self to right.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1]); Q = E([0,0])
            sage: P + Q
            (1 : 0 : 1)
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
            try:
                m = (3*x1*x1 + 2*a2*x1 + a4 - a1*y1) / (2*y1 + a1*x1 + a3)
            except ZeroDivisionError:
                raise ZeroDivisionError, "Inverse of %s does not exist"%(2*y1 + a1*x1 + a3)
        else:
            try:
                m = (y1-y2)/(x1-x2)
            except ZeroDivisionError:
                raise ZeroDivisionError, "Inverse of %s does not exist"%(x1-x2)

        x3 = -x1 - x2 - a2 + m*(m+a1)
        y3 = -y1 - a3 - a1*x3 + m*(x1-x3)
        return E.point([x3, y3, 1], check=False)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1]); Q = E([0,0])
            sage: P - Q
            (4 : 8 : 1)
            sage: (P - Q) + Q
            (-1 : 1 : 1)
            sage: P
            (-1 : 1 : 1)
        """
        return self + (-right)

    def __neg__(self):
        """
        Return the additive inverse of this point.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: Q = -P; Q
            (-1 : -2 : 1)
            sage: Q + P
            (0 : 1 : 0)
        """
        if self.is_zero():
            return self
        E, x, y = self.curve(), self[0], self[1]
        return E.point([x, -y - E.a1()*x - E.a3(), 1], check=False)

    def xy(self):
        """
        Return the x and y coordinates of this point, as a 2-tuple.
        If this is the point at infinity a ZeroDivisionError is raised.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: P.xy()
            (-1, 1)
            sage: Q = E(0); Q
            (0 : 1 : 0)
            sage: Q.xy()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Rational division by zero
        """
        if self[2] == 1:
            return self[0], self[1]
        else:
            return self[0]/self[2], self[1]/self[2]

    def is_divisible_by(self, m):
        """
        Return True if there exists a point Q defined over the same
        field as self such that m*Q == self.

        INPUT:
            m -- a positive integer
        OUTPUT:
            bool

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: Q = 5*E(0,0); Q
            (-2739/1444 : -77033/54872 : 1)
            sage: Q.is_divisible_by(4)
            False
            sage: Q.is_divisible_by(5)
            True
        """
        g = self.division_points(m, poly_only=True)
        return len(g.roots(multiplicities=False)) > 0

    def division_points(self, m, poly_only=False):
        """
        Return *all* points Q in the base field on the elliptic curve
        that contains self such that m*Q == self.

        INPUT:
            m -- a positive integer
            poly_only -- bool (default: False); if True return polynomial
                         whose roots give all possible x-coordinates of
                         m-th roots of self.
        OUTPUT:
            a (possibly empty) list

        EXAMPLES:
        We find the five 5-torsion points on an elliptic curve.
            sage: E = EllipticCurve('11a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: P = E(0); P
            (0 : 1 : 0)
            sage: P.division_points(5)
            [(0 : 1 : 0), (5 : -6 : 1), (5 : 5 : 1), (16 : -61 : 1), (16 : 60 : 1)]


        Note above that 0 is included since [5]*0 = 0.

        We create a curve of rank 1 with no torsion and do a consistency check:
            sage: E = EllipticCurve('11a').quadratic_twist(-7)
            sage: Q = E([44,-270])
            sage: (4*Q).division_points(4)
            [(44 : -270 : 1)]

        We create a curve over a non-prime finite field with group of order $18$:
            sage: k.<a> = GF(25)
            sage: E = EllipticCurve(k, [1,2+a,3,4*a,2])
            sage: P = E([3,3*a+4])
            sage: factor(E.order())
            2 * 3^2
            sage: P.order()
            9

        We find the $1$-division points as a consistency check -- there
        is just 1 of course:
            sage: P.division_points(1)
            [(3 : 3*a + 4 : 1)]

        The point $P$ has order coprime to $2$ but divisible by $3$ so:
            sage: P.division_points(2)
            [(2*a + 1 : 3*a + 4 : 1), (3*a + 1 : a : 1)]

        We check that each of the $2$ division points works as claimed:
            sage: [2*Q for Q in P.division_points(2)]
            [(3 : 3*a + 4 : 1), (3 : 3*a + 4 : 1)]

        Some other checks:
            sage: P.division_points(3)
            []
            sage: P.division_points(4)
            [(0 : 3*a + 2 : 1), (1 : 0 : 1)]
            sage: P.division_points(5)
            [(1 : 1 : 1)]

        An example over a number field (see trac #3383):
        sage: E = EllipticCurve('19a1')
        sage: K.<t> = NumberField(x^9-3*x^8-4*x^7+16*x^6-3*x^5-21*x^4+5*x^3+7*x^2-7*x+1)
        sage: EK = E.base_extend(K)
        sage: E(0).division_points(3)
        [(0 : 1 : 0), (5 : -10 : 1), (5 : 9 : 1)]
        sage: EK(0).division_points(3)
        [(0 : 1 : 0), (5 : 9 : 1), (5 : -10 : 1)]
        sage: E(0).division_points(9)
        [(0 : 1 : 0), (5 : -10 : 1), (5 : 9 : 1)]
        sage: EK(0).division_points(9)
        [(0 : 1 : 0), (5 : 9 : 1), (5 : -10 : 1), (-150/121*t^8 + 414/121*t^7 + 1481/242*t^6 - 2382/121*t^5 - 103/242*t^4 + 629/22*t^3 - 367/242*t^2 - 1307/121*t + 625/121 : 35/484*t^8 - 133/242*t^7 + 445/242*t^6 - 799/242*t^5 + 373/484*t^4 + 113/22*t^3 - 2355/484*t^2 - 753/242*t + 1165/484 : 1), (-150/121*t^8 + 414/121*t^7 + 1481/242*t^6 - 2382/121*t^5 - 103/242*t^4 + 629/22*t^3 - 367/242*t^2 - 1307/121*t + 625/121 : -35/484*t^8 + 133/242*t^7 - 445/242*t^6 + 799/242*t^5 - 373/484*t^4 - 113/22*t^3 + 2355/484*t^2 + 753/242*t - 1649/484 : 1), (-1383/484*t^8 + 970/121*t^7 + 3159/242*t^6 - 5211/121*t^5 + 37/484*t^4 + 654/11*t^3 - 909/484*t^2 - 4831/242*t + 6791/484 : 927/121*t^8 - 5209/242*t^7 - 8187/242*t^6 + 27975/242*t^5 - 1147/242*t^4 - 1729/11*t^3 + 1566/121*t^2 + 12873/242*t - 10871/242 : 1), (-1383/484*t^8 + 970/121*t^7 + 3159/242*t^6 - 5211/121*t^5 + 37/484*t^4 + 654/11*t^3 - 909/484*t^2 - 4831/242*t + 6791/484 : -927/121*t^8 + 5209/242*t^7 + 8187/242*t^6 - 27975/242*t^5 + 1147/242*t^4 + 1729/11*t^3 - 1566/121*t^2 - 12873/242*t + 10629/242 : 1), (-4793/484*t^8 + 6791/242*t^7 + 10727/242*t^6 - 18301/121*t^5 + 2347/484*t^4 + 2293/11*t^3 - 7311/484*t^2 - 17239/242*t + 26767/484 : 30847/484*t^8 - 21789/121*t^7 - 34605/121*t^6 + 117164/121*t^5 - 10633/484*t^4 - 29437/22*t^3 + 39725/484*t^2 + 55428/121*t - 176909/484 : 1), (-4793/484*t^8 + 6791/242*t^7 + 10727/242*t^6 - 18301/121*t^5 + 2347/484*t^4 + 2293/11*t^3 - 7311/484*t^2 - 17239/242*t + 26767/484 : -30847/484*t^8 + 21789/121*t^7 + 34605/121*t^6 - 117164/121*t^5 + 10633/484*t^4 + 29437/22*t^3 - 39725/484*t^2 - 55428/121*t + 176425/484 : 1)]

        """
        # Coerce the input m to an integer
        m = rings.Integer(m)
        # Check for trivial cases of m = 1, -1 and 0.
        if m == 1 or m == -1:
            return [self]
        if m == 0:
            if self == 0: # then every point Q is a solution, but...
                return [self]
            else:
                return []

        # ans will contain the list of division points.
        ans = []

        # We compute a polynomial g whose roots give all possible x
        # coordinates of the m-division points.  The number of
        # solutions (over the algebraic closure) is m^2, assuming that
        # the characteristic does not divide m.

        E = self.curve()
        P = self
        nP = -P
        P_is_2_torsion = (P==nP)

        # If self is the 0, then self is a solution, and the correct
        # poly is the m'th division polynomial
        if P == 0:
            ans.append(P)
            g = E.division_polynomial(m)
        else:
            # The poly g here is 0 at x(Q) iff x(m*Q) = x(P).
            g = E._multiple_x_numerator(m) - P[0]*E._multiple_x_denominator(m)

            # When 2*P=0, then -Q is a solution iff Q is.  For even m,
            # no 2-torsion point is a solution, so that g is the
            # square of a polynomial g1 of degree m^2/2, and each root
            # of g1 leads to a pair of solutions Q, -Q to m*Q=P.  For
            # odd m, P itself is the only 2-torsion solution, so g has
            # the form (x-x(P))*g1(x)^2 where g1 has degree (m^2-1)/2
            # and each root of g1 leads to a pair Q, -Q.

            if P_is_2_torsion:
                if m%2==0:
                    # This computes g.sqrt() which is not implemented
                    g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()

            # When 2*P!=0, then for each solution Q to m*Q=P, -Q is
            # not a solution (and points of order 2 are not
            # solutions).  Hence the roots of g are distinct and each
            # gives rise to precisely one solution Q.

                else:
                    g0 = g.variables()[0] - P[0]
                    g = g // g0
                    g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()
                    g = g0*g

        if poly_only:
            return g

        for x in g.roots(multiplicities=False):
            if E.is_x_coord(x):
                # Make a point on the curve with this x coordinate.
                Q = E.lift_x(x)
                nQ = -Q
                mQ = m*Q
                # if P==-P then Q works iff -Q works, so we include
                # both unless they are equal:
                if P_is_2_torsion:
                    if mQ == P:
                        ans.append(Q)
                        if nQ != Q:
                            ans.append(nQ)
                else:
                    # P is not 2-torsion so at most one of Q, -Q works
                    # and we must try both:
                    if mQ == P:
                        ans.append(Q)
                    elif mQ == nP:
                        ans.append(nQ)

        # Finally, sort and return
        ans.sort()
        return ans

    def _divide_out(self,p):
        """
        Return (Q,k) where p^k*Q==self and Q cannot be divided by p.

        WARNING: It is up to the caller to make sure that this does
        not loop endlessly.  It is used in
        EllipticCurve_generic._p_primary_torsion_basis(), when self
        will always have (finite) order which is a power of p, so that
        the order of Q increases by a factor of p at each stage.

        Since it will clearly be in danger of looping when
        self.is_zero(), this case is caught, but otherwise caveat
        user.

        EXAMPLES:
            sage: E = EllipticCurve('37a1')
            sage: P = E([0, 0])
            sage: R = 12*P
            sage: R._divide_out(2)
            ((-1 : -1 : 1), 2)
            sage: R._divide_out(3)
            ((2 : -3 : 1), 1)
            sage: R._divide_out(5)
            ((1357/841 : 28888/24389 : 1), 0)
            sage: R._divide_out(12)
            Traceback (most recent call last):
            ...
            ValueError: p (=12) should be prime.
        """
        p = rings.Integer(p)
        if not p.is_prime():
            raise ValueError, "p (=%s) should be prime."%p

        if self.is_zero():
            raise ValueError, "self must not be 0."

        k=0; Q=self
        pts = Q.division_points(p)
        while len(pts) > 0:
            Q = pts[0]
            k += 1
            pts = Q.division_points(p)
        return (Q,k)


    ##############################  end  ################################

class EllipticCurvePoint_number_field(EllipticCurvePoint_field):
    """
    A point on an elliptic curve over a number field.

    Most of the functionality is derived from the parent class
    EllipticCurvePoint_field.  In addition we have support for the
    order of a point, and heights (curently only implemented over QQ).

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
        Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    TESTS:
        sage: loads(S.dumps()) == S
        True
        sage: P = E(0,0); P
        (0 : 0 : 1)
        sage: loads(P.dumps()) == P
        True
        sage: T = 100*P
        sage: loads(T.dumps()) == T
        True

    Test pickling an elliptic curve that has known points on it.
        sage: e = EllipticCurve([0, 0, 1, -1, 0]); g = e.gens(); loads(dumps(e)) == e
        True
    """

    def order(self):
        """
        Return the order of this point on the elliptic curve.

        If the point has infinite order, returns +Infinity.  For
        curves defined over Q, we call pari.  For curves over other
        number fields we implement this here.

        EXAMPLE:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.order()
            +Infinity

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.order()
            2
        """
        try:
            return self._order
        except AttributeError:
            pass

        if self.is_zero():
            self._order = rings.Integer(1)
            return self._order

        E = self.curve()

        # Special code for curves over Q, calling pari
        from sage.libs.pari.gen import PariError
        try:
            n = int(E.pari_curve().ellorder([self[0], self[1]]))
            if n == 0: n = oo
            self._order = n
            return n
        except PariError:
            pass

        # Get the torsion order if known, else a bound on (multiple
        # of) the order.  We do not compute the torsion if it is not
        # already known, since computing the bound is faster (and is
        # also cached).

        try:
            N = E._torsion_order
        except AttributeError:
            N = E._torsion_bound()

        # Now self is a torsion point iff it is killed by N:
        if not (N*self).is_zero():
            self._order = oo
            return self._order

        # Finally we find the exact order using the generic code:
        self._order = generic.order_from_multiple(self,N,operation='+')
        return self._order

    def has_finite_order(self):
        """
        Return True iff this point has finite order on the elliptic curve.

        EXAMPLE:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.has_finite_order()
            False

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.has_finite_order()
            True
        """
        if self.is_zero(): return True
        return self.order() != oo

    def has_infinite_order(self):
        """
        Return True iff this point has infinite order on the elliptic curve.

        EXAMPLE:
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.has_infinite_order()
            True

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.has_infinite_order()
            False
        """
        if self.is_zero(): return False
        return self.order() == oo

    def is_on_identity_component(self, embedding=None):
        """
        Returns True iff this point is on the identity component of
        its curve with respect to a given (real or complex) embedding.

        INPUT:
            self -- a point on a curve over any ordered field (e.g. QQ)

            embedding -- an embedding from the base_field of the
                         point's curve into RR or CC; if None (the
                         default) it uses the first embedding of the
                         base_field into RR if any, else the first
                         embedding into CC.

        OUTPUT:
            True iff the point is on the identity component of the
            curve.  (If the point is zero then the result is True.)

        EXAMPLES:
        For K=QQ there is no need to specify an embedding:
            sage: E=EllipticCurve('5077a1')
            sage: [E.lift_x(x).is_on_identity_component() for x in range(-3,5)]
            [False, False, False, False, False, True, True, True]

        An example over a field with two real embeddings:
            sage: L.<a> = QuadraticField(2)
            sage: E=EllipticCurve(L,[0,1,0,a,a])
            sage: P=E(-1,0)
            sage: [P.is_on_identity_component(e) for e in L.embeddings(RR)]
            [False, True]

        We can check this as follows:
            sage: [e(E.discriminant())>0 for e in L.embeddings(RR)]
            [True, False]
            sage: e = L.embeddings(RR)[0]
            sage: E1 = EllipticCurve(RR,[e(ai) for ai in E.ainvs()])
            sage: e1,e2,e3 = E1.two_division_polynomial().roots(RR,multiplicities=False)
            sage: e1 < e2 < e3 and e(P[0]) < e3
            True
        """
        if self.is_zero():       # trivial case
            return True

        e = embedding
        # It is also trivially true if we have a complex embedding
        if not e is None:
            if not rings.is_RealField(e.codomain()):
                return True

        # find a suitable embedding if none was supplied:
        E = self.curve()
        K = E.base_field()
        if e is None:
            try:
                e = K.embeddings(rings.RealField())[0]
            except IndexError:
                e = K.embeddings(rings.ComplexField())[0]

        # If there is only one component, the result is True:
        if not rings.is_RealField(e.codomain()): # complex embedding
            return True
        if e(E.discriminant()) < 0: # only one component
            return True

        # Now we have a real embedding and two components and have to work:
        gx = E.two_division_polynomial()
        gxd = gx.derivative()
        gxdd = gxd.derivative()
        return ( e(gxd(self[0])) > 0 and e(gxdd(self[0])) > 0)

    def height(self, precision=None):
        """
        The Neron-Tate canonical height of the point.

        Currently only implemented for curves defined over QQ (using
        pari); and for points of finite order.

        INPUT:
            self -- a point on a curve over QQ

            precision -- int or None (default: None): the precision
                         in bits of the result (default real precision
                         if None)

        OUTPUT:
            the rational number 0 or a nonzero real field element

        EXAMPLES:
            sage: E = EllipticCurve('11a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: P = E([5,5]); P
            (5 : 5 : 1)
            sage: P.height()
            0
            sage: Q = 5*P
            sage: Q.height()
            0

            sage: E = EllipticCurve('37a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: P = E([0,0])
            sage: P.height()
            0.0511114082399688
            sage: P.order()
            +Infinity
            sage: E.regulator()      # slightly random output
            0.051111408239968840

            sage: E = EllipticCurve('4602a1'); E
            Elliptic Curve defined by y^2 + x*y  = x^3 + x^2 - 37746035*x - 89296920339 over Rational Field
            sage: x = 77985922458974949246858229195945103471590
            sage: y = 19575260230015313702261379022151675961965157108920263594545223
            sage: d = 2254020761884782243
            sage: E([ x / d^2,  y / d^3 ]).height()
            86.7406561381275

            sage: E = EllipticCurve([17, -60, -120, 0, 0]); E
            Elliptic Curve defined by y^2 + 17*x*y - 120*y = x^3 - 60*x^2 over Rational Field
            sage: E([30, -90]).height()
            0

            sage: E = EllipticCurve('389a1'); E
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
            sage: [P,Q] = [E(-1,1),E(0,-1)]
            sage: P.height(precision=100)
            0.68666708330558658572355210295
            sage: (3*Q).height(precision=100)/Q.height(precision=100)
            9.0000000000000000000000000000
            sage: _.parent()
            Real Field with 100 bits of precision

        Unfortunately, canonical height isn't implemented in general.
            sage: E = EllipticCurve('5077a1').change_ring(QuadraticField(-3,'a'))
            sage: P = E([-2,3,1])
            sage: P.height()
            Traceback (most recent call last):
            ...
            NotImplementedError: canonical height not yet implemented over general number fields.
        """
        if self.has_finite_order():
            return rings.QQ(0)

        if precision is None:
            precision = rings.RealField().precision()

        try:
            h = self.curve().pari_curve(prec=precision).ellheight([self[0], self[1]])
            return rings.RealField(precision)(h)
        except AttributeError:
            raise NotImplementedError, "canonical height not yet implemented over general number fields."

    def elliptic_logarithm(self, embedding=None, precision=100, algorithm='pari'):
        """
        Returns the elliptic logarithm of this elliptic curve point.

        An embedding of the base field into RR (with any precision)
        may be given; otherwise the first real embedding is used (with
        the specified precision), if there are any; otherwise a
        NotImplementedError is raised, since we have not yet
        implemented the complex elliptic logarithm.

        INPUT:
            embedding: an embedding of the base field into RR
            precision: a positive integer (default 100) setting the
                       number of bits of precision for the computation
            algorithm: either 'pari' (default) to use Pari's
                       ellpointtoz, or 'sage' for a native
                       implementation (currently not very accurate)

        ALGORITHM: See -[Co2] Cohen H., A Course in Computational
                        Algebraic Number Theory GTM 138, Springer 1996

        AUTHORS:
            - Michael Mardaus (2008-07) }
            - Tobias Nagel (2008-07)    } original version from [Co2]
            - John Cremona (2008-07)    revision following eclib code

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: E.discriminant() > 0
            True
            sage: P = E([-1,1])
            sage: P.is_on_identity_component ()
            False
            sage: P.elliptic_logarithm (precision=96)
            0.4793482501902193161295330101 + 0.985868850775824102211203849...*I
            sage: Q=E([3,5])
            sage: Q.is_on_identity_component()
            True
            sage: Q.elliptic_logarithm (precision=96)
            1.931128271542559442488585220

        An example with negative discriminant, and a torsion point:
            sage: E = EllipticCurve('11a1')
            sage: E.discriminant() < 0
            True
            sage: P = E([16,-61])
            sage: P.elliptic_logarithm (precision=96)
            0.2538418608559106843377589234
            sage: E.period_lattice().basis(prec=96)[0] / P.elliptic_logarithm (precision=96)
            5.000000000000000000000000000

        A larger example, illustrating some precision issues.  The
        default algorithm uses Pari and makes sure the result has the
        requested precision.  Note that Pari may complain if the
        requested precision is too low:
            sage: E = EllipticCurve([1, 0, 1, -85357462, 303528987048]) #18074g1
            sage: P = E([4458713781401/835903744, -64466909836503771/24167649046528, 1])
            sage: P.elliptic_logarithm(precision=64)
            Traceback (most recent call last):
            ...
            PariError: precision too low (18)
            sage: P.elliptic_logarithm(precision=65)
            0.2765620401410706146
            sage: P.elliptic_logarithm()  # 100 bits
            0.27656204014107061464076203097

        Note that the native algorithm 'sage' has trouble with
        precision in this example:
            sage: P.elliptic_logarithm(algorithm='sage')  # 100 bits
            0.2765620401410710087007...
        """
        from sage.rings.number_field.number_field import refine_embedding

        emb = embedding

        # find a suitable embedding if none was supplied:
        E = self.curve()
        K = E.base_field()
        if emb is None:
            try:
                emb = K.embeddings(rings.RealField(precision))[0]
            except IndexError:
                raise NotImplementedError, "elliptic logarithm not yet implemented for complex embeddings."
        else:
        # Check if we have been given a complex embedding
            if not rings.is_RealField(emb.codomain()):
                raise NotImplementedError, "elliptic logarithm not yet implemented for complex embeddings."
            else:
                # Get the precision of the supplied embedding
                prec = emb.codomain().precision()
                # if the precision parameter is greater, refine the embedding:
                if precision > prec:
                    emb = refine_embedding(emb,precision)

        # From now on emb() is a real embedding of K into RealField(precision)
        RR = rings.RealField(precision)
        # Check the trivial case: we do not put this earlier since we
        # want to return zero to the correct precision.
        if self.is_zero():
            return rings.ComplexField(precision)(0)

        #Initialize

        if algorithm == 'pari':
            from sage.libs.pari.all import pari
            from sage.libs.pari.gen import prec_words_to_bits
            E = self.curve()
            ai = [emb(a) for a in E.a_invariants()]
            ER = EllipticCurve(ai) # defined over RR
            E_pari = ER.pari_curve(prec=precision)
            pt_pari = [pari(emb(self[0])), pari(emb(self[1]))]
            log_pari = E_pari.ellpointtoz(pt_pari, precision=precision)
            while prec_words_to_bits(log_pari.precision()) < precision:
                working_prec = 2*precision
                emb = K.embeddings(rings.RealField(working_prec))[0]
                ai = [emb(a) for a in E.a_invariants()]
                ER = EllipticCurve(ai) # defined over RR
                E_pari = ER.pari_curve(prec=working_prec)
                pt_pari = [pari(emb(self[0])), pari(emb(self[1]))]
                log_pari = E_pari.ellpointtoz(pt_pari, precision=working_prec)
            C = rings.ComplexField(precision)
            return C(log_pari)

        if algorithm <> 'sage':
            raise ValueError, "algorithm must be either 'pari' or 'sage'"

        E = self.curve()
        ai = [emb(a) for a in E.a_invariants()]
        ER = EllipticCurve(ai) # defined over RR
        real_roots = ER.two_division_polynomial().roots(RR,multiplicities=False)
        xP, yP = self[0], self[1]
        wP = 2*yP+E.a1()*xP+E.a3()
        connected = emb(E.discriminant()) < 0

        if connected: #Connected Case
            # Here we use formulae equivalent to those in Cohen, but better
            # behaved when roots are close together

            # For some curves (e.g. 3314b3) the default precision is not enough!
            while len(real_roots)!=1:
                precision*=2
                RR=rings.RealField(precision)
                emb = refine_embedding(emb,precision)
                ai = [emb(a) for a in E.a_invariants()]
                ER = EllipticCurve(ai) # defined over RR
                real_roots = ER.two_division_polynomial().roots(RR,multiplicities=False)
            try:
                assert len(real_roots) == 1
            except:
                raise ValueError, ' none or more than one real root despite disc < 0'
            x,y,w = emb(xP), emb(yP), emb(wP)
            e1 = real_roots[0]
            roots = ER.two_division_polynomial().roots(rings.ComplexField(precision),multiplicities=False)
            roots.remove(e1)
            e2,e3 = roots
            a1,a2,a3 = ai[:3]

            zz = (e1-e2).sqrt() # complex
            beta = (e1-e2).abs()
            a = 2*zz.abs()
            b = 2*zz.real();
            c = (x-e1+beta)/((x-e1).sqrt())
            while (a - b)/a > 0.5**(precision-1):
                a,b,c = (a + b)/2, (a*b).sqrt(), (c + (c**2 + b**2 - a**2).sqrt())/2
            z = (a/c).arcsin()
            if w*((x-e1)*(x-e1)-beta*beta) >= 0:
                z = RR.pi() - z
            if w>0:
                z = z + RR.pi()
            z /= a
            return z

        else:                    #Disconnected Case, disc > 0
            # For some curves (e.g. 2370i5) the default precision is not enough!
            while len(real_roots)!=3:
                precision*=2
                RR=rings.RealField(precision)
                emb = refine_embedding(emb,precision)
                ai = [emb(a) for a in E.a_invariants()]
                ER = EllipticCurve(ai) # defined over RR
                real_roots = ER.two_division_polynomial().roots(RR,multiplicities=False)

            real_roots.sort() # increasing order
            real_roots.reverse() # decreasing order e1>e2>e3
            try:
                assert len(real_roots) == 3
            except:
                raise ValueError, ' none or more than one real root despite disc < 0'
            e1, e2, e3 = real_roots
            w1, w2 = E.period_lattice(emb).basis(precision)
            x,y,w = emb(xP), emb(yP), emb(wP)
            on_egg = (x < e1)
            a1,a2,a3 = ai[:3]

            # if P is on the "egg" (connected component), replace it by P+T3
            # where T3=(e3,y3) is a 2-torsion point on
            # the egg coming from w2/2 on the lattice

            if on_egg:
                y3 = -(a1*e3+a3)/2
                lam = (y-y3)/(x-e3)
                x3 = lam*(lam+a1)-a2-x-e3
                y = lam*(x-x3)-y-a1*x3-a3
                x = x3
                w = 2*y+a1*x+a3

            a = (e1 - e3).sqrt()
            b = (e1 - e2).sqrt()
            c = (x - e3).sqrt()
            while (a - b)/a > 0.5**(precision-1):
                a,b,c = (a + b)/2, (a*b).sqrt(), (c + (c**2 + b**2 - a**2).sqrt())/2

            z = (a/c).arcsin()/a
            if w > 0:
                z = w1 - z
            if on_egg:
                z = z + w2/2
            return z

    ##############################  end  ################################

class EllipticCurvePoint_finite_field(EllipticCurvePoint_field):

    def _magma_init_(self):
        """
        Return a string representation of self that MAGMA can
        understand.
        """
        E = self.curve()._magma_().name()
        x,y = map(lambda x: x._magma_().name(), self.xy())
        return "%s![%s,%s]"%(E,x,y)

    def discrete_log(self, Q, ord=None):
        """
        Returns discrete log of Q with respect to self, i.e. an
        integer m with 0<=m<order(self) such that m*self==Q, if one
        exists; otherwise raise an error.  The order of self is
        computed if not supplied

        AUTHOR: John Cremona. Adapted to use generic functions 2008-04-05

        EXAMPLE:
        sage: F=GF(3^6,'a')
        sage: a=F.gen()
        sage: E= EllipticCurve([0,1,1,a,a])
        sage: E.cardinality()
        762
        sage: A,G=E.abelian_group() ## set since this E is cyclic
        sage: P=G[0]
        sage: Q=400*P
        sage: P.discrete_log(Q)
        400
        """
        if ord==None: ord=self.order()
        try:
            return generic.discrete_log(Q,self,ord,operation='+')
        except:
            raise ValueError, "ECDLog problem has no solution"


    def order(self):
        """
        Return the order of this point on the elliptic curve.

        EXAMPLES:
            sage: k.<a> = GF(5^5)
            sage: E = EllipticCurve(k,[2,4]); E
            Elliptic Curve defined by y^2  = x^3 + 2*x + 4 over Finite Field in a of size 5^5
            sage: P = E(3*a^4 + 3*a , 2*a + 1 )
            sage: P.order()
            3227
            sage: Q = E(0,2)
            sage: Q.order()
            7

            In the next example, the cardinality of E will be computed
            (using SEA) and cached:

            sage: p=next_prime(2^150)
            sage: E=EllipticCurve(GF(p),[1,1])
            sage: P=E(831623307675610677632782670796608848711856078, 42295786042873366706573292533588638217232964)
            sage: P.order()
            1427247692705959881058262545272474300628281448
            sage: P.order()==E.cardinality()
            True


        ALGORITHM: Use generic functions.  If the group order is
        known, use order_from_multiple(), otherwise use
        order_from_bounds() with the Hasse bounds.  In the latter case
        we might find that we have a generator, in which case it is
        cached.

        We do not cause the group order to be calculated when not
        known, since this function is used in determining the group
        order via computation of several random points and their
        orders.

        AUTHOR: John Cremona, 2008-02-10, adapted 2008-04-05 to use
        generic functions

        """
        try:
            return self._order
        except AttributeError:
            pass
        if self.is_zero():
            return rings.Integer(1)
        E = self.curve()
        K = E.base_ring()
        bounds = ell_generic.Hasse_bounds(K.order())

        try:
            M = E._order
            try:
                plist = E._prime_factors_of_order
            except:
                plist = M.prime_divisors()
                E._prime_factors_of_order = plist
            N = generic.order_from_multiple(self,M,plist,operation='+')
        except:
            if K.is_prime_field():
                M = E.cardinality() # computed and cached
                plist = M.prime_divisors()
                E._prime_factors_of_order = plist
                N = generic.order_from_multiple(self,M,plist,operation='+')
            else:
                N = generic.order_from_bounds(self,bounds,operation='+')

        if 2*N>bounds[1]: # then we have a generator, so cache this
            try:
                dummy = E._order
            except:
                E._order = N
            try:
                dummy = E.__abelian_group
            except:
                E.__abelian_group = AbelianGroup([N]), (self,)

        self._order = N
        return self._order

