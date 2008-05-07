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

from math import ceil, floor, sqrt
from sage.structure.element import AdditiveGroupElement, RingElement

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
        This is synonmous with curve() which is perhaps more
        intuituve.

        Technically, points on curves in Sage are scheme maps from the
        domain Spec(F) where F is the base field of the curve to the
        codomain which is the curve.  See also domain() and codomain().

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
        on.  Synonymous with curve() which is perhaps more intuituve.

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
        Return the order of this point on the elliptic curve.  If the
        point has infinite order, returns 0.  This is only implemented
        here for curves defined over Q, by calling pari.  For curves
        over finite fields, see below.

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
        """
        Return the curve that this point is a point on.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: P.curve()
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return self.scheme()

    def __nonzero__(self):
        """
        Return True if this is the zero point on the curve.

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

    def is_finite_order(self):
        """
        Return True if this point has finite additive order as an element
        of the group of points on this curve.

        EXAMPLES:
            sage: E = EllipticCurve('11a')
            sage: P = E([5,5])
            sage: P.is_finite_order()
            True
            sage: Q = 5*P; Q
            (0 : 1 : 0)
            sage: Q.is_finite_order()
            True
            sage: E = EllipticCurve('37a')
            sage: P = E([0,0])
            sage: P.is_finite_order()
            False
        """
        return self[2] == 0 or self.order() != oo

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
            sage: Q = 5*E.1; Q
            (-2739/1444 : 22161/54872 : 1)
            sage: Q.is_divisible_by(4)
            False
            sage: Q.is_divisible_by(5)
            True
        """
        # For now we simply compute all the m division poins
        # and see if there are any.  There are many potential
        # tricks to speed this up that work usually but fall
        # down like crazy in corner cases.  If you really need
        # to check this condition much more quickly, maybe
        # copy some code out of the division_points function
        # and think hard about what you are doing.
        return len(self.division_points(m)) > 0

    def _short_division_polynomial(self, m):
        r"""
        Return the m-division polynomial of this point, which is a
        polynomial in x that captures all the x-coordinates of points
        Q such that m*Q = self ON THE SHORT WEIERSTRASS MODEL.  Except
        in rare cases involving cancellation, these are exactly the
        $x$-coordinates of the m-division points on the short model.

        WARNING: It is highly recommended that you read the source code for
        \code{self._division_points} before using the division polynomial.

        INPUT:
            m -- positive integer
        OUTPUT:
            a polynomial in a single variable x over the base field.

        EXAMPLES:
        This function is incredibly useful for, e.g., quickly deciding
        for numerous primes p whether a point has any m-th roots
        modulo p.
            sage: E = EllipticCurve('389a')
            sage: P = E.0
            sage: f = P._short_division_polynomial(3)
            sage: [p for p in primes(100) if len(f.change_ring(GF(p)).roots()) > 0]
            [2, 3, 7, 11, 13, 23, 31, 37, 43, 47, 61, 67, 71, 89, 97]
            sage: f = E.gen(1)._short_division_polynomial(3)
            sage: [p for p in primes(100) if len(f.change_ring(GF(p)).roots()) > 0]
            [2, 3, 7, 11, 13, 19, 23, 31, 37, 43, 47, 59, 61, 67, 71, 89, 97]
        """
        return self._division_points(m, poly_only=True)

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
        """
        return self._division_points(m, poly_only=False)

    def _division_points(self, m, poly_only):
        """
        Used internally by both division_points and
        division_polynomial.  Returns either all the m-division points
        of a point or if poly_only=True, return the polynomial whose
        roots 'determine' all the x-coordinates of m-division points.

        """
        # Coerce the input m to an integer
        m = rings.Integer(m)
        # Check for trivial case of m = 1.
        if m == 1:
            return [self]

        # Get the curve and convert it to short Weierstrass form so
        # that we can compute the multiplication by m map explicitly.
        # Also compute explicit maps back and forth.
        E = self.curve()
        F = E.short_weierstrass_model()
        F_to_E = F.isomorphism_to(E)
        E_to_F = E.isomorphism_to(F)
        f = F.multiplication_by_m(m, x_only=True)

        # Map self (our point) over to the short Weierstrass model.
        W = E_to_F(self)

        # ans will contain the list of division points.
        ans = []

        # We compute a polynomial g whose roots give all
        # possible x coordinates of the m-division points.
        # This polynomial can have even more roots in some
        # corner cases, but that's OK since we double check
        # ach below.
        if self == 0:
            # If self is the 0, then m*self = self automatically.
            ans.append(self)
            # Also the correct poly to choose is just the denominator.
            g = f.denominator()
        else:
            # The poly g below is 0 at x if f(x) = W[0].
            g = f.numerator() - W[0]*f.denominator()
        R = E.base_ring()['x']
        h = R(g)
        if poly_only:
            return h
        for x,_ in h.roots():
            # We use *no* special tricks or shortcuts here.
            # Every time I thought of one I discovered corner
            # cases where it didn't work. So we have none here
            # at all.
            try:
                # Make a point on the curve with this x coordinate.
                Z = F.lift_x(x)
                if m*Z == W:
                    # if it works, take it.
                    ans.append(F_to_E(Z))
                # Now try the other possible point with the same x coordinate.
                # Negation keeps the same x coordinate because F is in
                # short weierstrass form.
                nZ = -Z
                if m*nZ == W:
                    ans.append(F_to_E(nZ))
            except ValueError:
                # This can happen in some exceptional case.
                pass
        # Finally remove possible duplicates and sort.  (It's possible I guess
        # that there would be duplicates in some weird corner case.)
        ans = list(set(ans))
        ans.sort()
        return ans

    def height(self):
        """
        The Neron-Tate canonical height of the point.

        INPUT:
            self -- a point on a curve over QQ

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
        """
        if self.is_finite_order():
            return rings.QQ(0)
        h = self.curve().pari_curve().ellheight([self[0], self[1]])
        return rings.RR(h)

    def xy(self):
        """
        Return the x and y coordinates of this point, as a 2-tuple.
        If this is point at infinity a ZeroDivisionError is raised.

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


