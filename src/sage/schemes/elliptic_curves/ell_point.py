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

    def __list__(self):
        """
        Return the coordinates of this point as a list
        """
        return list(self._coords)

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

    def _bsgs(self, Q, lower, upper):
        """
        Implementation of baby-step-giant-step algorithm

        returns n with lower<=n<=upper such that n*self==Q, raising an
        error if non exists

        Mainly for internal use:  users should use discrete_log()

        AUTHOR: John Cremona, following his C++ implementation which
        in turn is closely based on LiDIA code

        EXAMPLE:
        sage: F=GF(3^6,'a')
        sage: a=F.gen()
        sage: E= EllipticCurve([0,1,1,a,a])
        sage: E.cardinality()
        762
        sage: A,G=E.abelian_group() ## set since this E is cyclic
        sage: P=G[0]
        sage: Q=400*P
        sage: P._bsgs(Q,1,E.cardinality())
        400
        """
        if not self.curve() == Q.curve():
            raise ValueError, "ecdlog requires points to be on the same curve"

        if lower<0 or upper<lower:
            raise ValueError, "_bsgs() requires 0<=lower<=upper"

        if self.is_zero() and not Q.is_zero():
            raise ValueError, "No solution in _bsgs()"

        if lower==0 and Q.is_zero():
            return rings.Integer(0)

#        print "In bsgs"
        ran = upper - lower
        if (ran < 30):    # for very small intervals
#            print "interval length=",ran," so using naive method"
            h=lower
            H=lower*self
            if H == Q: return rings.Integer(h);
            while h <= upper:
                h+=1
                H+=self
                if H == Q:  return rings.Integer(h)
            raise ValueError, "No solution in _bsgs()"

        # now we use the Babystep Giantstep algorithm

        # compute number of babysteps
        number_baby = 1+ran.isqrt()
        if number_baby > 3000000:  number_baby = 3000000
#        print "In bsgs: number_baby=",number_baby

        HT = dict()  # will store table of pairs (i*self,i)

        H2 = Q-lower*self
        H = self.curve()(0)

        # do the babysteps

        for i in range(1,number_baby+1):
            H+=self          # so H = i*self and H2 = Q-lower*self
            if H == H2:   # if i*self = Q-lower*self then solution = lower+i
                return rings.Integer(lower + i)
            if not H.is_zero(): # store [H,i] in table
                HT[H]=i
#        print "In bsgs: finished baby steps"

        # Now we have a table of pairs [i*self,i] for i in (1..number_baby)
        # and H  =  number_baby*self
        assert H==number_baby*self
        # and H2 =  Q-lower*self

        # We will subtract H from H2 repeatedly, so
        # H2=Q-lower*self-j*H in the loop

        # Giantsteps

        number_giant = 1+((upper - lower)//(number_baby)) # rounded!
#        print "In bsgs: number_giant=",number_giant

        step_size = number_baby;

        for j in range(number_giant+1):
            # Here H2=Q-(lower+j*step_size)*self
            if H2.is_zero(): # on the nail, no need to check table
                h = lower + j * step_size;
                if h <= upper:
                    return rings.Integer(h)
                else:
                    raise ValueError, "No solution in _bsgs()"

            # look in table to see if H2= i*self for a suitable i
            i = HT.get(H2)
            if not i==None:
                h = lower + i + j * step_size;
                if h <= upper:
                    return rings.Integer(h)
                else:
                    raise ValueError, "No solution in _bsgs()"
            H2-=H
        raise ValueError, "No solution in _bsgs()"

    def discrete_log(self, Q, ord=None):
        """
        Returns discrete log of Q with respect to self, i.e. an
        integer m with 0<=m<order(self) such that m*self==Q, if one
        exists; otherwise raise an error.  The order of self is
        computed if not supplied

        AUTHOR: John Cremona, following his C++ implementation which
        in turn is closely based on LiDIA code

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
            return self._bsgs(Q,0,ord-1)
        except:
            raise ValueError, "ECDLog problem has no solution"


    def order(self):
        """
        Return the order of this point on the elliptic curve.
        If the point has infinite order, returns 0.

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


        ALGORITHM: first finds a multiple of the order: either the
        group order, if this is known, or uses baby-step-giant-step
        algorithm.  If the group order is known, its factorization is
        cached.  From a multiple of the order and its factorization it
        is then easy to determine the exact order.  We must not cause
        the group order to be calculated when not known since this
        function is used in determining the group order via
        computation of several random points and their orders.

        AUTHOR: John Cremona, 2008-02-10.

        """
        try:
            return self._order
        except AttributeError:
            pass
        if self.is_zero():
            return rings.Integer(1)
        E = self.curve()
        K = E.base_ring()
        lb,ub = ell_generic.Hasse_bounds(K.order())

        try:
            M = E._order
            try:
                plist = E._prime_factors_of_order
            except:
                plist = M.prime_divisors()
                E._prime_factors_of_order = plist
        except:
            if K.is_prime_field():
                M = E.cardinality() # computed and cached
                plist = M.prime_divisors()
                E._prime_factors_of_order = plist
            else:
                M = self._bsgs(E(0),lb,ub)
                plist = M.prime_divisors()

        # Now M is a multiple of the order and plist is a list of
        # its prime factors

        # For each p in plist we determine the power of p dividing
        # the order, accumulating the order in N

        N=rings.Integer(1)
        for p in plist:
            Q=M.prime_to_m_part(p)*self   # so Q has p-power order
            while not Q.is_zero():
                Q=p*Q
                N*=p

        # now N is the exact order of self

        if 2*N>ub: # then we have a generator, so cache this
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

    # returns (m,a) where m>0 is minimal s.t. m*Q is in <P> with m*Q=a*P.
    # Special case: if <Q> and <P> are disjoint, then returns m=order(Q)
    # and a=0.

    def linear_relation(self, Q):
        """
        Function which solves the equation m*Q=a*self form (m,a) with
        minimal m>0.  Special case: if <self> and <Q> intersect only
        in {0} then (m,a)=(Q.order(),0).

        Used in determining group structure.  Applicable in general
        finite abelian groups: just uses the bsgs() function.

        EXAMPLE:
        sage: F.<a>=GF(3^6,'a')
        sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
        sage: P=E(a^5 + a^4 + a^3 + a^2 + a + 2 , 0)
        sage: Q=E(2*a^3 + 2*a^2 + 2*a , a^3 + 2*a^2 + 1)
        sage: P.linear_relation(Q)
        (2, 1)
        sage: 2*Q == P
        True
        """
        P = self
        n = P.order()
        m = Q.order()
        g=arith.gcd(n,m)
        if g==1: return (m,rings.Integer(0))
        n1=n//g
        m1=m//g
        P1=n1*P  # both of exact order g
        Q1=m1*Q  #
        # now see if Q1 is a multiple of P1; the only multiples we
        # need check are h*Q1 where h|g
        for h in g.divisors(): # positive divisors!
            try:
                a = P1._bsgs(h*Q1,0,g-1)
                a = a*n1
                m = h*m1
                assert m*Q==a*P        # debugging:
                return (m,a)
            except:
                pass # to next h
        raise ValueError, "No solution found in linear_relation!"

    # old version of order function

    def order_old(self):
        """
        Return the order of this point on the elliptic curve.
        If the point has infinite order, returns 0.

        EXAMPLE:
            sage: k.<a> = GF(5^5)
            sage: E = EllipticCurve(k,[2,4]); E
            Elliptic Curve defined by y^2  = x^3 + 2*x + 4 over Finite Field in a of size 5^5
            sage: P = E(3*a^4 + 3*a , 2*a + 1 )
            sage: P.order_old()
            3227
            sage: Q = E(0,2)
            sage: Q.order_old()
            7


        ALGORITHM: uses PARI's \code{ellzppointorder} if base ring is
        prime or baby-step-giant-step algorithm as presented in

        Washington, Lawrence C.; 'Elliptic Curves: Number Theory and
        Cryptography', Boca Raton 2003

        REMARKS (John Cremona, 2008-02-10): if we know the group order
        already then there is no need to use BSGS to compute a
        multiple of the order.  It would also be good to cache the
        factorization of the order.  But we don't want to rely on this
        since this function is used in the group order computation!

        """
        try:
            return self._order
        except AttributeError:
            pass
        if self.is_zero():
            return rings.Integer(1)
        E = self.curve()
        K = E.base_ring()
        if K.is_prime_field():
            e = E._gp()
            self._order = rings.Integer(e.ellzppointorder(list(self.xy())))
            return self._order
        else:
            P = self
            E = P.curve()
            k = E.base_ring()
            q = k.order()

            if q < 256: # TODO: check this heuristc
                n = 1
                while not P.is_zero():
                    n += 1
                    P += self
                self._order = rings.Integer(n)
                return self._order

            try:
                M=E._order
                try:
                    plist=E._prime_factors_of_order
                except:
                    plist = M.prime_divisors()
                    E._prime_factors_of_order=plist
            except:

                # 1. Compute Q = (q+1)P
                Q = (q+1) * P

                # 2. Choose an integer m with m > q^{1/4}. Compute and
                # store the points jP for j = 0,1,2,...,m

                m = rings.Integer((q**rings.RR(0.25)).floor() + 1) #

                l = dict()
                X = E(0)
                for j in range(0,m+1):
                    l[X] = j
                    X = P + X

                    # 3. Compute the points Q + k(2mP) for k = -m,
                    # -(m+1), ..., m until there is a match Q + k(2mP)
                    # = +- jP with a point or its negative on the list

                    twomP = (2*m*P)
                    for k in range(-m,m+1):
                        W =  Q + k*twomP
                        if W in l:
                            # 4a. Conclude that (q + 1 + 2mk - j)P =
                            # oo. Let M = q + 1 + 2mk - j
                            M = q + 1 + 2*m*k - l[W]
                            break
                        elif -W in l:
                            # 4b. Conclude that (q + 1 + 2mk + j)P =
                            # oo. Let M = q + 1 + 2mk + j
                            M = q + 1 + 2*m*k + l[-W]
                            break

                # 5. Now M is a multiple of the point's order and
                # plist is a list of the distinct prime factors of M.

                # Change by John Cremona 2008-02-10: avoid repeated
                # factorization of M

                plist = M.prime_divisors()

            # 6. Compute (M/p)P for p in plist. If (M/p)P = 0 for some
            # p, replace M with M/p and repeat. If (M/pi)P != 0 for
            # all p then M is the order of the point P.

            while True:
                N = M
                plist = filter(lambda p: p.divides(M), plist)
                for p in plist:
                    if (int(M/p) * P).is_zero():
                        M = M/p
                        break
                if N == M:
                    self._order = rings.Integer(M)
                    return self._order

