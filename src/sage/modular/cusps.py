r"""
The set $\PP^1(\Q)$ of cusps

EXAMPLES:
    sage: Cusps()
    The set P^1(Q) of cusps.

    sage: Cusp(oo)
    Infinity
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

from sage.rings.all import infinity, is_Infinity
from sage.rings.integer_ring import IntegerRing
from sage.rings.rational_field import RationalField
import sage.rings.rational as rational
import sage.rings.integer as integer
from sage.ext.sage_object import SageObject

ZZ = IntegerRing(); QQ = RationalField()

class Cusps(SageObject):
    """
    The set of cusps.

    EXAMPLES:
        sage: C = Cusps(); C
        The set P^1(Q) of cusps.
        sage: loads(C.dumps()) == C
        True
    """
    def __init__(self):
        pass

    def __cmp__(self, right):
        if isinstance(right, Cusps):
            return 0
        return -1

    def __repr__(self):
        return "The set P^1(Q) of cusps."

    def __call__(self, x):
        if isinstance(x, Cusp):
            return x
        return Cusp(x)


class Cusp(SageObject):
    """
    A cusp.

    A cusp is either a rational number or infinity, i.e., an element
    of the projective line over Q.  A Cusp is stored as a pair (a,b),
    where gcd(a,b)=1 and a,b are of type Integer.
    """

    def __init__(self, a, b=ZZ(1), construct=False):
        r"""
        Create the cusp a/b in $\PP^1(\Q)$, where if b=0 this is the
        cusp at infinity.

        EXAMPLES:
            sage: Cusp(2,3)
            2/3
            sage: Cusp(3,6)
            1/2
            sage: Cusp(1,0)
            Infinity
            sage: Cusp(infinity)
            Infinity
            sage: Cusp(5)
            5
            sage: Cusp(QQ("1/2"))             # rational number
            1/2
            sage: Cusp('2/3')             # rational number
            2/3
            sage: Cusp(1.5)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce 1.50,1 to a Cusp


            sage: a = Cusp(2,3)
            sage: loads(a.dumps()) == a
            True
        """
        if construct:
            self.__a = a; self.__b = b
            return

        if is_Infinity(a):
            self.__a = ZZ(1)
            self.__b = ZZ(0)
            return

        elif isinstance(a, Cusp):
            self.__a = a.__a
            self.__b = a.__b
            return

        elif isinstance(a, rational.Rational):
            a = a/b
            self.__a = a.numer()
            self.__b = a.denom()
            return

        elif isinstance(a, (int, long, integer.Integer)) and \
                 isinstance(b, (int, long, integer.Integer)):
            a = ZZ(a)
            b = ZZ(b)

        elif isinstance(a, integer.Integer) and isinstance(b, integer.Integer):
            pass

        elif isinstance(a, str):
            a = RationalField()(a)/b
            self.__a = a.numer()
            self.__b = a.denom()
            return

        else:

            raise TypeError, "Unable to coerce %s,%s to a Cusp"%(a,b)


        # Now a, b are both of type ZZ.
        if b < 0:
            b *= -1
            a *= -1
        g = a.gcd(b)
        self.__a = a//g
        self.__b = b//g
        return

    def __cmp__(self, right):
        if not isinstance(right, Cusp):
            return -1
        return cmp([self.__a, self.__b], [right.__a, right.__b])

    def is_infinity(self):
        """
        Returns True if this is the cusp infinity.
        """
        return self.__b == 0

    def numerator(self):
        """
        Return the numerator of the cusp a/b.

        sage: x=Cusp(6,9); x
        2/3
        sage: x.numerator()
        2
        """
        return self.__a

    def denominator(self):
        """
        Return the denominator of the cusp a/b.

        sage: x=Cusp(6,9); x
        2/3
        sage: x.denominator()
        3
        """
        return self.__b

    def _rational_(self):
        return self.__a / self.__b

    def _integer_(self):
        if self.__b != 1:
            raise TypeError, "cusp %s is not an integer"%self
        return self.__a

    def parent(self):
        """
        Return the set of all cusps.
        """
        return Cusps()

    def __repr__(self):
        if self.__b.is_zero():
            return "Infinity"
        if self.__b != 1:
            return "%s/%s"%(self.__a,self.__b)
        else:
            return str(self.__a)

    def _latex_(self):
        if self.__b.is_zero():
            return "\\infty"
        if self.__b != 1:
            return "\\frac{%s}{%s}"%(self.__a,self.__b)
        else:
            return str(self.__a)

    def __neg__(self):
        return Cusp(-self.__a, self.__b)

    def is_gamma0_equiv(self, other, N, transformation = False):
        """
        Return whether self and other are equivalent modulo the action
        of Gamma_0(N) via linear fractional transformations.

        INPUT:
            other -- Cusp
            N -- an integer (specifies the group Gamma_0(N))

            transformation -- bool (default: False), if True, also
                              return upper left entry of a matrix in
                              Gamma_0(N) that sends self to other.

        OUTPUT:
            bool -- True if self and other are equivalent
            integer -- returned only if transformation is True

        EXAMPLES:
            sage: x = Cusp(2,3)
            sage: y = Cusp(4,5)
            sage: x.is_gamma0_equiv(y, 2)
            True
            sage: x.is_gamma0_equiv(y, 2, True)
            (True, 1)
            sage: x.is_gamma0_equiv(y, 3)
            False
            sage: x.is_gamma0_equiv(y, 3, True)
            (False, None)
            sage: Cusp(1,0)
            Infinity
            sage: z = Cusp(1,0)
            sage: x.is_gamma0_equiv(z, 3, True)
            (True, 2)

        ALGORITHM:
            See Proposition 2.2.3 of Cremona's book "Algorithms for Modular
            Elliptic Curves", or Prop 2.27 of Stein's Ph.D. thesis.
        """
        N = ZZ(N)
        u1 = self.__a
        v1 = self.__b
        u2 = other.__a
        v2 = other.__b
        if (u1,v1) != (ZZ(0),ZZ(1)):
            if v1 in [ZZ(0),ZZ(1)]:
                s1 = ZZ(1)
            else:
                s1 = u1.inverse_mod(v1)
        else:
            s1 = 0
        if (u2,v2) != (ZZ(0),ZZ(1)):
            if v2 in [ZZ(0),ZZ(1)]:
                s2 = 1
            else:
                s2 = u2.inverse_mod(v2)
        else:
            s2 = 0
        g = (v1*v2).gcd(N)
        a = s1*v2 - s2*v1
        if a%g != ZZ(0):
            if transformation:
                return False, None
            else:
                return False

        if not transformation:
            return True

        # Now we know the cusps are equivalent.  Use the proof of Prop 2.2.3
        # of Cremona to find a matrix in Gamma_0(N) relating them.
        dum, s2, r2 = u2.xgcd(-v2)
        assert dum.is_one()
        dum, s1, r1 = u1.xgcd(-v1)
        assert dum.is_one()
        a = s1*v2 - s2*v1
        assert (a%g).is_zero()
        # solve x*v1*v2 + a = 0 (mod N).
        d,x0,y0 = (v1*v2).xgcd(N)          # x0*v1*v2 + y0*N = d = g.
        # so x0*v1*v2 - g = 0 (mod N)
        x = -x0 * (a/g)
        # now  x*v1*v2 + a = 0 (mod N)
        s1p = s1+x*v1
        return (True, (u2*s1p-r2*v1)%N)

    def is_gamma1_equiv(self, other, N):
        """
        Return whether self and other are equivalent modulo the action
        of Gamma_1(N) via linear fractional transformations.

        INPUT:
            other -- Cusp
            N -- an integer (specifies the group Gamma_1(N))

        OUTPUT:
            bool -- True if self and other are equivalent
            int -- 0, 1 or -1, gives further information
                   about the equivalence:  If the two cusps
                   are u1/v1 and u2/v2, then they are equivalent
                   if and only if
                        v1 = v2 (mod N) and u1 = u2 (mod gcd(v1,N))
                   or
                        v1 = -v2 (mod N) and u1 = -u2 (mod gcd(v1,N))
                   The sign is +1 for the first and -1 for the second.
                   If the two cusps are not equivalent then 0 is returned.
        """
        N = ZZ(N)
        u1 = self.__a
        v1 = self.__b
        u2 = other.__a
        v2 = other.__b
        g = v1.gcd(N)
        if ((v2 - v1) % N == 0 and (u2 - u1)%g== 0):
            return True, 1
        elif ((v2 + v1) % N == 0 and (u2 + u1)%g== 0):
            return True, -1
        return False, 0

    def apply(self, g):
        """
        Return g(self), where g=[a,b,c,d] is a list of length 4, which
        we view as a linear fractional transformation.
        """
        return Cusp(g[0]*self.__a + g[1]*self.__b, g[2]*self.__a + g[3]*self.__b)


