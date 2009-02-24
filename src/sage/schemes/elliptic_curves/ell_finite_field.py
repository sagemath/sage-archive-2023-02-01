"""
Elliptic curves over finite fields

AUTHORS:

- William Stein (2005): Initial version

- Robert Bradshaw et al....

- John Cremona (2008-02): Point counting and group structure for
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

from sage.misc.randstate import current_randstate
import sys
from math import ceil, floor, sqrt

from ell_generic import Hasse_bounds
from ell_field import EllipticCurve_field
from constructor import EllipticCurve
from sage.schemes.hyperelliptic_curves.hyperelliptic_finite_field import HyperellipticCurve_finite_field
import sage.rings.ring as ring
from sage.rings.all import Integer, ZZ, PolynomialRing, ComplexField, FiniteField, GF, polygen
import gp_cremona
import sea
from sage.groups.all import AbelianGroup
import sage.groups.generic as generic
import ell_point
from sage.calculus.calculus import log
from sage.rings.arith import integer_ceil, integer_floor, gcd
from sage.structure.sequence import Sequence

import sage.plot.all as plot

import sage.libs.pari
pari = sage.libs.pari.all.pari

class EllipticCurve_finite_field(EllipticCurve_field, HyperellipticCurve_finite_field):
    """
    Elliptic curve over a finite field.
    """
    def __init__(self, x, y=None):
        """
        Special constructor for elliptic curves over a finite field

        EXAMPLES::

            sage: EllipticCurve(GF(101),[2,3])
            Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field of size 101

        ::

            sage: F=GF(101^2, 'a')
            sage: EllipticCurve([F(2),F(3)])
            Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field in a of size 101^2
        """
        if isinstance(x, list):
            seq = Sequence(x)
        else:
            seq = Sequence(y, universe=x)
        ainvs = list(seq)
        field = seq.universe()
        if not isinstance(field, ring.Ring):
            raise TypeError

        EllipticCurve_field.__init__(self, ainvs)

        self._point_class = ell_point.EllipticCurvePoint_finite_field

    def _pari_(self):
        """
        Return a GP/PARI elliptic curve

        EXAMPLES::

            sage: EllipticCurve(GF(41),[2,5])._pari_()
            [Mod(0, 41), Mod(0, 41), Mod(0, 41), Mod(2, 41), Mod(5, 41), Mod(0, 41), Mod(4, 41), Mod(20, 41), Mod(37, 41), Mod(27, 41), Mod(26, 41), Mod(4, 41), Mod(11, 41), 0, 0, 0, 0, 0, 0]
        """
        try:
            return self.__pari
        except AttributeError:
            pass
        F = self.base_ring()
        self.__pari = pari('ellinit(Mod(1,%s)*%s)'%(F.characteristic(), [b._pari_() for b in self.ainvs()]))
        return self.__pari

    def _gp(self):
        """
        Return an elliptic curve in a GP/PARI interpreter with all
        Cremona's code for finite fields preloaded. This includes
        generators, which will vary from run to run.

        The base field must have prime order.

        EXAMPLES::

            sage: EllipticCurve(GF(41),[2,5])._gp()
            [Mod(0, 41), Mod(0, 41), Mod(0, 41), Mod(2, 41), Mod(5, 41), Mod(0, 41), Mod(4, 41), Mod(20, 41), Mod(37, 41), Mod(27, 41), Mod(26, 41), Mod(4, 41), Mod(11, 41), 44, [2, 2; 11, 1], [22, 2], ...
        """
        try:
            return self.__gp
        except AttributeError:
            pass
        F = self.base_ring()
        if not F.is_prime_field():
            raise NotImplementedError
        self.__gp = gp_cremona.ellinit(self.a_invariants(), F.characteristic())
        return self.__gp

    def plot(self, *args, **kwds):
        """
        Draw a graph of this elliptic curve over a prime finite field.

        INPUT:


        -  ``*args, **kwds`` - all other options are passed
           to the circle graphing primitive.


        EXAMPLES::

            sage: E = EllipticCurve(FiniteField(17), [0,1])
            sage: P = plot(E, rgbcolor=(0,0,1))
        """
        R = self.base_ring()
        if not R.is_prime_field():
            raise NotImplementedError

        G = plot.Graphics()
        G += plot.points([P[0:2] for P in self.points() if not P.is_zero()], *args, **kwds)

        return G

    def _points_via_group_structure(self):
        """
        Return a list of all the points on the curve, for prime fields only
        (see points() for the general case)

        EXAMPLES::

            sage: S=EllipticCurve(GF(97),[2,3])._points_via_group_structure()
            sage: len(S)
            100

        See trac #4687, where the following example did not work::

            sage: E=EllipticCurve(GF(2),[0, 0, 1, 1, 1])
            sage: E.points()
            [(0 : 1 : 0)]

        ::

            sage: E=EllipticCurve(GF(2),[0, 0, 1, 0, 1])
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (1 : 1 : 1)]

        ::

            sage: E=EllipticCurve(GF(4,'a'),[0, 0, 1, 0, 1])
            sage: E.points()
            [(0 : 1 : 0), (0 : a : 1), (0 : a + 1 : 1), (1 : 0 : 1), (1 : 1 : 1), (a : 0 : 1), (a : 1 : 1), (a + 1 : 0 : 1), (a + 1 : 1 : 1)]
        """
        # TODO, eliminate when polynomial calling is fast
        G, pts = self.abelian_group()

        ni = G.invariants()
        ngens = G.ngens()

        H0=[self(0)]
        if ngens == 0:    # trivial group
            return H0
        for m in range(1,ni[0]): H0.append(H0[-1]+pts[0])
        if ngens == 1:    # cyclic group
            return H0

        # else noncyclic group
        H1=[self(0)]
        for m in range(1,ni[1]): H1.append(H1[-1]+pts[1])
        return [P+Q for P in H0 for Q in H1]

    def points(self):
        r"""
        All the points on this elliptic curve. The list of points is cached
        so subsequent calls are free.

        EXAMPLES::

            sage: p = 5
            sage: F = GF(p)
            sage: E = EllipticCurve(F, [1, 3])
            sage: a_sub_p = E.change_ring(QQ).ap(p); a_sub_p
            2

        ::

            sage: len(E.points())
            4
            sage: p + 1 - a_sub_p
            4
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (4 : 1 : 1), (4 : 4 : 1)]

        ::

            sage: K = GF(p**2,'a')
            sage: E = E.change_ring(K)
            sage: len(E.points())
            32
            sage: (p + 1)**2 - a_sub_p**2
            32
            sage: w = E.points(); w
            [(0 : 1 : 0), (0 : 2*a + 4 : 1), (0 : 3*a + 1 : 1), (1 : 0 : 1), (2 : 2*a + 4 : 1), (2 : 3*a + 1 : 1), (3 : 2*a + 4 : 1), (3 : 3*a + 1 : 1), (4 : 1 : 1), (4 : 4 : 1), (a : 1 : 1), (a : 4 : 1), (a + 2 : a + 1 : 1), (a + 2 : 4*a + 4 : 1), (a + 3 : a : 1), (a + 3 : 4*a : 1), (a + 4 : 0 : 1), (2*a : 2*a : 1), (2*a : 3*a : 1), (2*a + 4 : a + 1 : 1), (2*a + 4 : 4*a + 4 : 1), (3*a + 1 : a + 3 : 1), (3*a + 1 : 4*a + 2 : 1), (3*a + 2 : 2*a + 3 : 1), (3*a + 2 : 3*a + 2 : 1), (4*a : 0 : 1), (4*a + 1 : 1 : 1), (4*a + 1 : 4 : 1), (4*a + 3 : a + 3 : 1), (4*a + 3 : 4*a + 2 : 1), (4*a + 4 : a + 4 : 1), (4*a + 4 : 4*a + 1 : 1)]

        Note that the returned list is an immutable sorted Sequence::

            sage: w[0] = 9
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        try:
            return self.__points
        except AttributeError: pass

        from sage.structure.sequence import Sequence
        if self.base_ring().is_prime_field():
            v = self._points_via_group_structure()
        else:
            v =self._points_fast_sqrt()
        v.sort()
        self.__points = Sequence(v, immutable=True)
        return self.__points

    def random_element(self):
        """
        Returns a random point on this elliptic curve.

        Returns the point at infinity with probability `1/(q+1)`
        where the base field has cardinality `q`.

        EXAMPLES::

            sage: k = GF(next_prime(7^5))
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element(); P
            (751 : 6230 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        ::

            sage: k.<a> = GF(7^5)
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element(); P
            (a^4 + a + 5 : 6*a^4 + 3*a^3 + 2*a^2 + 4 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        ::

            sage: k.<a> = GF(2^5)
            sage: E = EllipticCurve(k,[a^2,a,1,a+1,1])
            sage: P = E.random_element(); P
            (a^4 : 0 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True
        """
        random = current_randstate().python_random().random
        k = self.base_field()
        # The following allows the origin self(0) to be picked
        if random() <= 1/float(k.order()+1):
            return self(0)

        while True:
            try:
                return self.lift_x(k.random_element())
            except:
                pass

    random_point = random_element


    def trace_of_frobenius(self):
        r"""
        Return the trace of Frobenius acting on this elliptic curve.

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        EXAMPLES::

            sage: E=EllipticCurve(GF(101),[2,3])
            sage: E.trace_of_frobenius()
            6
            sage: E=EllipticCurve(GF(11^5,'a'),[2,5])
            sage: E.trace_of_frobenius()
            802

        The following shows that the issue from trac #2849 is fixed::

            sage: E=EllipticCurve(GF(3^5,'a'),[-1,-1])
            sage: E.trace_of_frobenius()
            -27
        """
        return 1 + self.base_field().order() - self.cardinality()

    def _cardinality_with_j_invariant_1728(self):
        r"""
        Special function to compute cardinality when j=1728.

        EXAMPLES: An example with q=p=1 (mod 4)

        ::

            sage: F=GF(10009)
            sage: [EllipticCurve(F,[0,0,0,11^i,0])._cardinality_with_j_invariant_1728() for i in range(4)]
            [10016, 10210, 10004, 9810]

        An example with q=p=3 (mod 4)

        ::

            sage: F=GF(10007)
            sage: [EllipticCurve(F,[0,0,0,5^i,0])._cardinality_with_j_invariant_1728() for i in range(4)]
            [10008, 10008, 10008, 10008]

        An example with `q=p^2`, p=3 (mod 4)

        ::

            sage: F.<a>=GF(10007^2,'a')
            sage: [EllipticCurve(F,[0,0,0,a^i,0])._cardinality_with_j_invariant_1728() for i in range(4)]
            [100160064, 100140050, 100120036, 100140050]

        Examples with `q=2^d`, d odd (3 isomorphism classes)::

            sage: F.<a> = GF(2**15,'a')
            sage: ais = [[0,0,1,0,0],[0,0,1,1,0],[0,0,1,1,1]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [e._cardinality_with_j_invariant_1728() for e in curves]
            [32769, 33025, 32513]

        Examples with `q=2^d`, d even (7 isomorphism classes)::

            sage: F.<a> = GF(2**16,'a')
            sage: b = a^11 # trace 1
            sage: ais = [[0,0,1,0,0],[0,0,1,0,b],[0,0,1,b,0],[0,0,a,0,0],[0,0,a,0,a^2*b],[0,0,a^2,0,0],[0,0,a^2,0,a^4*b]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [e._cardinality_with_j_invariant_1728() for e in curves]
            [65025, 66049, 65537, 65793, 65281, 65793, 65281]

        Examples with `q=3^d`, d odd (4 isomorphism classes)::

            sage: F.<a> = GF(3**15,'a')
            sage: b=a^7  # has trace 1
            sage: ais=[[0,0,0,1,0],[0,0,0,-1,0],[0,0,0,-1,b],[0,0,0,-1,-b]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [e._cardinality_with_j_invariant_1728() for e in curves]
            [14348908, 14348908, 14342347, 14355469]

        Examples with `q=3^d`, d even (6 isomorphism classes)::

            sage: F.<g>=GF(3^18,'g')
            sage: i=F(-1).sqrt()
            sage: a=g^8  # has trace 1
            sage: ais= [[0,0,0,1,0],[0,0,0,1,i*a],[0,0,0,g,0],[0,0,0,g^3,0],[0,0,0,g^2,0], [0,0,0,g^2,i*a*g^3]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [E._cardinality_with_j_invariant_1728() for E in curves]
            [387459856, 387400807, 387420490, 387420490, 387381124, 387440173]
        """
        try:
            return self._order
        except AttributeError:
            pass

        k = self.base_ring()
        assert self.j_invariant()==k(1728)
        q = k.cardinality()
        p = k.characteristic()
        d = k.degree()
        x=polygen(ZZ)

# p=2, j=0=1728
#
# Number of isomorphism classes is 3 (in odd degree) or 7 (in even degree)
#
        if p==2:
            if d%2==1:
                # The 3 classes are represented, independently of d,
                # by [0,0,1,0,0], [0,0,1,1,0], [0,0,1,1,1]
                E=EllipticCurve(k,[0,0,1,0,0])
                if self.is_isomorphic(E):
                    N = q+1
                else:
                    n = (d+1)//2
                    t = 2**n
                    n = n%4
                    if n==0 or n==1: t=-t
                    E=EllipticCurve(k,[0,0,1,1,1])
                    if self.is_isomorphic(E): t=-t
                    N = q+1-t
            else:
                # The 7 classes are represented by E1=[0,0,1,0,0],
                # E2=[0,0,1,0,b], E3=[0,0,1,b,0], E4=[0,0,a,0,0],
                # E4=[0,0,a,0,a^2*b], E6=[0,0,a^2,0,0],
                # E7=[0,0,a^2,0,a^4*b], where a is a non-cube and b
                # has trace 1.  E1's Frobenius is pi=(-2)**(d//2); the
                # Frobeniuses are then pi, -pi, 0; w*pi, -w*pi;
                # w^2*pi, -w^2*pi where w is either cube root of
                # unity, so the traces are 2*pi, -2*pi, 0, -pi, +pi;
                # -pi, +pi.
                delta = self.discriminant()
                discube = (delta**((q-1)//3) == k(1))
                pi = (-2)**(d//2)
                if discube:
                    a = k.gen()
                    b = a
                    while b.trace()==0: b*=a
                    if self.is_isomorphic(EllipticCurve(k,[0,0,1,b,0])):
                        t = 0
                    else:
                        t = 2*pi
                        if not self.is_isomorphic(EllipticCurve(k,[0,0,1,0,0])):
                            t = -t

                else:
                    t = pi
                    if self.is_isomorphic(EllipticCurve(k,[0,0,delta,0,0])):
                        t = -t
                N = q+1-t


# p=3, j=0=1728
#
# Number of isomorphism classes is 4 (odd degree) or 6 (even degree)
#
        elif p==3:
            if d%2==1:
                # The 4 classes are represented by [0,0,0,1,0],
                # [0,0,0,-1,0], [0,0,0,-1,a], [0,0,0,-1,-a] where a
                # has trace 1
                delta = self.discriminant()
                if (-delta).is_square():
                    t = 0
                else:
                    u = delta.sqrt()
                    if not u.is_square(): u=-u
                    tr = ((self.a3()**2+self.a6())/u).trace()
                    if tr==0:
                        t = 0
                    else:
                        d2 = (d+1)//2
                        t = 3**d2
                        if d2%2==1: t = -t
                        if tr==-1:  t = -t
                N = q+1-t
            else:
                # The 6 classes are represented by: [0,0,0,1,0],
                # [0,0,0,1,i*a]; [0,0,0,g,0], [0,0,0,g^3,0];
                # [0,0,0,g^2,0], [0,0,0,g^2,i*a*g^3]; where g
                # generates the multiplicative group modulo 4th
                # powers, and a has nonzero trace.
                A4 = self.a4()-self.a1()*self.a3()
                i = k(-1).sqrt()
                t = 0
                if A4.is_square():
                    u = A4.sqrt()
                    t = (-3)**(d//2)
                    if u.is_square():
                        A6 = (self.a3()**2+self.a6())/u
                        if (i*A6).trace()==0:
                            t = 2*t
                        else:
                            t = -t
                    else:
                        A6 = (self.a3()**2+self.a6())/(u*A4)
                        if (i*A6).trace()==0:
                            t = -2*t
                N = q+1-t

# p>3, j=1728
#
# Number of isomorphism classes is 4 if q=1 (mod 4), else 2
#
        elif p%4==3:
            if d%2==1:
                t = 0
            else:
                t  = (-p)**(d//2)
                w = (self.c4()/k(48))**((q-1)//4)
                if   w==1:    t =  2*t
                elif w==-1:   t = -2*t
                else: t = 0

            N = q+1-t

# p=1 (mod 4).  First find Frobenius pi=a+b*i for [0,0,0,-1,0] over GF(p):
# N(pi)=p and N(pi-1)=0 (mod 8).
#
        else: # p%4==1
            R = ZZ.extension(x**2+1,'i')
            i = R.gen(1)
            pi = R.fraction_field().factor(p)[0][0].gens_reduced()[0]
            a,b = pi.list()
            if a%2==0:
                a,b = -b,a
            if (a+b+1)%4==0:
                a,b = -a,-b
            pi = a+b*i        # Now pi=a+b*i with (a,b)=(1,0),(3,2) mod 4

        # Lift to Frobenius for [0,0,0,-1,0] over GF(p^d):

            if d>1:
                pi = pi**d
                a,b = pi.list()

        # Compute appropriate quartic twist:

            w = (self.c4()/k(48))**((q-1)//4)
            if w==1:
                t = 2*a
            elif w==-1:
                t = -2*a
            elif k(b)==w*k(a):
                t = 2*b
            else:
                t = -2*b
            N = q+1-t

        self._order = Integer(N)
        return self._order

    def _cardinality_with_j_invariant_0(self):
        r"""
        Special function to compute cardinality when j=0.

        EXAMPLES: An example with q=p=1 (mod 6)

        ::

            sage: F=GF(1009)
            sage: [EllipticCurve(F,[0,0,0,0,11^i])._cardinality_with_j_invariant_0() for i in range(6)]
            [948, 967, 1029, 1072, 1053, 991]

        An example with q=p=5 (mod 6)

        ::

            sage: F=GF(1013)
            sage: [EllipticCurve(F,[0,0,0,0,3^i])._cardinality_with_j_invariant_0() for i in range(6)]
            [1014, 1014, 1014, 1014, 1014, 1014]

        An example with `q=p^2`, p=5 (mod 6)

        ::

            sage: F.<a>=GF(1013^2,'a')
            sage: [EllipticCurve(F,[0,0,0,0,a^i])._cardinality_with_j_invariant_0() for i in range(6)]
            [1028196, 1027183, 1025157, 1024144, 1025157, 1027183]

        For examples in characteristic 2 and 3, see the function
        _cardinality_with_j_invariant_1728()
        """

        try:
            return self._order
        except AttributeError:
            pass

        k = self.base_ring()
        assert self.j_invariant()==k(0)
        p = k.characteristic()
        if p==2 or p==3:  # then 0==1728
            return self._cardinality_with_j_invariant_1728()

        q = k.cardinality()
        d = k.degree()
        x=polygen(ZZ)

# p>3, j=0
#
# Number of isomorphism classes is 4 if q=1 (mod 4), else 2
#
        if p%6==5:
            if d%2==1:
                t = 0
            else:
                t  = (-p)**(d//2)
                w = (self.c6()/k(-864))**((q-1)//6)
                if   w==1:    t =  2*t
                elif w==-1:   t = -2*t
                elif w**3==1: t = -t

            N = q+1-t

# p=1 (mod 6).  First find Frobenius pi=a+b*w for [0,0,0,0,1] over GF(p):
# N(pi)=p and N(pi-1)=0 (mod 12).
#
        else: # p%6==1
            R = ZZ.extension(x**2-x+1,'zeta6')
            zeta6 = R.gen(1)
            pi = R.fraction_field().factor(p)[0][0].gens_reduced()[0]
            while (pi-1).norm()%12 !=0:  pi*=zeta6
            a,b = pi.list()
            z = k(-b)/k(a)  # a *specific* 6th root of unity in k

# Now pi=a+b*zeta6 with N(pi-1)=0 (mod 12)

# Lift to Frobenius for [0,0,0,0,1] over GF(p^d):

            if d>1:
                pi = pi**d
                a,b = pi.list()

# Compute appropriate sextic twist:

            w = (self.c6()/k(-864))**((q-1)//6)

            if   w==1:       t = 2*a+b  # = Trace(pi)
            elif w==-1:      t = -2*a-b # = Trace(-pi)
            elif w==z:       t = a-b    # = Trace(pi*zeta6)
            elif w==z**2:    t = -a-2*b # = Trace(pi*zeta6**2)
            elif w==z**4:    t = b-a    # = Trace(pi*zeta6**4)
            elif w==z**5:    t = a+2*b  # = Trace(pi*zeta6**5)

            N = q+1-t

        self._order = Integer(N)
        return self._order

    def cardinality(self, algorithm='heuristic', extension_degree=1):
        r"""
        Return the number of points on this elliptic curve over an
        extension field (default: the base field).

        INPUT:


        -  ``algorithm`` - string (default: 'heuristic'), used
           only for point counting over prime fields

            -  ``'heuristic'`` - use a heuristic to choose between
               pari, bsgs and sea.

            -  ``'pari'`` - use the baby step giant step method as
               implemented in PARI via the C-library function ellap.

            -  ``'sea'`` - use sea.gp as implemented in PARI by
               Christophe Doche and Sylvain Duquesne.

            -  ``bsgs`` - use the baby step giant step method as
               implemented in Sage, with the Cremona -
               Sutherland version of Mestre's trick.

            - ``all`` - (over prime fields only) compute
              cardinality with all of pari, sea and
              bsgs; return result if they agree or raise
              a RuntimeError if they do not.

        -  ``early_abort`` - bool (default: False); this is
           used only by sea. if True, stop early if a small factor of the
           order is found.

        -  ``extension_degree`` - int (default: 1); if the
           base field is `k=GF(p^n)` and extension_degree=d, returns
           the cardinality of `E(GF(p^{n d}))`.


        OUTPUT: an integer

        The cardinality is cached.

        Over prime fields, one of the above algorithms is used. Over
        non-prime fields, the serious point counting is done on a standard
        curve with the same j-invariant over the field GF(p)(j), then
        lifted to the base_field, and finally account is taken of twists.

        For j=0 and j=1728 special formulas are used instead.

        EXAMPLES::

            sage: EllipticCurve(GF(4,'a'),[1,2,3,4,5]).cardinality()
            8
            sage: k.<a> = GF(3^3)
            sage: l = [a^2 + 1, 2*a^2 + 2*a + 1, a^2 + a + 1, 2, 2*a]
            sage: EllipticCurve(k,l).cardinality()
            29

        ::

            sage: l = [1, 1, 0, 2, 0]
            sage: EllipticCurve(k,l).cardinality()
            38

        An even bigger extension (which we check against Magma)::

            sage: EllipticCurve(GF(3^100,'a'),[1,2,3,4,5]).cardinality()
            515377520732011331036459693969645888996929981504
            sage: magma.eval("Order(EllipticCurve([GF(3^100)|1,2,3,4,5]))")    # optional - magma
            '515377520732011331036459693969645888996929981504'

        ::

            sage: EllipticCurve(GF(10007),[1,2,3,4,5]).cardinality()
            10076
            sage: EllipticCurve(GF(10007),[1,2,3,4,5]).cardinality(algorithm='sea')
            10076
            sage: EllipticCurve(GF(10007),[1,2,3,4,5]).cardinality(algorithm='pari')
            10076
            sage: EllipticCurve(GF(next_prime(10**20)),[1,2,3,4,5]).cardinality(algorithm='sea')
            100000000011093199520

        The cardinality is cached::

            sage: E = EllipticCurve(GF(3^100,'a'),[1,2,3,4,5])
            sage: E.cardinality() is E.cardinality()
            True
            sage: E=EllipticCurve(GF(11^2,'a'),[3,3])
            sage: E.cardinality()
            128
            sage: EllipticCurve(GF(11^100,'a'),[3,3]).cardinality()
            137806123398222701841183371720896367762643312000384671846835266941791510341065565176497846502742959856128
        """
        if extension_degree>1:
            # A recursive call to cardinality() with
            # extension_degree=1, which will cache the cardinality, is
            # made by the call to frobenius_order() here:
            R=self.frobenius_order()
            if R.degree()==1:
                return (self.frobenius()**extension_degree-1)**2
            else:
                return (self.frobenius()**extension_degree-1).norm()

        # Now extension_degree==1
        if algorithm != 'all':
            try:
                return self._order
            except AttributeError:
                pass

        k = self.base_ring()

        # use special code for j=0, 1728 (for any field)
        j = self.j_invariant()
        if j==k(0):
            return self._cardinality_with_j_invariant_0()
        if j==k(1728):
            return self._cardinality_with_j_invariant_1728()

        N = 0
        q = k.cardinality()
        p = k.characteristic()
        d = k.degree()

        # Over prime fields, we have a variety of algorithms to choose from:

        if d == 1:
            if algorithm == 'heuristic':
                algorithm = 'sea' if p > 10**18 else 'pari'
            if algorithm == 'pari':
                N = self.cardinality_pari()
            elif algorithm == 'sea':
                N = self.cardinality_sea()
            elif algorithm == 'bsgs':
                N = self.cardinality_bsgs()
            elif algorithm == 'all':
                N1 = self.cardinality_pari()
                N2 = self.cardinality_sea()
                N3 = self.cardinality_bsgs()
                if N1 == N2 and N2 == N3:
                    N = N1
                else:
                    if N1!=N2:
                        raise RuntimeError, "BUG! Cardinality with pari=%s but with sea=%s"%(N1, N2)
                    if N1!=N3:
                        raise RuntimeError, "BUG! Cardinality with pari=%s but with bsgs=%s"%(N1, N3)
            self._order = Integer(N)
            return self._order

        # now k is not a prime field and j is not 0, 1728

        # we count points on a standard curve with the same
        # j-invariant, defined over the field it generates, then lift
        # to the curve's own field, and finally allow for twists

        # Since j is not 0, 1728 the only twists are quadratic

        j_pol=j.minimal_polynomial()
        j_deg=j_pol.degree()

        # if not possible to work over a smaller field:
        if d==j_deg:
            self._order = self.cardinality_bsgs()
            return self._order

        kj=GF(p**j_deg,name='a',modulus=j_pol)
        jkj=kj.gen() if j_deg>1 else j_pol.roots(multiplicities=False)[0]

        # recursive call which will do all the real work:
        Ej = EllipticCurve(jkj)
        N=Ej.cardinality(extension_degree=d//j_deg)

        # if curve ia a (quadratic) twist of the "standard" one:
        if not self.is_isomorphic(EllipticCurve(j)): N=2*(q+1)-N

        self._order = N
        return self._order

    order = cardinality # alias

    def frobenius_polynomial(self):
        r"""
        Return the characteristic polynomial of Frobenius.

        The Frobenius endomorphism of the elliptic curve has quadratic
        characteristic polynomial. In most cases this is irreducible and
        defines an imaginary quadratic order; for some supersingular
        curves, Frobenius is an integer a and the polynomial is
        `(x-a)^2`.

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_polynomial()
            x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z and the polynomial
        is a square::

            sage: E=EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: E.frobenius_polynomial().factor()
            (x + 5)^2
        """
        x=polygen(ZZ)
        return x**2-self.trace_of_frobenius()*x+self.base_field().cardinality()

    def frobenius_order(self):
        r"""
        Return the quadratic order Z[phi] where phi is the Frobenius
        endomorphism of the elliptic curve

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_order()
            Order in Number Field in phi with defining polynomial x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z and the Frobenius
        order is Z::

            sage: E=EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: R=E.frobenius_order()
            sage: R
            Order in Number Field in phi with defining polynomial x + 5
            sage: R.degree()
            1
        """
        f = self.frobenius_polynomial().factor()[0][0]
        return ZZ.extension(f,names='phi')

    def frobenius(self):
        r"""
        Return the frobenius of self as an element of a quadratic order

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        Frobenius is only determined up to conjugacy.

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius()
            phi
            sage: E.frobenius().minpoly()
            x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z::

            sage: E=EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: E.frobenius()
            -5
        """
        R = self.frobenius_order()
        if R.degree()==1:
            return self.frobenius_polynomial().roots(multiplicities=False)[0]
        else:
            return R.gen(1)

    def cardinality_exhaustive(self):
        r"""
        Return the cardinality of self over the base field. Simply adds up
        the number of points with each x-coordinate: only used for small
        field sizes!

        EXAMPLES::

            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_exhaustive()
            1020
            sage: E=EllipticCurve(GF(3^4,'a'),[1,1])
            sage: E.cardinality_exhaustive()
            64
        """
        self._order = Integer(1+sum([len(self.lift_x(x,all=True)) for x in self.base_field()]))
        return self._order

    def cardinality_pari(self):
        r"""
        Return the cardinality of self over the (prime) base field using pari.

        The result is not cached.

        EXAMPLES::

            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_pari()
            1020
            sage: K=GF(next_prime(10^6))
            sage: E=EllipticCurve(K,[1,0,0,1,1])
            sage: E.cardinality_pari()
            999945

        TESTS::

            sage: K.<a>=GF(3^20)
            sage: E=EllipticCurve(K,[1,0,0,1,a])
            sage: E.cardinality_pari()
            Traceback (most recent call last):
            ...
            ValueError: cardinality_pari() only works over prime fields.
            sage: E.cardinality()
            3486794310

        """
        k = self.base_ring()
        p = k.characteristic()
        if k.degree()==1:
            return ZZ(p + 1 - int(self._pari_().ellap(p)))
        else:
            raise ValueError, "cardinality_pari() only works over prime fields."

    def cardinality_sea(self, early_abort=False):
        r"""
        Return the cardinality of self over the (prime) base field using sea.

        INPUT:

        - ``early_abort`` - bool (default: False).  if True, an early
          abort technique is used and the computation is interrupted
          as soon as a small divisor of the order is detected.  The
          function then returns 0.  This is useful for ruling out
          curves whose cardinality is divisible by a small prime.

        The result is not cached.

        EXAMPLES::

            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_sea()
            1020
            sage: K=GF(next_prime(10^6))
            sage: E=EllipticCurve(K,[1,0,0,1,1])
            sage: E.cardinality_sea()
            999945

        TESTS::

            sage: K.<a>=GF(3^20)
            sage: E=EllipticCurve(K,[1,0,0,1,a])
            sage: E.cardinality_sea()
            Traceback (most recent call last):
            ...
            ValueError: cardinality_sea() only works over prime fields.
            sage: E.cardinality()
            3486794310

        """
        k = self.base_ring()
        p = k.characteristic()
        if k.degree()==1:
            return sea.ellsea(self.a_invariants(), p, early_abort=early_abort)
        else:
            raise ValueError, "cardinality_sea() only works over prime fields."

    def cardinality_bsgs(self, verbose=False):
        r"""
        Return the cardinality of self over the base field. Will be called
        by user function cardinality only when necessary, i.e. when the
        j_invariant is not in the prime field.

        ALGORITHM: A variant of "Mestre's trick" extended to all finite
        fields by Cremona and Sutherland, 2008.

        .. note::

           1. The Mestre-Schoof-Cremona-Sutherland algorithm may fail for
              a small finite number of curves over `F_q` for `q` at most 49, so
              for `q<50` we use an exhaustive count.

           2. Quadratic twists are not implemented in characteristic 2
              when `j=0 (=1728)`; but this case is treated separately.

       EXAMPLES::

            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_bsgs()
            1020
            sage: E=EllipticCurve(GF(3^4,'a'),[1,1])
            sage: E.cardinality_bsgs()
            64
            sage: F.<a>=GF(101^3,'a')
            sage: E=EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
            sage: E.cardinality_bsgs()
            1031352
        """
        E1 = self
        k = self.base_field()
        q = k.order()
        if q<50:
            if verbose:
                print "q=",q,"< 50 so using exhaustive count"
            return self.cardinality_exhaustive()

        # Construct the quadratic twist:
        E2 = E1.quadratic_twist()
        if verbose:
            print "Quadratic twist is ",E2.ainvs()

        bounds = Hasse_bounds(q)
        lower, upper = bounds
        B = upper-q-1 # = floor(2*sqrt(q))
        a = ZZ(0)
        N1 = N2 = M = ZZ(1)
        kmin = -B
        kmax = B
        q1 = q+1
        # Throughout, we have #E=q+1-t where |t|<=B and t=a+k*M = a
        # (mod M) where kmin <= k <= kmax.

        # M is the lcm of the orders of all the points found on E1 and
        # E2, which will eventually exceed 2*B, at which point
        # kmin=kmax.

        if q > 2**10:
            N1 = ZZ(2)**sum([e for P,e in E1._p_primary_torsion_basis(2)])
            N2 = ZZ(2)**sum([e for P,e in E2._p_primary_torsion_basis(2)])
            if q > 2**20:
                N1 *= ZZ(3)**sum([e for P,e in E1._p_primary_torsion_basis(3)])
                N2 *= ZZ(3)**sum([e for P,e in E2._p_primary_torsion_basis(3)])
                if q > 2**40:
                    N1 *= ZZ(5)**sum([e for P,e in E1._p_primary_torsion_basis(5)])
                    N2 *= ZZ(5)**sum([e for P,e in E2._p_primary_torsion_basis(5)])
            # We now know that t=q+1 (mod N1) and t=-(q+1) (mod N2)
            a = q1
            M = N1
            g,u,v = M.xgcd(N2) # g==u*M+v*N2
            if N2>g:
                a = (a*v*N2-q1*u*M)//g
                M *= (N2//g) # = lcm(M,N2)
                a = a%M
                if verbose:
                    print "(a,M)=",(a,M)
                kmin = ((-B-a)/M).ceil()
                kmax = ((B-a)/M).floor()
                if kmin==kmax:
                    self._order = q1-a-kmin*M
                    if verbose: print "no random points were needed"
                    return self._order
            if verbose: print "(2,3,5)-torsion subgroup gives M=",M

        # N1, N2 are divisors of the orders of E1, E2 separately,
        # which are used to speed up the computation of the orders of
        # points on each curve.  For large q it is worth initializing
        # these with the full order of the (2,3,5)-torsion which are
        # often non-trivial.

        while kmax!=kmin:
            # Get a random point on E1 and find its order, using the
            # Hasse bounds and the fact that we know that the group
            # order is a multiple of N1:
            n = generic.order_from_bounds(E1.random_point(),bounds,N1,operation='+')
            if verbose: print "New point on E has order ",n
            # update N1 and M
            N1 = N1.lcm(n)
            g,u,v = M.xgcd(n) # g==u*M+v*n
            if n>g:
                # update congruence a (mod M) with q+1 (mod n)
                a = (a*v*n+q1*u*M)//g
                M *= (n//g) # = lcm(M,n)
                a = a%M
                if verbose: print "(a,M)=",(a,M)
                kmin = ((-B-a)/M).ceil()
                kmax = ((B-a)/M).floor()
                if kmin==kmax:
                    self._order = q1-a-kmin*M
                    return self._order
                if verbose: print "number of possibilities is now ",kmax-kmin+1

            # Get a random point on E2 and find its order, using the
            # Hasse bounds and the fact that we know that the group
            # order is a multiple of N2:
            n = generic.order_from_bounds(E2.random_point(),bounds,N2,operation='+')
            if verbose:  print "New point on E' has order ",n
            # update N2 and M
            N2 = N2.lcm(n)
            g,u,v = M.xgcd(n) # g==u*M+v*n
            if n>g:
                # update congruence a (mod M) with -(q+1) (mod n)
                a = (a*v*n-q1*u*M)//g
                M *= (n//g) # = lcm(M,n)
                a = a%M
                if verbose: print "(a,M)=",(a,M)
                kmin = ((-B-a)/M).ceil()
                kmax = ((B-a)/M).floor()
                if kmin==kmax:
                    self._order = q1-a-kmin*M
                    return self._order
                if verbose: print "number of possibilities is now ",kmax-kmin+1

    def gens(self):
        """
        Returns a tuple of length up to 2 of points which generate the
        abelian group of points on this elliptic curve. See
        abelian_group() for limitations.

        The algorithm uses random points on the curve, and hence the
        generators are likely to differ from one run to another; but they
        are cached so will be consistent in any one run of Sage.

        AUTHORS:

        - John Cremona

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[2,5]) # random output
            sage: E.gens()
            ((0 : 4 : 1),)
            sage: EllipticCurve(GF(41),[2,5]).gens() # random output
            ((21 : 1 : 1), (8 : 0 : 1))
            sage: F.<a>=GF(3^6,'a')
            sage: E=EllipticCurve([a,a+1])
            sage: pts=E.gens()
            sage: len(pts)
            1
            sage: pts[0].order()==E.cardinality()
            True
        """
        try:
            A, G =  self.abelian_group()
            return G
        except AttributeError:
            pass

    def __getitem__(self, n):
        """
        Return the n'th point in self's __points list. This enables users
        to iterate over the curve's point set.

        EXAMPLE::

            sage: E=EllipticCurve(GF(97),[2,3])
            sage: S=E.points()
            sage: E[10]
            (10 : 76 : 1)
            sage: E[15]
            (17 : 10 : 1)
            sage: for P in E: print P.order()
            1
            50
            50
            50
            50
            5
            5
            50
            ...
        """
        return self.points()[n]

    def abelian_group(self, debug=False):
        r"""
        Returns the abelian group structure of the group of points on this
        elliptic curve.

        .. warning::

           The algorithm is definitely *not* intended for use with
           *large* finite fields! The factorization of the orders of
           elements must be feasible. Also, baby-step-giant-step
           methods are used which have space and time requirements
           which are `O(\sqrt{q})`.

        Also, the algorithm uses random points on the curve and hence the
        generators are likely to differ from one run to another; but the
        group is cached so the generators will not change in any one run of
        Sage.

        .. note::

           This function applies to elliptic curves over arbitrary
           finite fields. The related function
           abelian_group_prime_field() uses the pari script, for prime
           fields only; it is now obsolete

        INPUT:


        -  ``debug`` - (default: False): if True, print
           debugging messages


        OUTPUT:

        - an abelian group

        - tuple of images of each of the generators of the abelian
          group as points on this curve

        AUTHORS:

        - John Cremona

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C10, ...

        ::

            sage: E=EllipticCurve(GF(41),[2,5])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C22 x C2, ...

        ::

            sage: F.<a>=GF(3^6,'a')
            sage: E=EllipticCurve([a^4 + a^3 + 2*a^2 + 2*a, 2*a^5 + 2*a^3 + 2*a^2 + 1])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C26 x C26, ...

        ::

            sage: F.<a>=GF(101^3,'a')
            sage: E=EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C1031352, ...

        The group can be trivial::

            sage: E=EllipticCurve(GF(2),[0,0,1,1,1])
            sage: E.abelian_group()
            (Trivial Abelian Group, ())

        Of course, there are plenty of points if we extend the field::

            sage: E.cardinality(extension_degree=100)
            1267650600228231653296516890625

        This tests the patch for trac #3111, using 10 primes randomly
        selected::

            sage: E = EllipticCurve('389a')
            sage: for p in [5927, 2297, 1571, 1709, 3851, 127, 3253, 5783, 3499, 4817]:
            ...       G = E.change_ring(GF(p)).abelian_group()
            sage: for p in prime_range(10000):           #long time (~20s)
            ...       if p != 389:
            ...           G=E.change_ring(GF(p)).abelian_group()

        This tests that the bug reported in trac #3926 has been fixed::

            sage: K.<i> = QuadraticField(-1)
            sage: OK = K.ring_of_integers()
            sage: P=K.factor(10007)[0][0]
            sage: OKmodP = OK.residue_field(P)
            sage: E = EllipticCurve([0,0,0,i,i+3])
            sage: Emod = E.change_ring(OKmodP); Emod
            Elliptic Curve defined by y^2  = x^3 + ibar*x + (ibar+3) over Residue field in ibar of Fractional ideal (10007)
            sage: Emod.abelian_group() #random generators
            (Multiplicative Abelian Group isomorphic to C50067594 x C2,
            ((3152*ibar + 7679 : 7330*ibar + 7913 : 1), (8466*ibar + 1770 : 0 : 1)))
        """
        if not debug:
            # if we're in debug mode, always recalculate
            try:
                return self.__abelian_group
            except AttributeError:
                pass

        k = self.base_field()
        q = k.order()
        p = k.characteristic()
        d = k.degree()
        j = self.j_invariant()
        if d>1:
            d = j.minimal_polynomial().degree()


        # Before computing the group structure we compute the
        # cardinality.  While this is not strictly necessary, it makes
        # the code simpler and also makes computation of orders of
        # points faster.

        # j=0,1728

        if j==k(0):
            N = self._cardinality_with_j_invariant_0()
        if j==k(1728):
            N = self._cardinality_with_j_invariant_1728()

        bounds = Hasse_bounds(q)
        lower, upper = bounds
        if debug:
            print "Lower and upper bounds on group order: [",lower,",",upper,"]"

        try:
            N=self._order
            if debug:
                print "Group order already known to be ",N
        except:
            if (q<50):
                if debug:
                    print "Computing group order naively"
                N=self.cardinality_exhaustive()
            elif d==1:
                if debug:
                    print "Computing group order using SEA"
                N=self.cardinality(algorithm='sea')
            else:
                if debug:
                    print "Computing group order using bsgs"
                N=self.cardinality_bsgs()
            if debug:
                print "... group order = ",N

        self._order=N
        plist = N.prime_factors()
        P1=self(0)
        P2=self(0)
        n1= Integer(1)
        n2= Integer(1)
        P1._order=n1
        P2._order=n2
        npts = 0

        # At all stages the current subgroup is generated by P1, P2 with
        # orders n1,n2 which are disjoint.  We stop when n1*n2 == N

        if debug:
            "About to start generating random points"

        while n1*n2 != N:
            if debug:
                "Getting a new random point"
            Q = self.random_point()
            while Q.is_zero(): Q = self.random_point()
            npts += 1
            if debug:
                print "Q = ",Q,":",
                print " Order(Q) = ", Q.order()

            Q1=n1*Q;

            if Q1.is_zero() and npts>=10: # then P1,n1 will not change but we may increase n2
                if debug: print "Case 2: n2 may increase"
                n1a = 1; n1b = n1
                P1a = P1
                n1a = n1.prime_to_m_part(N//n1)
                n1b = n1//n1a
                Q = n1a*Q       # has order | n1b
                P1a = n1a*P1    # has order = n1b
                if debug: print "n1a=",n1a
                a = None
                for m in n1b.divisors():
                    try:
                        a = generic.bsgs(m*P1a,m*Q,(0,(n1b//m)-1),operation='+')
                        break
                    except ValueError:
                        pass
                assert a != None
                a *= (m*n1a)
                if debug: print "linear relation gives m=",m,", a=",a
                if debug: assert m*Q==a*P1
                if m>1: # else Q is in <P1>
                    Q=Q-(a//m)*P1; # has order m and is disjoint from P1
                    if debug: assert Q.order()==m
                    Q._order=m
                    if n2==1: # this is our first nontrivial P2
                        P2=Q
                        n2=m
                        if debug:
                            print "Adding second generator ",P2," of order ",n2
                            print "Subgroup order now ",n1*n2,"=",n1,"*",n2
                    else:     # we must merge P2 and Q:
                        oldn2=n2 # holds old value
                        P2,n2=generic.merge_points((P2,n2),(Q,m),operation='+');
                        if debug: assert P2.order()==n2
                        P2._order=n2
                        if debug:
                            if n2>oldn2:
                                print "Replacing second generator by ",P2,
                                print " of order ",n2, "  gaining index ",n2//oldn2
                                print "Subgroup order now ",n1*n2,"=",n1,"*",n2
            elif not Q1.is_zero(): # Q1 nonzero: n1 will increase
                if debug:  print "Case 1: n1 may increase"
                oldn1=n1
                if n2>1:
                    P3=(n1//n2)*P1  # so P2,P3 are a basis for n2-torsion
                    if debug: assert P3.order()==n2
                    P3._order=n2
                    if debug: print "storing generator ",P3," of ",n2,"-torsion"
                m = generic.order_from_multiple(Q,N,plist,operation='+')
                P1,n1=generic.merge_points((P1,n1),(Q,m))
                if debug: assert P1.order()==n1
                P1._order=n1
                if debug:
                    print "Replacing first  generator by ",P1," of order ",
                    print n1,", gaining index ",n1//oldn1
                    print "Subgroup order now ",n1*n2,"=",n1,"*",n2

                # Now replace P2 by a point of order n2 s.t. it and
                # (n1//n2)*P1 are still a basis for n2-torsion:
                if n2>1:
                    a,m = generic.linear_relation(P1,P3,operation='+')
                    if debug: print "linear relation gives m=",m,", a=",a
                    P3=P3-(a//m)*P1
                    if debug: assert P3.order()==m
                    P3._order=m
                    if debug: print "First  P2 component =",P3
                    if m==n2:
                        P2=P3
                    else:
                        a,m = generic.linear_relation(P1,P2,operation='+')
                        if debug: print "linear relation gives m=",m,", a=",a
                        P2=P2-(a//m)*P1;
                        if debug: assert P2.order()==m
                        P2._order=m
                        if debug: print "Second  P2 component =",P2
                        P2,n2=generic.merge_points((P2,n2),(P3,m))
                        if debug: assert P2.order()==n2
                        P2._order=n2
                        if debug: print "Combined P2 component =",P2

            if debug:
                if P1.order()!=n1:
                    print "Generator P1 = ",P1," has order ",P1.order(),
                    print " and not ",n1
                    raise ValueError
                if P2.order()!=n2:
                    print "Generator P2 = ",P2," has order ",P2.order()
                    print " and not ",n2
                    raise ValueError
                if n2>1:
                    if generic.linear_relation(P1,P2,operation='+')[1]!=n2:
                        print "Generators not independent!"
                        raise ValueError
                print "Generators: P1 = ",P1," of order ",n1,
                print ", P2 = ",P2," of order ",n2
                print "Subgroup order is now ",n1*n2,"=",n1,"*",n2

        # Finished: record group order, structure and generators

        self._order = n1*n2
        if n1==1:
            self.__abelian_group = AbelianGroup([]), ()
        else:
            if n2==1:
                self.__abelian_group = AbelianGroup([n1]), (P1,)
            else:
                self.__abelian_group = AbelianGroup([n1,n2]), (P1,P2)
        return self.__abelian_group
