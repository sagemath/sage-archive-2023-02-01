"""
Elliptic curves over finite fields

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

import random
from sage.misc.misc import cputime
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

        EXAMPLES:
            sage: EllipticCurve(GF(101),[2,3])
            Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field of size 101

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

        EXAMPLES:
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

    def _magma_init_(self):
        """
        Return a Magma command that creates this curve.

        EXAMPLES:
            sage: EllipticCurve(GF(41),[2,5])._magma_init_() # optional -- requires Magma
            'EllipticCurve([_sage_[1]|GF(41)!0,GF(41)!0,GF(41)!0,GF(41)!2,GF(41)!5])'
            sage: magma(E) # optional -- requires Magma
            Elliptic Curve defined by y^2 = x^3 + 2*x + 5 over GF(41)
       """
        k = self.base_ring()
        kmn = k._magma_().name()
        return 'EllipticCurve([%s|%s])'%(kmn,','.join([x._magma_init_() for x in self.ainvs()]))

    def _gp(self):
        """
        Return an elliptic curve in a GP/PARI interpreter with all
        Cremona's code for finite fields preloaded.  This includes
        generators, which will vary from run to run.

        The base field must have prime order.

        EXAMPLES:
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
            *args, **kwds -- all other options are passed to the circle
                      graphing primitive.

        EXAMPLES:
            sage: E = EllipticCurve(FiniteField(17), [0,1])
            sage: P = plot(E, rgbcolor=(0,0,1))
        """
        R = self.base_ring()
        if not R.is_prime_field():
            raise NotImplementedError

        G = plot.Graphics()
        for P in self.points():
            if not P.is_zero():
                G += plot.point(P, *args, **kwds)
        return G

    def _points_via_group_structure(self):
        """
        Return a list of all the points on the curve, for prime fields
        only (see points() for the general case)

        EXAMPLES:
            sage: S=EllipticCurve(GF(97),[2,3])._points_via_group_structure()
            sage: len(S)
            100
        """
        # TODO, eliminate when polynomial calling is fast
        G, pts = self.abelian_group()

        ni = G.invariants()

        def multiples(P,n):
            H=[self(0)]
            for m in range(1,n):
                H.append(H[-1]+P)
            return H

        H0=multiples(pts[0], ni[0])
        if len(ni)==1:   # cyclic case
            return H0
        else:            # noncyclic
            H1=multiples(pts[1], ni[1])
            return [P+Q for P in H0 for Q in H1]

    def points(self):
        r"""
        All the points on this elliptic curve.  The list of points is
        cached so subsequent calls are free.

        EXAMPLES:
            sage: p = 5
            sage: F = GF(p)
            sage: E = EllipticCurve(F, [1, 3])
            sage: a_sub_p = E.change_ring(QQ).ap(p); a_sub_p
            2

            sage: len(E.points())
            4
            sage: p + 1 - a_sub_p
            4
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (4 : 1 : 1), (4 : 4 : 1)]


            sage: K = GF(p**2,'a')
            sage: E = E.change_ring(K)
            sage: len(E.points())
            32
            sage: (p + 1)**2 - a_sub_p**2
            32
            sage: w = E.points(); w
            [(0 : 1 : 0), (0 : 2*a + 4 : 1), (0 : 3*a + 1 : 1), (1 : 0 : 1), (2 : 2*a + 4 : 1), (2 : 3*a + 1 : 1), (3 : 2*a + 4 : 1), (3 : 3*a + 1 : 1), (4 : 1 : 1), (4 : 4 : 1), (a : 1 : 1), (a : 4 : 1), (a + 2 : a + 1 : 1), (a + 2 : 4*a + 4 : 1), (a + 3 : a : 1), (a + 3 : 4*a : 1), (a + 4 : 0 : 1), (2*a : 2*a : 1), (2*a : 3*a : 1), (2*a + 4 : a + 1 : 1), (2*a + 4 : 4*a + 4 : 1), (3*a + 1 : a + 3 : 1), (3*a + 1 : 4*a + 2 : 1), (3*a + 2 : 2*a + 3 : 1), (3*a + 2 : 3*a + 2 : 1), (4*a : 0 : 1), (4*a + 1 : 1 : 1), (4*a + 1 : 4 : 1), (4*a + 3 : a + 3 : 1), (4*a + 3 : 4*a + 2 : 1), (4*a + 4 : a + 4 : 1), (4*a + 4 : 4*a + 1 : 1)]

        Note that the returned list is an immutable sorted Sequence:
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

        Returns the point at infinity with probability $1/(q+1)$
        where the base field has cardinality $q$.

        EXAMPLES:
            sage: k = GF(next_prime(7^5))
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element()
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

            sage: k.<a> = GF(7^5)
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element()
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

            sage: k.<a> = GF(2^5)
            sage: E = EllipticCurve(k,[a^2,a,1,a+1,1])
            sage: P = E.random_element()
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True
        """
        k = self.base_field()
        # The following allows the origin self(0) to be picked
        if random.random() <= 1/float(k.order()+1):
            return self(0)

        while True:
            try:
                return self.lift_x(k.random_element())
            except:
                pass

    random_point = random_element


    def trace_of_frobenius(self):
        """
        Return the trace of Frobenius acting on this elliptic curve.

        NOTE:
            This computes the curve cardinality, which may be time-consuming.

        EXAMPLES:
            sage: E=EllipticCurve(GF(101),[2,3])
            sage: E.trace_of_frobenius()
            6
            sage: E=EllipticCurve(GF(11^5,'a'),[2,5])
            sage: E.trace_of_frobenius()
            802
        """
        return 1 + self.base_field().order() - self.cardinality()

    def _cardinality_with_jinvariant_0_or_1728(self):
        r"""
        Helper function to handle cardinality when supersingular.

        EXAMPLES:
            sage: E = EllipticCurve(GF(101)(1728))
            sage: E.j_invariant()
            11
            sage: E.cardinality() # indirect doctest
            100

            sage: E = EllipticCurve(GF(103)(1728))
            sage: E.j_invariant()
            80
            sage: E.cardinality() # indirect doctest
            104

            sage: E = EllipticCurve(GF(103^2, 'a')(1728))
            sage: E.j_invariant()
            80
            sage: E.cardinality() # indirect doctest
            10816

            sage: E = EllipticCurve(GF(101)(0))
            sage: E.j_invariant()
            0
            sage: E.cardinality() # indirect doctest
            102

            sage: E = EllipticCurve(GF(103^2, 'a')(0))
            sage: E.j_invariant()
            0
            sage: E.cardinality() # indirect doctest
            10416
        """
        j = self.j_invariant()
        k = self.base_ring()
        q = k.cardinality()
        # j=0, 1728 cases not properly implemented yet

        # Two easy cases:
        if j==k(1728) and q%4==3:
            self._order=Integer(q+1)
            return self._order
        if j==k(0) and q%6==5:
            self._order=Integer(q+1)
            return self._order

        # A quick test to see if the curve's coefficients are all
        # in the prime field:
        if not all(a in k.prime_subfield() for a in self.a_invariants()):
            # resort to basic algorithm:
            self._order = self.cardinality_from_group()
            return self._order

        # OK, curve's coefficients are in prime field:
        E = EllipticCurve_finite_field(k.prime_subfield(), self.a_invariants())
        self._order = E.cardinality(extension_degree=k.degree())
        return self._order

    def cardinality(self, algorithm='heuristic', early_abort=False, disable_warning=False, extension_degree=1):
        r"""
        Return the number of points on this elliptic curve over an
        extension field (default: the base field).

        INPUT:
            algorithm    -- string (default: 'heuristic')
                         -- used only for point counting over prime fields

                  'heuristic' -- use a heuristic to choose between bsgs and sea.
                  'bsgs' -- use the baby step giant step method as implemented in
                            PARI via the C-library function ellap.
                  'sea'  -- use sea.gp as implemented in PARI by Christophe
                            Doche and Sylvain Duquesne.
                  'all'  -- compute cardinality with both bsgs and sea and
                            return result if they agree or raise a RuntimeError
                            if they do not.

            early_abort -- bool (default: False); this is used only by
                            sea.  if True, stop early if a small
                            factor of the order is found.

            extension_degree -- int (default: 1); if the base field is
                            $k=GF(p^n)$ and extension_degree=d, returns
                            the cardinality of $E(GF(p^{n d}))$.

        OUTPUT: an integer

        The cardinality is cached.

        Over prime fields, one of the above algorithms is used.  Over
        non-prime fields, the serious point counting is done on a
        standard curve with the same j-invariant over the field
        GF(p)(j), then lifted to the base_field, and finally account
        is taken of twists.

        EXAMPLES:
            sage: EllipticCurve(GF(4,'a'),[1,2,3,4,5]).cardinality()
            8
            sage: k.<a> = GF(3^3)
            sage: l = [a^2 + 1, 2*a^2 + 2*a + 1, a^2 + a + 1, 2, 2*a]
            sage: EllipticCurve(k,l).cardinality()
            29

            sage: l = [1, 1, 0, 2, 0]
            sage: EllipticCurve(k,l).cardinality()
            38

            An even bigger extension (which we check against Magma):

            sage: EllipticCurve(GF(3^100,'a'),[1,2,3,4,5]).cardinality()
            515377520732011331036459693969645888996929981504
            sage: magma.eval("Order(EllipticCurve([GF(3^100)|1,2,3,4,5]))")    # optional -- requires magma
            '515377520732011331036459693969645888996929981504'


            sage: EllipticCurve(GF(10007),[1,2,3,4,5]).cardinality()
            10076
            sage: EllipticCurve(GF(10007),[1,2,3,4,5]).cardinality(algorithm='sea')
            10076
            sage: EllipticCurve(GF(10007),[1,2,3,4,5]).cardinality(algorithm='bsgs')
            10076
            sage: EllipticCurve(GF(next_prime(10**20)),[1,2,3,4,5]).cardinality(algorithm='sea')
            100000000011093199520

            The cardinality is cached:
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
        try:
            return self._order
        except AttributeError:
            pass

        N = 0
        k = self.base_ring()
        q = k.cardinality()
        p = k.characteristic()
        d = k.degree()

        if d == 1:             # prime field
            if algorithm == 'heuristic':
                if p > 10**18:
                    algorithm = 'sea'
                else:
                    algorithm = 'bsgs'
            if algorithm == 'bsgs':
                E = self._pari_()
                N = p+1 - int(E.ellap(p))
            elif algorithm == 'sea':
                N = sea.ellsea(self.a_invariants(), self.base_ring().characteristic(), \
                               early_abort=early_abort)
            elif algorithm == 'all':
                N1 = self.cardinality('bsgs')
                N2 = self.cardinality('sea')
                if N1 == N2:
                    N = N1
                else:
                    raise RuntimeError, "BUG! Cardinality with bsgs=%s but with sea=%s"%(N1, N2)
            self._order = Integer(N)
            return self._order

        # now k is not a prime field

        # we count points on a standard curve with the same
        # j-invariant defined over the field it generates, then
        # lift to the curve's own field, and finally allow for twists
        j=self.j_invariant()

        if j==k(0) or j==k(1728):
            return self._cardinality_with_jinvariant_0_or_1728()

        # Now the only twists are quadratic, which is simpler

        j_pol=j.minimal_polynomial()
        j_deg=j_pol.degree()

        # if not possible to work over a smaller field:
        if d==j_deg:
            self._order = self.cardinality_from_group()
            return self._order

        kj=GF(p**j_deg,name='a',modulus=j_pol)
        jkj=kj.gen() if j_deg>1 else j_pol.roots(multiplicities=False)[0]

        # recursive call which will do all the real work:
        Ej = EllipticCurve(jkj)
        N=Ej.cardinality(extension_degree=d//j_deg)
        # is curve a (quadratic) twist of the "standard" one?
        if not self.is_isomorphic(EllipticCurve(j)): N=2*(q+1)-N

        self._order = N
        return self._order

    order = cardinality # alias

    def frobenius_polynomial(self):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism of the elliptic curve

        NOTE:
            This computes the curve cardinality, which may be time-consuming.

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_polynomial()
            x^2 - 4*x + 11
        """
        x=polygen(ZZ)
        return x**2-self.trace_of_frobenius()*x+self.base_field().cardinality()

    def frobenius_order(self):
        r"""
        Return the quadratic order Z[phi] where phi is the Frobenius
        endomorphism of the elliptic curve

        NOTE:
            This computes the curve cardinality, which may be time-consuming.

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_order()
            Order in Number Field in phi with defining polynomial x^2 - 4*x + 11
        """
        f = self.frobenius_polynomial().factor()[0][0]
        return ZZ.extension(f,names='phi')

    def frobenius(self):
        r"""
        Return the frobenius of self as an element of a quadratic order

        NOTES:
            This computes the curve cardinality, which may be time-consuming.

            Frobenius is only determined up to conjugacy.

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius()
            phi
            sage: E.frobenius().minpoly()
            x^2 - 4*x + 11
        """
        return self.frobenius_order().fraction_field().gen()

    def cardinality_exhaustive(self):
        r"""
        Return the cardinality of self over the base field.  Simply
        adds up the number of points with each x-coordinate: only used
        for small field sizes!

        EXAMPLES:
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

    def _merge_points(self,P1,n1,P2):
        """
        Given a point P1 of order n1 and another point P2, returns
        (Q,m) where Q is a point Q of order m=lcm(n1,n2) where n2 is
        the order of P2

        EXAMPLES:
            sage: F.<a>=GF(3^6,'a')
            sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
            sage: P=E(2*a^5 + 2*a^4 + a^3 + 2 , a^4 + a^3 + a^2 + 2*a + 2)
            sage: P.order()
            7
            sage: Q=E(2*a^5 + 2*a^4 + 1 , a^5 + 2*a^3 + 2*a + 2 )
            sage: Q.order()
            4
            sage: E._merge_points(P,7,Q)
            ((2*a^2 + 2*a + 1 : a^5 + 2*a^4 + a^3 + 2*a^2 + a + 1 : 1), 28)
        """
        if n1*P2==self(0):
            return (P1,n1)
        n2=P2.order()
        if n1.divides(n2):
            return (P2,n2)
        m,k1,k2 = tidy_lcm(n1,n2);
        m1=n1//k1
        m2=n2//k2
        return (m1*P1+m2*P2, m)

    def cardinality_from_group(self):
        r"""
        Return the cardinality of self over the base field.  Will be
        called by user function cardinality only when necessary,
        i.e. when the j_invariant is not in the prime field.

        This function just calls abelian_group(), so results in the
        group structure and generators being cached as well as the
        group order.

        EXAMPLES:
            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_from_group()
            1020
            sage: E=EllipticCurve(GF(3^4,'a'),[1,1])
            sage: E.cardinality_from_group()
            64
        """

        A, gens = self.abelian_group()
        return self._order


    def _cremona_abgrp_data(self):
        """
        Interface with Cremona's gp program for finding the group
        structure over finite prime fields

        Returns a list whose first element is the list of structure
        constants ([n] or [n1,n2] with n2|n1) and the second is the
        corresponding list of generators

        Users are recommended to use abelian_group() instead.

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E._cremona_abgrp_data()                       # generators random
            [[10], [[0, 7]]]
            sage: EllipticCurve(GF(41),[2,5])._cremona_abgrp_data() # generators random
            [[22, 2], [[35, 8], [23, 0]]]
        """
        E = self._gp()
        gp = E.parent()
        return eval(gp.eval('[%s.isotype, lift(%s.generators)]'%(E.name(), E.name())))

    def abelian_group_prime_field(self):
        """
        Returns the abelian group structure of the group of points on
        this elliptic curve.  Only for prime fields, since is uses an
        underlying gp program which has this restriction.  This is now
        essentially obsolete, replaced by the Sage function
        abelian_group() for arbitrary finite fields.

        WARNING:
            The underlying algorithm is definitely *not* intended for use with
            *large* finite fields!  The factorization of order of elements
            must be feasible.

        NOTE:
            The algorithm uses random points on the curve and hence the
            generators are likely to differ from one run to another; but the
            group is cached so the generators will not change in any one run
            of Sage.

        OUTPUT:
            -- an abelian group
            -- tuple of images of each of the generators of the
               abelian group as points on this curve

        AUTHOR: John Cremona

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E.abelian_group_prime_field()
            (Multiplicative Abelian Group isomorphic to C10, ...
            sage: EllipticCurve(GF(41),[2,5]).abelian_group_prime_field()
            (Multiplicative Abelian Group isomorphic to C22 x C2, ...
        """
        try:
            return self.__abelian_group
        except AttributeError:
            pass

        if self.base_ring().degree() == 1:
            I, G = self._cremona_abgrp_data()
            A = AbelianGroup(I)
            G = tuple([self(P) for P in G])
        else:
            raise NotImplementedError
        self.__abelian_group = A, G
        return A, G

    def gens(self):
        """
        Returns a tuple of length up to 2 of points which generate the
        abelian group of points on this elliptic curve.  See
        abelian_group() for limitations.

        The algorithm uses random points on the curve, and hence the
        generators are likely to differ from one run to another; but
        they are cached so will be consistent in any one run of Sage.

        AUTHOR: John Cremona

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E.gens()                           # random
            sage: EllipticCurve(GF(41),[2,5]).gens() # random
            ...
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
        Return the n'th point in self's __points list.  This enables
        users to iterate over the curve's point set.

        EXAMPLE:
            sage: E=EllipticCurve(GF(97),[2,3])
            sage: S=E.points()
            sage: E[10]           # random
            (29 : 43 : 1)
            sage: E[15]          # random
            (56 : 8 : 1)
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
        Returns the abelian group structure of the group of points on
        this elliptic curve.

        WARNING: The algorithm is definitely *not* intended for use
            with *large* finite fields!  The factorization of the
            orders of elements must be feasible.  Also,
            baby-step-giant-step methods are used which have space and
            time requirements which are $O(\sqrt{q})$.

        Also, the algorithm uses random points on the curve and hence
        the generators are likely to differ from one run to another;
        but the group is cached so the generators will not change in
        any one run of Sage.

        Note: This function applies to elliptic curves over arbitrary
        finite fields.  The related function
        abelian_group_prime_field() uses the pari script,  for prime
        fields only; it is now obsolete

        INPUT:
            -- debug (default: False): if True, print debugging messages

        OUTPUT:
            -- an abelian group
            -- tuple of images of each of the generators of the
               abelian group as points on this curve

        AUTHOR: John Cremona

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C10, ...

            sage: E=EllipticCurve(GF(41),[2,5])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C22 x C2, ...

            sage: F.<a>=GF(3^6,'a')
            sage: E=EllipticCurve([a^4 + a^3 + 2*a^2 + 2*a, 2*a^5 + 2*a^3 + 2*a^2 + 1])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C26 x C26, ...

            sage: F.<a>=GF(101^3,'a')
            sage: E=EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C1031352, ...

        The group can be trivial (but this is the only example of that
        up to isomorphism!)
            sage: E=EllipticCurve(GF(2),[0,0,1,1,1])
            sage: E.abelian_group()
            (Trivial Abelian Group, ())

        Of course, there are plenty of points if we extend the field:
            sage: E.cardinality(extension_degree=100)
            1267650600228231653296516890625
        """
        if not debug:
            # if we're in debug mode, always recalculate
            try:
                return self.__abelian_group
            except AttributeError:
                pass

        k=self.base_field()
        q=k.order()
        lower, upper = Hasse_bounds(q)
        if debug:
            print "Lower and upper bounds on group order: [",lower,",",upper,"]"

        group_order_known = False
        try:
            N=self._order
            if debug:
                print "Group order alrady known to be ",N
            group_order_known = True
            lower=N
            upper=N
        except:
            if (q<75):
                if debug:
                    print "Computing group order naively"
                N=self.cardinality_exhaustive()
                if debug:
                    print "... group order = ",N
                group_order_known = True
                lower=N
                upper=N
                self._order=N
            elif k.is_prime_field() and q>10**18:
                if debug:
                    print "Computing group order using SEA"
                N=self.cardinality(algorithm='sea')
                if debug:
                    print "... group order = ",N
                group_order_known = True
                lower=N
                upper=N
                self._order=N

        if group_order_known and debug:
            print "Lower and upper bounds on group order adjusted ",
            print "to actual order ",lower

        P1=self(0)
        P2=self(0)
        n1= Integer(1)
        n2= Integer(1)
        P1._order=n1
        P2._order=n2
        npts = 0

        # At all stages the current subgroup is generated by P1, P2 with
        # orders n1,n2 which are disjoint.  We stop when n1*n2 >= lower

        if debug:
            "About to start generating random points"
            sys.stdout.flush()

        while n1*n2<lower:
            if debug:
                "Getting a new random point"
                sys.stdout.flush()
            t=cputime()
            Q = self.random_point()
            npts += 1
            if debug:
                print "Q = ",Q,":",
                t=cputime()
                print " Order(Q) = ", Q.order()
                sys.stdout.flush()

            t=cputime()
            Q1=n1*Q;
            sys.stdout.flush()

            if Q1.is_zero() and npts>=10: # then P1,n1 will not change but we may increase n2
                if debug: print "Case 2: n2 may increase"
                m,a = P1.linear_relation(Q)
                if debug: print "linear relation gives m=",m,", a=",a
                if m>1: # else Q is in <P1>
                    Q=Q-(a//m)*P1; # has order m and is disjoint from P1
                    Q._order=m
                    if n2==1: # this is our first nontrivial P2
                        P2=Q
                        n2=m
                        if debug:
                            print "Adding second generator ",P2," of order ",n2
                            print "Group order now ",n1*n2,"=",n1,"*",n2
                    else:     # we must merge P2 and Q:
                        oldn2=n2 # holds old value
                        P2,n2=self._merge_points(P2,n2,Q);
                        P2._order=n2
                        if debug:
                            if n2>oldn2:
                                print "Replacing second generator by ",P2,
                                print " of order ",n2, "  gaining index ",n2//a
                                print "Group order now ",n1*n2,"=",n1,"*",n2
            elif not Q1.is_zero(): # Q1 nonzero: n1 will increase
                if debug:  print "Case 1: n1 may increase"
                oldn1=n1
                if n2>1:
                    P3=(n1//n2)*P1  # so P2,P3 are a basis for n2-torsion
                    P3._order=n2
                    if debug: print "storing generator ",P3," of ",n2,"-torsion"
                P1,n1=self._merge_points(P1,n1,Q)
                P1._order=n1
                if debug:
                    print "Replacing first  generator by ",P1," of order ",
                    print n1,", gaining index ",n1//oldn1
                    print "Group order now ",n1*n2,"=",n1,"*",n2
                # Now replace P2 by a point of order n2 s.t. it and
                # (n1//n2)*P1 are still a basis for n2-torsion:
                if n2>1:
                    m,a = P1.linear_relation(P3)
                    if debug: print "linear relation gives m=",m,", a=",a
                    P3=P3-(a//m)*P1
                    P3._order=m
                    if debug: print "First  P2 component =",P3
                    if m==n2:
                        P2=P3
                    else:
                        m,a = P1.linear_relation(P2)
                        if debug: print "linear relation gives m=",m,", a=",a
                        P2=P2-(a//m)*P1;
                        P2._order=m
                        if debug: print "Second  P2 component =",P2
                        P2,n2=self._merge_points(P2,n2,P3)
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
                if P1.linear_relation(P2)[0]!=n2:
                    print "Generators not independent!"
                    raise ValueError
                print "Generators: P1 = ",P1," of order ",n1,
                print ", P2 = ",P2," of order ",n2
                print "Group order now ",n1*n2,"=",n1,"*",n2

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


# utility function
def tidy_lcm(m,n):
    """
    Utility function: given two positive integers m,n, returns a
    triple (l,m1,n1 such that l=lcm(m,n)=m1*n1 where m1|m, n1|n and
    gcd(m1,n1)=1.  All with no factorization.

    Used to construct an element of order l from elements of orders
    m,n in any group.

    EXAMPLES:
        sage: from sage.schemes.elliptic_curves.ell_finite_field import tidy_lcm
        sage: tidy_lcm(120,36)
        (360, 40, 9)
    """
    m0=m; n0=n
    g=gcd(m,n)
    l=m*n//g      # = lcm(m,n)
    g=gcd(m,n//g) # divisible by those primes which divide n to a
                  # higher power than m

    while not g==1:
        m//=g
        g=gcd(m,g)

    n=l//m;
    assert l==m*n
    assert gcd(m,n)==1
    assert m.divides(m0)
    assert n.divides(n0)
    return (l,m,n)
