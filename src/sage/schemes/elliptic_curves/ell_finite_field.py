"""
Elliptic curves over finite fields
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

from ell_field import EllipticCurve_field
from sage.schemes.hyperelliptic_curves.hyperelliptic_finite_field import HyperellipticCurve_finite_field
import sage.rings.ring as ring
from sage.rings.all import Integer, PolynomialRing, ComplexField
import gp_cremona
import sea
from sage.groups.all import AbelianGroup
import ell_point
from sage.calculus.calculus import log
from sage.rings.arith import integer_ceil, integer_floor

import sage.plot.all as plot

import ell_point

import sage.libs.pari
pari = sage.libs.pari.all.pari

class EllipticCurve_finite_field(EllipticCurve_field, HyperellipticCurve_finite_field):
    """
    Elliptic curve over a finite field.
    """
    def __init__(self, x, y=None):
        if isinstance(x, list):
            ainvs = x
            field = ainvs[0].parent()
        else:
            field = x
            ainvs = y
        if not (isinstance(field, ring.Ring) and isinstance(ainvs,list)):
            raise TypeError

        EllipticCurve_field.__init__(
            self, [field(x) for x in ainvs])

        self._point_class = ell_point.EllipticCurvePoint_finite_field

    def _pari_(self):
        try:
            return self.__pari
        except AttributeError:
            pass
        F = self.base_ring()
        self.__pari = pari('ellinit(Mod(1,%s)*%s)'%(F.characteristic(), [b._pari_() for b in self.ainvs()]))
        return self.__pari

    def _magma_init_(self):
        k = self.base_ring()
        kmn = k._magma_().name()
        return 'EllipticCurve([%s|%s])'%(kmn,','.join([x._magma_init_() for x in self.ainvs()]))

    def _gp(self):
        """
        Return an elliptic curve in a GP/PARI interpreter with all Cremona's code
        for finite fields preloaded.

        The base field must have prime order.
        """
        try:
            return self.__gp
        except AttributeError:
            pass
        F = self.base_ring()
        if not F.is_prime():
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
        #G.xmin(0)
        #G.ymin(0)
        return G

    def _points_over_prime_field(self):
        # TODO, eliminate when polynomial calling is fast
        G, pts = self.abelian_group()
        points = [self(0)]
        if len(pts) == 0:
            return points
        P = pts[0]
        Q = P
        while Q[2] != 0:
            points.append(Q)
            Q += P
        n = len(points)

        if len(pts) > 1:
            P = pts[1]
            Q = P
            while Q[2] != 0:
                for R in points[:n]:
                    points.append(R + Q)
                Q += P
        return points

    def points(self):
        r"""
        All the points on this elliptic curve.

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
          [(0 : 1 : 0), (4 : 1 : 1), (1 : 0 : 1), (4 : 4 : 1)]

          sage: K = GF(p**2,'a')
          sage: E = E.change_ring(K)
          sage: len(E.points())
          32
          sage: (p + 1)**2 - a_sub_p**2
          32
          sage: E.points()
          [(0 : 1 : 0), (0 : 2*a + 4 : 1), (0 : 3*a + 1 : 1), (4*a : 0 : 1), (a + 3 : 4*a : 1), (a + 3 : a : 1), (a + 2 : 4*a + 4 : 1), (a + 2 : a + 1 : 1), (a + 4 : 0 : 1), (2 : 2*a + 4 : 1), (2 : 3*a + 1 : 1), (2*a + 4 : 4*a + 4 : 1), (2*a + 4 : a + 1 : 1), (4*a + 4 : a + 4 : 1), (4*a + 4 : 4*a + 1 : 1), (4 : 1 : 1), (4 : 4 : 1), (a : 1 : 1), (a : 4 : 1), (4*a + 3 : a + 3 : 1), (4*a + 3 : 4*a + 2 : 1), (4*a + 1 : 1 : 1), (4*a + 1 : 4 : 1), (3 : 2*a + 4 : 1), (3 : 3*a + 1 : 1), (2*a : 3*a : 1), (2*a : 2*a : 1), (3*a + 1 : a + 3 : 1), (3*a + 1 : 4*a + 2 : 1), (3*a + 2 : 2*a + 3 : 1), (3*a + 2 : 3*a + 2 : 1), (1 : 0 : 1)]
        """
        try:
            return self.__points
        except AttributeError: pass

        if self.base_ring().is_prime():
            self.__points = self._points_over_prime_field() # _points_cache_sqrt
        else:
            self.__points = self._points_fast_sqrt()

        return self.__points

    def random_element(self):
        """
        Returns a random point on this elliptic curve.

        Returns the point at infinity with probability $1/(\#k+1)$
        where $k$ is the base field.

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
        if random.random() <= 1/float(k.order()+1):
            return self(0)
        a1, a2, a3, a4, a6 = self.ainvs()

        if k.characteristic() == 2:
            P = PolynomialRing(k,'y')
            y = P.gen()
            while True:
                x = k.random_element()
                f = y**2 + a1*x*y + a3*y - x**3 + a2*x**2 + a4*x + a6
                roots = f.roots()
                n = len(roots)
                if n:
                    y = roots[random.randint(0,n-1)][0]
                    return self([x,y])
        else:
            while True:
                x = k.random_element()
                d = 4*x**3 + (a1**2 + 4*a2)*x**2 + (2*a3*a1 + 4*a4)*x + (a3**2 + 4*a6)
                try:
                    m = d.sqrt(extend=False)
                    if random.random() < 0.5:
                        m = -m
                    y = (-(a1*x + a3) + m) / k(2)
                    return self([x,y])
                except ValueError:
                    pass

    def trace_of_frobenius(self):
        """
        Return the trace of Frobenius acting on this elliptic curve.
        """
        q = self.base_field().order()
        return q + 1 - self.cardinality()

    def cardinality(self, algorithm='heuristic', early_abort=False, disable_warning=False):
        r"""
        Return the number of points on this elliptic curve over a
        finite extension of the base field.

        \note{If the cardinality of the base field is not prime and
        the coefficients of the Weierstrass form of self are not in
        the prime subfield this function literally enumerates the
        points and counts them. It's so stupid, it prints a
        warning. You can disable the warning with the disable_warning
        flag.}

        INPUT:
            algorithm    -- string (default: 'heuristic')
                  'heuristic' -- use a heuristic to choose between bsgs and sea.
                  'bsgs' -- use the baby step giant step method as implemented in
                            PARI via the C-library function ellap.
                  'sea'  -- use sea.gp as implemented in PARI by Christophe
                            Doche and Sylvain Duquesne.
                  'all'  -- compute cardinality with both bsgs and sea and
                            return result if they agree or raise a RuntimeError
                            if they do not.
            early_abort  -- bool (default: False); this is used only by sea.
                            if True, stop early if a small factor of the order is found.
        OUTPUT:
            an integer

        The cardinality is \emph{not} cached.

        EXAMPLES:
            sage: EllipticCurve(GF(4,'a'),[1,2,3,4,5]).cardinality()
            8
            sage: k.<a> = GF(3^3)
            sage: l = [a^2 + 1, 2*a^2 + 2*a + 1, a^2 + a + 1, 2, 2*a]
            sage: EllipticCurve(k,l).cardinality()
            WARNING: Using very very stupid algorithm for counting
            points over non-prime finite field. Please rewrite.
            See the file ell_finite_field.py.
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
        """
        N = 0
        if self.base_ring().degree() == 1:
            p = self.base_ring().cardinality()
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
        else:
            # TODO: Instead of only computing in p^1 if possible, also
            # try to compute in p^m if m|n and q is p^n.
            k = self.base_ring()
            p = k.characteristic()
            degree = k.degree()
            q = k.order()
            try:
                A = []
                for a in map(int, self.a_invariants()):
                    if a%p == a:
                        A.append(a)
                    else:
                        raise TypeError
                E = EllipticCurve_finite_field(self.base_ring().prime_subfield(), A)
                prec = integer_ceil(degree*log(p)/log(2))+3
                C = ComplexField(prec)
                t = C(p+1-E.cardinality(algorithm=algorithm,early_abort=early_abort))
                alphap = t/2 + (-p+t**2/4).sqrt()
                apd = alphap**degree
                bpd = alphap.conjugate()**degree
                M = q + 1 - apd - bpd
                N = integer_floor(M.real())

            except TypeError:
                pass # cannot compute over GF(p)

        if N == 0:
            if not disable_warning:
                print "WARNING: Using very very stupid algorithm for counting "
                print "points over non-prime finite field. Please rewrite."
                print "See the file ell_finite_field.py."
            p = self.base_field().cardinality()
            N = len(self.points())

        return N

    order = cardinality # alias

    def _cremona_abgrp_data(self):
        E = self._gp()
        gp = E.parent()
        return eval(gp.eval('[%s.isotype, lift(%s.generators)]'%(E.name(), E.name())))

    def abelian_group(self):
        """
        Returns the abelian group structure of the group of points on
        this elliptic curve.

        OUTPUT:
            -- an abelian group
            -- tuple of images of each of the generators of the
               abelian group as points on this curve

        AUTHOR: John Cremona

        NOTE: And finally...from what I remember of those gp programs,
        the algorithm I was using (which is not exactly the same as
        others use, for example it is not the same as described in
        Cohen's book) is entirely valid for general finite fields.  It
        would just be a question of replacing Mod(a,p) by the
        appropriate alternative.  As I understand it, pari has no
        trouble doing arithmetic over general finite fileds, with
        elements represented as Mod(a,b) where b is a suitable
        irreducible polynomial mod p; though the output format is hard
        to read then (usually improved by allying lift() twice).

        One word of warning though: the underlying algorithm is
        definitely *not* intended for use with *large* finite fields!
        For example, I would regard the factorization of the order of
        an element as inexpensive.
        """
        try:
            return self.__abelian_group
        except AttributeError:
            pass
        if self.base_ring().degree() == 1:
            I, G = self._cremona_abgrp_data()
            A = AbelianGroup(I)
            G = tuple(reversed([self(P) for P in G]))
        else:
            raise NotImplementedError
        self.__abelian_group = A, G
        return A, G

    def __getitem__(self, n):
        return self.points()[n]


