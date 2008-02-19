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
        """
        Special constructor for elliptic curves over a finite field

        EXAMPLES:
            sage: EllipticCurve(GF(101),[2,3])
            Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field of size 101
            sage: F=GF(101)
            sage: EllipticCurve([F(2),F(3)])
            Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field of size 101
        """
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
        Return a string suitable for passing to Magma (if installed)
        to create the same curve

        EXAMPLES:
            sage: EllipticCurve(GF(41),[2,5])._magma_init_()
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
        #G.xmin(0)
        #G.ymin(0)
        return G

    def _points_over_prime_field(self):
        """
        Return a list of all the points on the curve, for prime fields
        only (see points() for the general case)

        EXAMPLES:
            sage: S=EllipticCurve(GF(97),[2,3])._points_over_prime_field()
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

        if self.base_ring().is_prime_field():
            self.__points = self._points_over_prime_field() # _points_cache_sqrt
        else:
            self.__points = self._points_fast_sqrt()

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
        # The following allows the origin self(0) to be picked with
        # approproate probability:
        if random.random() <= 1/float(k.order()+1):
            return self(0)

        while True:
            try:
                return self.lift_x(k.random_element())
            except:
                pass


    def trace_of_frobenius(self):
        """
        Return the trace of Frobenius acting on this elliptic curve.

        EXAMPLES:
            sage: E=EllipticCurve(GF(101),[2,3])
            sage: E.trace_of_frobenius()
            6
            sage: E=EllipticCurve(GF(11^5,'a'),[2,5])
            sage: E.trace_of_frobenius()
            802
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

        The cardinality is cached:
            sage: E = EllipticCurve(GF(3^100,'a'),[1,2,3,4,5])
            sage: E.cardinality() is E.cardinality()
            True
        """
        try:
            return self.__order
        except AttributeError:
            pass
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

        self.__order = N
        return N

    order = cardinality # alias

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

    def abelian_group(self):
        """
        Returns the abelian group structure of the group of points on
        this elliptic curve.  Only implemented for prime fields, since
        the underlying gp program is thus, but the algorithm there is
        general.

        A word of warning: the underlying algorithm is definitely
        *not* intended for use with *large* finite fields!  The
        factorization of the orders of an elements must be feasible.

        Also, the algorithm uses random points on the curve and hence
        the generators are likely to differ from one run to another.

        OUTPUT:
            -- an abelian group
            -- tuple of images of each of the generators of the
               abelian group as points on this curve

        AUTHOR: John Cremona

        EXAMPLES:
            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E.abelian_group()
            (Multiplicative Abelian Group isomorphic to C10, ...
            sage: EllipticCurve(GF(41),[2,5]).abelian_group()
            (Multiplicative Abelian Group isomorphic to C22 x C2,
            ...

        """
        try:
            return self.__abelian_group
        except AttributeError:
            pass
        if self.base_ring().degree() == 1:

            # In the non-cyclic case, Sage orders the invariants as
            # n1,n2 with n1 a multiple of n2 (true on 2008-01-27!)
            # which agrees with orders of the generators returned by
            # _cremona_abgrp_data(), so there is no need to switch the
            # points as in the past!

            I, G = self._cremona_abgrp_data()
            A = AbelianGroup(I)
            G = tuple([self(P) for P in G])
        else:
            raise NotImplementedError
        self.__abelian_group = A, G
        return A, G

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
            2
            50
            50
            25
            ...
        """
        return self.points()[n]


