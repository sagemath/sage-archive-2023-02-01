"""
Elliptic curves over finite fields
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

import random

from ell_field import EllipticCurve_field
import sage.rings.ring as ring
from sage.rings.all import Integer
import gp_cremona
import sea
from sage.groups.all import AbelianGroup
import ell_point

import sage.plot.all as plot

import ell_point

import sage.libs.pari
pari = sage.libs.pari.all.pari

class EllipticCurve_finite_field(EllipticCurve_field):
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

    def _plot_(self, **args):
        """
        Draw a graph of this elliptic curve over a prime finite field.

        INPUT:
            **args -- all other options are passed to the circle
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
                G += plot.point(P, **args)
        #G.xmin(0)
        #G.ymin(0)
        return G

    def __points_over_prime_field(self):
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


    def __points_over_arbitrary_field(self):
        # TODO -- rewrite this insanely stupid implementation!!
        print "WARNING: Using very very stupid algorithm for finding points over"
        print "non-prime finite field.  Please rewrite.  See the file ell_finite.field.py."
        # The best way to rewrite is to extend Cremona's code (either gp or mwrank) so
        # it works over non-prime fields (should be easy), then generate up the group.
        points = [self(0)]
        for x in self.base_field():
            for y in self.base_field():
                try:
                    points.append(self([x,y]))
                except TypeError:
                    pass
        return points

    def points(self):
        r"""
        All the points on this elliptic curve.

        If the base field is another other than $\FF_p$, this function currently
        uses a VERY VERY naive algorithm.   The problem is that Cremona's gp
        script \code{ell_zp.gp} currently only treats the case of prime fields.
        If you need more general functionality, please volunteer to implement
        it, either by extending \code{ell_zp.gp} (if you know PARI/GP well) or
        by writing new \sage code.
        """
        try:
            return self.__points
        except AttributeError: pass

        if self.base_ring().is_prime():
            self.__points = self.__points_over_prime_field()
        else:
            self.__points = self.__points_over_arbitrary_field()

        return self.__points

    def random_element(self):
        """
        Returns a random point on this elliptic curve.

        Returns the point at infinity with probability $1/(\#k+1)$
        where $k$ is the base field.
        """
        k = self.base_field()
        if random.random() <= 1/float(k.order()+1):
            return self(0)
        a1, a2, a3, a4, a6 = self.ainvs()
        while True:
            x = k.random_element()
            d = 4*x**3 + (a1**2 + 4*a2)*x**2 + (2*a3*a1 + 4*a4)*x + (a3**2 + 4*a6)
            try:
                m = d.square_root()
                y = (-(a1*x + a3) + m) / k(2)
                return self([x,y])
            except ValueError:
                pass

    def cardinality(self, algorithm='heuristic', early_abort=False):
        r"""
        Return the number of points on this elliptic curve over this
        finite field.

        \note{If the cardinality of the base field is not prime, this
        function currently uses a very very naive enumeration of all
        points.  It's so stupid, that it prints a warning.}

        INPUT:
            algorithm -- string (default: 'heuristic')
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

        \note{'sea' doesn't work in Windows XP under Cygwin (as of 2005-12-06).}

        The cardinality is \emph{not} cached.

        EXAMPLES:
            sage: EllipticCurve(GF(4),[1,2,3,4,5]).cardinality()
            WARNING: Using very very stupid algorithm for finding points over
            non-prime finite field.  Please rewrite.  See the file ell_finite.field.py.
            8
            sage: EllipticCurve(GF(9),[1,2,3,4,5]).cardinality()
            WARNING: Using very very stupid algorithm for finding points over
            non-prime finite field.  Please rewrite.  See the file ell_finite.field.py.
            16
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

        if N == 0:
            N = len(self.points())
        self.__cardinality = Integer(N)
        return N

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


