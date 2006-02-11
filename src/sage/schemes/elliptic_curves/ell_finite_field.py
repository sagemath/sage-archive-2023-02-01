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

from ell_field import EllipticCurve_field
import sage.rings.ring as ring
from sage.rings.all import Integer
import gp_cremona
import sea
from sage.groups.all import AbelianGroup
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

    def x_pari_(self):
        try:
            return self.__pari
        except AttributeError:
            pass
        F = self.base_ring()
        self.__pari = pari('ellinit(Mod(1,%s)*%s)'%(F.characteristic(), [b._pari_() for b in self.ainvs()]))
        return self.__pari

    def points(self):
        """
        All the points on this elliptic curve.

        TODO: currently VERY VERY stupid implementation
        """
        try:
            return self.__points
        except AttributeError: pass

        points = [self(0)]
        for x in self.base_field():
            for y in self.base_field():
                try:
                    points.append(self([x,y]))
                except TypeError:
                    pass
        self.__points = points
        return self.__points

    def cardinality(self, algorithm='heuristic', early_abort=False):
        r"""
        Return the number of points on this elliptic curve over this
        finite field.

        \note{If the cardinality of the base field is not prime, this
        function currently uses a very naive enumeration of all points.}

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
            8
            sage: EllipticCurve(GF(9),[1,2,3,4,5]).cardinality()
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

    def abelian_group(self):
        """
        Returns the abelian group structure of the group of points on
        this elliptic curve.

        OUTPUT:
            -- an abelian group
            -- tuple of images of each of the generators of the
               abelian group as points on this curve
        """
        try:
            return self.__abelian_group
        except AttributeError:
            pass
        if self.base_ring().degree() == 1:
            I, G = eval(gp_cremona.ellzp(self.a_invariants(), self.base_ring().characteristic()))
            A = AbelianGroup(I)
            G = tuple(reversed([self(P) for P in G]))
        else:
            raise NotImplementedError
        self.__abelian_group = A, G
        return A, G

    def __getitem__(self, n):
        return self.points()[n]
