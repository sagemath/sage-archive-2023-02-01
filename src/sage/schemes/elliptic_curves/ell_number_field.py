"""
Elliptic curves over number fields

EXAMPLES:
    sage: k.<i> = NumberField(x^2+1)
    sage: E = EllipticCurve([i,2])
    sage: E.j_invariant()
    -23328/365*i + 864/365
    sage: E.simon_two_descent()
    (1, 1, [(2*i : -2*i + 2 : 1)])
    sage: P = E([2*i,-2*i+2])
    sage: P+P
    (15/32*i + 3/4 : 139/256*i + 339/256 : 1)

"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                          William Stein   <wstein@gmail.com>
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


import sage.rings.ring as ring
import sage.rings.all as rings
from ell_field import EllipticCurve_field
import ell_point
import constructor
from sage.rings.all import PolynomialRing, ZZ, QQ
from sage.rings.arith import lcm

from sage.misc.misc import prod

import sage.databases.cremona

from gp_simon import simon_two_descent
from sage.libs.pari.all import pari


class EllipticCurve_number_field(EllipticCurve_field):
    """
    Elliptic curve over a padic field.
    """
    def __init__(self, x, y=None):
        if y is None:
            if isinstance(x, list):
                ainvs = x
                field = ainvs[0].parent()
        else:
            if isinstance(y, str):
                field = x
                X = sage.databases.cremona.CremonaDatabase()[y]
                ainvs = [field(a) for a in X.a_invariants()]
            else:
                field = x
                ainvs = y
        if not (isinstance(field, ring.Ring) and isinstance(ainvs,list)):
            raise TypeError

        EllipticCurve_field.__init__(self, [field(x) for x in ainvs])

        self._point_class = ell_point.EllipticCurvePoint_field
        self._genus = 1

    def simon_two_descent(self, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Computes (probably) the rank of the Mordell-Weil group,
        with certainty the rank of the 2-Selmer group, and a list
        of independent points on the Weierstrass model of self.

        If the curve has 2-torsion, only the probable rank is returned.

        \note{The points are not translated back to self only because
        I haven't written code to do this yet.}

        INPUT:
            verbose -- integer, 0,1,2,3; (default: 0), the verbosity level
            lim1    -- (default: 5) limite des points triviaux sur les quartiques
            lim3    -- (default: 50) limite des points sur les quartiques ELS
            limtriv -- (default: 10) limite des points triviaux sur la
                                     courbe elliptique
            maxprob -- (default: 20)
            limbigprime -- (default: 30)  to distinguish between small and large prime
                                          numbers. Use probabilistic tests for large
                                          primes. If 0, don't use probabilistic tests.

        OUTPUT:
            integer -- "probably" the rank of self
            integer -- the 2-rank of the Selmer group
            list    -- list of independent points on the Weierstrass model

        NOTE: For non-quadratic number fields, this code does return, but it takes a long time.

        IMPLEMENTATION: Uses {\bf Denis Simon's} GP/PARI scripts from
                         \url{http://www.math.unicaen.fr/~simon/}

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E.simon_two_descent()
            (2, 2, [(-4 : -4 : 1), (2*a - 10 : -4*a - 48 : 1)])

            sage: K.<a> = NumberField(x^2 + 7, 'a')
            sage: E = EllipticCurve(K, [0,0,0,1,a]); E
            Elliptic Curve defined by y^2  = x^3 + x + a over Number Field in a with defining polynomial x^2 + 7
            sage: v = E.simon_two_descent(verbose=1); v
            courbe elliptique : Y^2 = x^3 + Mod(3*y, y^2 + 7)*x^2 + Mod(-20, y^2 + 7)*x + Mod(-5*y, y^2 + 7)
            points triviaux sur la courbe = [[1, 1, 0], [Mod(-1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            #S(E/K)[2]    = 2
            #E(K)/2E(K)   = 2
            #III(E/K)[2]  = 1
            rang(E/K)     = 1
            listpointsmwr = [[Mod(-1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            (1, 1, [(1/2*a + 3/2 : -a - 2 : 1)])

        A curve with 2-torsion
            sage: K.<a> = NumberField(x^2 + 7, 'a')
            sage: E = EllipticCurve(K, '15a')
            sage: v = E.simon_two_descent(); v  # long time
            (1, -1, [])
        """

        F = self.integral_weierstrass_model()
        a1,a2,a3,a4,a6 = F.a_invariants()
        x = PolynomialRing(self.base_ring(), 'x').gen(0)
        t = simon_two_descent(a2,a4,a6,
                              verbose=verbose, lim1=lim1, lim3=lim3, limtriv=limtriv,
                              maxprob=maxprob, limbigprime=limbigprime)
        prob_rank = rings.Integer(t[0])
        two_selmer_rank = rings.Integer(t[1])
        prob_gens = [F(P) for P in t[2]]
        return prob_rank, two_selmer_rank, prob_gens

    def integral_weierstrass_model(self):
        a1,a2,a3,a4,a6 = self.weierstrass_model().a_invariants()
        # Find minimum d such that a4*d^4 and a6*d^6 in ZZ.
        d = lcm(prod([r**((e+3)//4) for r, e in a4.denominator().factor()]),
                prod([r**((e+5)//6) for r, e in a6.denominator().factor()]))
        # do transformation x -> x/d^2
        #                   y -> y/d^3
        return constructor.EllipticCurve([a4 * d**4, a6 * d**6])



