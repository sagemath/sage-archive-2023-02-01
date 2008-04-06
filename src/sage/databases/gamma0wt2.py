"""nodoctest
Table of arithmetic information about modular forms of
weight 2 on Gamma_0(N) for N <= 10000.
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
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

import os

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.arith import GCD
import sage.databases.db


# TODO:
# I messed up and now self.d is i and self.i is d,
# i.e., the degree and number are swapped.


class ModularForm:
    def __init__(self, N, d, i, Wq, r, charpolys, disc):
        self.N = N
        self.d = d
        self.i = i
        self.Wq = Wq
        self.r = r
        self.charpolys = charpolys
        self.disc = disc

    def __repr__(self):
        s = "Modular form: level %s, degree %s, number %s"%(
            self.N,self.d,self.i) + \
            ", Wq: %s, an rank bnd: %s"%(
            self.Wq, self.r)
        return s

    def congruence_multiple(self, f, anemic=True):
        """
        Return an integer C such that any prime of congruence between
        this modular form and f is a divisor of C.

        If anemic=True (the default), include only coefficients a_p
        with p coprime to the levels of self and f.

        This function returns the gcd of the resultants of each of the
        charpolys of coefficients of self and f that are both known.
        """
        C = 0
        N = self.N * f.N
        for p in self.charpolys.keys():
            if p in f.charpolys.keys():
                if (not anemic) or (anemic and N%p!=0):
                    ap = self.charpolys[p]
                    bp = f.charpolys[p]
                    C = GCD(C, ap.resultant(bp))
        return C

    def torsion_multiple(self):
        """
        Returns a multiple of the order of the torsion subgroup of the
        corresponding abelian variety.
        """
        T = 0
        for p in self.charpolys.keys():
            if self.N % p != 0:
                f = self.charpolys[p]
                T = GCD(T, long(f(p+1)))
        return T

class Gamma0Wt2Database(sage.databases.db.Database):
    def __init__(self, read_only=True):
        sage.databases.db.Database.__init__(self,\
             name="gamma0wt2", read_only=read_only)

    def __repr__(self):
        return "Table of arithmetic informationa about newforms of weight 2 on Gamma_0(N)"


def migrate(Nstart, Nstop):
    import conv
    d = Database(read_only=False)
    for N in xrange(Nstart, Nstop):
        print N
        C = conv.newforms(N)
        D = []
        for f in C:
            D.append(sage.databases.gamma0wt2.ModularForm(
                f.N, f.i, f.d,
                f.Wq, f.r,
                f.charpolys, f.disc))
        d[N] = D
    d.commit()

