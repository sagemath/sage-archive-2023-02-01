"""
Denis Simon's PARI scripts
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

from __future__ import with_statement


from sage.structure.parent_gens import localvars

from sage.interfaces.gp import Gp
from sage.misc.sage_eval import sage_eval
from sage.rings.all import PolynomialRing, ZZ, QQ

gp = None
def init(K):
    global gp
    if gp is None:
        gp = Gp(script_subdirectory='simon')
        gp.read("ell.gp")
        if K == QQ:
            gp.read("ellQ.gp")
        gp.read("qfsolve.gp")
        gp.read("resultant3.gp")


def simon_two_descent(A, B, C, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
    K = A.parent()
    init(K)
    shift = 0
    if K == QQ:
        cmd = 'main(%s,%s,%s);'%(A,B,C)
    else:
        # Simon's program requires that the coefficients do NOT all lie in a subfield of K
        if A.minpoly().degree() != K.degree():
            shift = a = K.gen()
            x = PolynomialRing(K, 'x').gen(0)
            shifted = (x+a)**3 + A*(x+a)**2 + B*(x+a) + C
            C, B, A, _ = shifted
        # Simon's program requires that this name be y.
        with localvars(K.polynomial().parent(), 'y'):
            gp.eval("K = bnfinit(%s);" % K.polynomial())
        gp.eval("%s = Mod(y,K.pol);" % K.gen())
        cmd = 'main(K, %s,%s,%s);'%(A,B,C)

    gp('DEBUGLEVEL=%s; LIM1=%s; LIM3=%s; LIMTRIV=%s; MAXPROB=%s; LIMBIGPRIME=%s;'%(
        verbose, lim1, lim3, limtriv, maxprob, limbigprime))
    s = gp.eval('ans=%s;'%cmd)
    if s.find("***") != -1:
        print s
        raise RuntimeError, "An error occured while running Simon's 2-descent program"
    if verbose > 0:
        print s
    v = gp.eval('ans')
    # pari represents field elements as Mod(poly, defining-poly)
    # so this function will return the respective elements of K
    def gp_mod(*args):
        return args[0]
    ans = sage_eval(v, {'Mod': gp_mod, 'y': K.gen(0)})
    if shift:
        # Undo the x -> x+a translation above
        ans[2] = [[P[0]+shift, P[1]] for P in ans[2]]
    return ans

