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

from sage.structure.parent_gens import localvars

from sage.interfaces.gp import Gp
from sage.misc.sage_eval import sage_eval
from sage.misc.randstate import current_randstate
from sage.rings.all import QQ
from constructor import EllipticCurve

gp = None
def init():
    """
    Function to initialize the gp process
    """
    global gp
    if gp is None:
        gp = Gp(script_subdirectory='simon')
        gp.read("ell.gp")
        gp.read("ellQ.gp")
        gp.read("qfsolve.gp")
        gp.read("resultant3.gp")


def simon_two_descent(E, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
    """
    Interface to Simon's gp script for two-descent.

    .. NOTE::

       Users should instead run E.simon_two_descent()

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.gp_simon
        sage: E=EllipticCurve('389a1')
        sage: sage.schemes.elliptic_curves.gp_simon.simon_two_descent(E)
        [2, 2, [(1 : 0 : 1), (-11/9 : 28/27 : 1)]]

    TESTS::

        sage: E = EllipticCurve('37a1').change_ring(QuadraticField(-11,'x'))
        sage: E.simon_two_descent()
        (1, 1, [(-1 : 0 : 1)])

    """
    init()

    current_randstate().set_seed_gp(gp)

    K = E.base_ring()
    K_orig = K
    # The following is to correct the bug at \#5204: the gp script
    # fails when K is a number field whose generator is called 'x'.
    if not K is QQ:
        K = K.change_names('a')
    E_orig = E
    E = EllipticCurve(K,[K(list(a)) for a in E.ainvs()])
    F = E.integral_model()

    if K != QQ:
        # Simon's program requires that this name be y.
        with localvars(K.polynomial().parent(), 'y'):
            gp.eval("K = bnfinit(%s);" % K.polynomial())
            if verbose >= 2:
                print "K = bnfinit(%s);" % K.polynomial()
        gp.eval("%s = Mod(y,K.pol);" % K.gen())
        if verbose >= 2:
            print "%s = Mod(y,K.pol);" % K.gen()

    if K == QQ:
        cmd = 'ellrank([%s,%s,%s,%s,%s]);' % F.ainvs()
    else:
        cmd = 'bnfellrank(K, [%s,%s,%s,%s,%s]);' % F.ainvs()

    gp('DEBUGLEVEL_ell=%s; LIM1=%s; LIM3=%s; LIMTRIV=%s; MAXPROB=%s; LIMBIGPRIME=%s;'%(
        verbose, lim1, lim3, limtriv, maxprob, limbigprime))

    if verbose >= 2:
        print cmd
    s = gp.eval('ans=%s;'%cmd)
    if s.find("***") != -1:
        raise RuntimeError, "\n%s\nAn error occurred while running Simon's 2-descent program"%s
    if verbose > 0:
        print s
    v = gp.eval('ans')
    if v=='ans': # then the call to ellrank() or bnfellrank() failed
        return 'fail'
    if verbose >= 2:
        print "v = ", v
    # pari represents field elements as Mod(poly, defining-poly)
    # so this function will return the respective elements of K
    def _gp_mod(*args):
        return args[0]
    ans = sage_eval(v, {'Mod': _gp_mod, 'y': K.gen(0)})
    inv_transform = F.isomorphism_to(E)
    ans[2] = [inv_transform(F(P)) for P in ans[2]]
    ans[2] = [E_orig([K_orig(list(c)) for c in list(P)]) for P in ans[2]]
    return ans

