"""
Cremona PARI Scripts

Access to Cremona's PARI scripts via SAGE.
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

from sage.interfaces.gp import Gp
from sage.rings.all import Integer, RealField
from sage.misc.randstate import current_randstate
from sage.misc.misc import verbose
R = RealField()

gp = None
def init():
    """
    Function to initialize the gp process
    """
    global gp
    if gp is None:
        gp = Gp(script_subdirectory='cremona')
        gp.read("bgc.gp")
        gp.read("ell_baby.gp")
        gp.read("ell_ff.gp")
        gp.read("ell_weil.gp")
        gp.read("ell_zp.gp")
        gp.eval('debug_group=0;')

def ellanalyticrank_prec(e,prec=None):
    """
    Try to compute analytic rank with precision set to prec.

    INPUT:
        e -- five-tuple of integers that define a minimal weierstrass equation
    OUTPUT:
        integer -- the ("computed") analytic rank r of E

    ALGORITHM:
        Uses Cremona's gp script

    NOTE:
        Users are commended to use EllipticCurve(ai).analytic_rank() instead.

    EXAMPLES:
        sage: sage.schemes.elliptic_curves.gp_cremona.ellanalyticrank_prec([0,-1,1,-10,-20])
        0
        sage: sage.schemes.elliptic_curves.gp_cremona.ellanalyticrank_prec([0,0,1,-1,0])
        1
        sage: sage.schemes.elliptic_curves.gp_cremona.ellanalyticrank_prec([0,1,1,-2,0])
        2
        sage: sage.schemes.elliptic_curves.gp_cremona.ellanalyticrank_prec([0,0,1,-7,6])
        3
        sage: sage.schemes.elliptic_curves.gp_cremona.ellanalyticrank_prec([0,0,1,-7,36])
        4
    """
    init()
    if prec: old_prec = gp.set_real_precision(prec)
    cmd = "ellanalyticrank(ellinit(%s),0)"%e
    x = gp.eval(cmd)
    if prec: gp.set_real_precision(old_prec)
    if x.find("***") != -1:
        raise RuntimeError, "error '%s' running '%s'"%(x,cmd)
    return Integer(x)

def ellanalyticrank(e):
    """
    Try to compute analytic rank.

    INPUT:
        e -- five-tuple of integers that define a minimal weierstrass equation

    OUTPUT:
        integer -- the ("computed") analytic rank r of E

    ALGORITHM:
        Uses Cremona's gp script

    NOTE:
        Users are commended to use EllipticCurve(ai).analytic_rank() instead.

    EXAMPLES:
        sage: sage.schemes.elliptic_curves.gp_cremona.ellanalyticrank([0,-1,1,-10,-20])
        0
    """
    prec = 16
    while True:
        try:
            return ellanalyticrank_prec(e, prec)
        except RuntimeError,msg:
            if 'precision too low' in str(msg):
                prec *= 2
            else:
                raise

def ellzp(e, p):
    """
    INPUT:
        e -- five-tuple of integers that define an elliptic curve over Z/pZ
        p -- prime
    OUTPUT:
        A string containing information about the elliptic curve
        modulo p: group structure and generators.

    NOTE: This is an internal function used in the function
    _abelian_group_data() for curves over finite (prime) field.  Users
    should instead use higher-level funtions -- see examples.

    WARNING: The algorithm uses random points, so the generators in
    the second part of the output will vary from run to run.

    EXAMPLES:
        sage: import sage.schemes.elliptic_curves.gp_cremona
        sage: sage.schemes.elliptic_curves.gp_cremona.ellzp([0,0,1,-7,6],97) #random
        '[[46, 2], [[58, 45], [45, 48]]]
        sage: EllipticCurve(GF(97),[0,0,1,-7,6]).abelian_group() #random
        (Multiplicative Abelian Group isomorphic to C46 x C2,
        ((52 : 13 : 1), (45 : 48 : 1)))
    """
    init()
    cmd = "e=ellzpinit(%s,%s); [e.isotype, lift(e.generators)]"%(e,p)
    x = gp.eval(cmd)
    if x.find("***") != -1:
        raise RuntimeError, "Error: '%s'"%x
    return x

def ellinit(e, p):
    """
    INPUT:
        e -- five-tuple of integers that define an elliptic curve over Z/pZ
        p -- prime
    OUTPUT:
        GP/PARI object representing the elliptic curve modulo p,
        including the following fields:

        E[14] : the group order n = #E(Fp)
        E[15] : factorization of the group order n (as a matrix)
        E[16] : the group structure: [n] if cyclic, or
                                     [n1,n2] with n=n1*n2 and n2|n1
        E[17] : the group generators: [P] is cyclic, or
                                      [P1,P2] with order(Pi)=ni


    EXAMPLES:
        sage: import sage.schemes.elliptic_curves.gp_cremona
        sage: E = sage.schemes.elliptic_curves.gp_cremona.ellinit([0,0,1,-7,6],97)
        sage: E # random generators
        [Mod(0, 97), Mod(0, 97), Mod(1, 97), Mod(90, 97), Mod(6, 97), Mod(0, 97), Mod(83, 97), Mod(25, 97), Mod(48, 97), Mod(45, 97), Mod(32, 97), Mod(33, 97), Mod(63, 97), 92, [2, 2; 23, 1], [46, 2], [[Mod(96, 97), Mod(3, 97)], [Mod(30, 97), Mod(48, 97)]], 0, 0]
        sage: type(E)
        <class 'sage.interfaces.gp.GpElement'>
    """
    init()
    current_randstate().set_seed_gp(gp)
    return gp("e=ellzpinit(%s,%s);"%(e,p))


################
# allisog.gp
################
gp_allisog = None
def p_isog(e, p):
    """
    Return a list of the elliptic curves p-isogenous to e.

    ALGORITHM:
        Uses Cremona's gp script

    EXAMPLES:
        sage: import sage.schemes.elliptic_curves.gp_cremona
        sage: E=EllipticCurve('11a1')
        sage: sage.schemes.elliptic_curves.gp_cremona.p_isog(E.pari_curve(),5)
        '[[0, -1, 1, -7820, -263580], [0, -1, 1, 0, 0]]'
    """
    global gp_allisog
    if gp_allisog is None:
        gp_allisog = Gp(script_subdirectory='cremona')
        gp_allisog.read("allisog.gp")

    x = gp_allisog.eval('lisogs(ellinit(%s),%s)'%(e,p))
    if x.find("***") != -1:
        raise RuntimeError, "Error: '%s'"%x
    return x

def allisog(e):
    """
    Return a list of the curves p-isogenous to e for some prime p.

    OUTPUT: A list of lists [p,curves] where curves is a list of
        elliptic curves p-isogenous to e.

    ALGORITHM:
        Uses Cremona's gp script

    EXAMPLES:
        sage: import sage.schemes.elliptic_curves.gp_cremona
        sage: E=EllipticCurve('14a1')
        sage: sage.schemes.elliptic_curves.gp_cremona.allisog(E.pari_curve())
        '[[2, [[1, 0, 1, -36, -70]]], [3, [[1, 0, 1, -171, -874], [1, 0, 1, -1, 0]]]]'
    """
    global gp_allisog
    if gp_allisog is None:
        gp_allisog = Gp(script_subdirectory='cremona')
        gp_allisog.read("allisog.gp")

    x = gp_allisog.eval('allisog(ellinit(%s))'%e)
    if x.find("***") != -1:
        raise RuntimeError, "Error: '%s'"%x
    return x

