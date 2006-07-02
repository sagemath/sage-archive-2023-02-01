"""
Cremona PARI Scripts

Access to Cremona's PARI scripts via SAGE.
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

from sage.interfaces.gp import Gp
from sage.rings.all import Integer, RealField
R = RealField()

gp = None
def init():
    global gp
    if gp is None:
        gp = Gp(script_subdirectory='cremona')
        gp.read("bgc.gp")
        gp.read("ell_baby.gp")
        gp.read("ell_ff.gp")
        gp.read("ell_weil.gp")
        gp.read("ell_zp.gp")
        gp.eval('debug_group=0;')

def ellanalyticrank(e):
    """
    INPUT:
        e -- five-tuple of integers that define a minimal weierstrass equation
    OUTPUT:
        integer -- the ("computed") analytic rank r of E
    """
    init()
    cmd = "ellanalyticrank(ellinit(%s),0)"%e
    x = gp.eval(cmd)
    if x.find("***") != -1:
        raise RuntimeError, "Error: '%s'"%x
    return Integer(x)

def ellzp(e, p):
    """
    INPUT:
        e -- five-tuple of integers that define an elliptic curve over Z/pZ
        p -- prime
    OUTPUT:
        ...
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
        GP/PARI object
    """
    init()
    return gp("e=ellzpinit(%s,%s);"%(e,p))


################
# allisog.gp
################
gp_allisog = None
def p_isog(e, p):
    global gp_allisog
    if gp_allisog is None:
        gp_allisog = Gp(script_subdirectory='cremona')
        gp_allisog.read("allisog.gp")

    x = gp_allisog.eval('lisogs(ellinit(%s),%s)'%(e,p))
    if x.find("***") != -1:
        raise RuntimeError, "Error: '%s'"%x
    return x

def allisog(e):
    global gp_allisog
    if gp_allisog is None:
        gp_allisog = Gp(script_subdirectory='cremona')
        gp_allisog.read("allisog.gp")

    x = gp_allisog.eval('allisog(ellinit(%s))'%e)
    if x.find("***") != -1:
        raise RuntimeError, "Error: '%s'"%x
    return x

