"""
SEA: Schoof, Elkies, Atkins point counting

Interface to the GP SEA implementation of Christophe Doche and Sylvain
Duquesne.
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
from sage.rings.all import Integer

gp = None
def ellsea(E, p, early_abort=False):
    """
    INPUT:
        E -- list of 5 integers that defines an elliptic curve
        p -- prime number
        early_abort -- bool (default: False); if True an early abort
                       technique is used and the computation is
                       interrupted as soon as a small divisor of the
                       order is detected.  The function then returns
                       0.  This is useful for ruling out curves whose
                       cardinality is divisible by a small prime.
    """
    global gp
    if gp is None:
        gp = Gp(script_subdirectory='SEA')
        gp.eval('allocatemem();allocatemem();allocatemem();allocatemem();allocatemem();allocatemem()')
        gp.read("sea.gp")

    gp.eval('E = ellinit(%s*Mod(1,%s));'%(E,p))
    N = gp.eval("ellsea(E,%s,0,%s)"%(p,int(early_abort)))
    if N.find("*") != -1:
        raise RuntimeError, "Error: '%s'"%N
    return Integer(N)
