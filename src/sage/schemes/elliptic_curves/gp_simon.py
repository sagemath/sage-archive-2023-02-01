"""
Denis Simon's PARI scripts
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
from sage.misc.preparser import preparse
from sage.rings.all import RationalField, IntegerRing
QQ = RationalField()
ZZ = IntegerRing()

gp = None
def init():
    global gp
    if gp is None:
        gp = Gp(script_subdirectory='simon')
        gp.read("ell.gp")
        gp.read("ellQ.gp")
        gp.read("qfsolve.gp")
        gp.read("resultant3.gp")


def simon_two_descent(A, B, C, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
    init()
    gp('DEBUGLEVEL=%s; LIM1=%s; LIM3=%s; LIMTRIV=%s; MAXPROB=%s; LIMBIGPRIME=%s;'%(
        verbose, lim1, lim3, limtriv, maxprob, limbigprime))
    s = gp.eval('ans=main(%s,%s,%s);'%(A,B,C))
    if s.find("***") != -1:
        print s
        raise RuntimeError, "An error occured while running Simon's 2-descent program"
    if verbose > 0:
        print s
    v = gp.eval('ans')
    v = preparse(v)
    return eval(v)

