"""
MAGMA 3 descent

Access to Michael Stoll's 3-descent via MAGMA, with a slight
change so that Proof is true by default.

This file fires up its on copy of MAGMA when first used.  This copy of
MAGMA attaches files that may change the default behavior of MAGMA
(hence a separate process).
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.interfaces.magma import Magma
from sage.rings.all import Integer

magma = None
def init():
    global magma
    if magma is None:
        magma = Magma(script_subdirectory='stoll')
        magma.attach('3descent.m')

def three_selmer_rank(E, bound=0, method=2):
    """
    INPUT:
        E -- SAGE elliptic curve over Q with irreducible mod 3 representation

        bound -- integer (default: 0); if 0 use a bound determined by Magma;
                                       otherwise, use this bound (use a small
                                       value to see what the answer probably
                                       is, but without proof)

        method -- integer (default: 2) this "parameter specifies how to deal
                  with the global restriction involving the algebra B. If
                  Method = 0, use class and unit groups and ideal
                  factorisation.  If Method = 1, use class and unit groups
                  and reductions mod primes. If Method = 2, use test on cubes."
                  (Quoting from Stoll's documentation for the Magma function.)

    OUTPUT:
        integer -- the rank of the 3-selmer group of E, i.e., the
                   dimension over F_3 of Sel^{(3)}(E/Q)

    EXAMPLES:
        sage: from sage.schemes.elliptic_curves.magma_3descent import three_selmer_rank
        sage: three_selmer_rank(EllipticCurve('11a'))   # optional - magma
        0
    """
    init()
    cmd = '_, d, _, _ := ThreeSelmerGroup(%s : Bound := %s, Method := %s)'%(
        E._magma_init_(), bound, method)
    magma.eval(cmd)
    x = magma.eval('d')
    return Integer(x)
