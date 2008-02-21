"""
Complex multiplication for elliptic curves
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

from sage.interfaces.all import magma
from sage.rings.all import (Integer,
                            RationalField,
                            IntegerRing,
                            is_fundamental_discriminant,
                            PolynomialRing)

def hilbert_class_polynomial(D):
    r"""
    Returns the Hilbert class polynomial of for the discriminant $D$,
    computed using \emph{Magma}.

    \note{This function will not work if Magma is not installed on
    your computer. The Magma function that it calls was implemented
    by David Kohel.}

    INPUT:
        D -- a negative integer congruent to 0 or 1 modulo 4.

    OUTPUT:
        A polynomial over the integers.

    EXAMPLES:
        sage: hilbert_class_polynomial(-4)      # optional MAGMA required
        x - 1728
        sage: hilbert_class_polynomial(-7)      # optional MAGMA required
        x + 3375
        sage: hilbert_class_polynomial(-23)     # optional MAGMA required
        x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
        sage: hilbert_class_polynomial(-37*4)   # optional MAGMA required
        x^2 - 39660183801072000*x - 7898242515936467904000000
    """
    D = Integer(D)
    if D >= 0:
        raise ValueError, "D (=%s) must be negative"%D
    #if not is_fundamental_discriminant(D):
    #    raise ValueError, "D (=%s) must be a fundamental discriminant"%D
    if not (D%4 in [0,1]):
         raise ValueError, "D (=%s) must be a discriminant"%D
    magma.eval("R<x> := PolynomialRing(IntegerRing())")
    f = str(magma.eval("HilbertClassPolynomial(%s)"%D))
    x = PolynomialRing(IntegerRing(), name='x').gen()
    f = f.replace('^','**')
    return eval(f)

def cm_j_invariants(K):
    r"""
    Return a list of all j-invariants of CM elliptic curves such that $j \in K$.

    INPUT:
        K -- a number field (but only implemented for K=QQ)

    OUTPUT:
        list -- of CM j-invariants

    \note{This is currently only implemented for the rationals.  David Kohel
    has large tables for other fields, but they are not in \sage yet.}

    EXAMPLE:
    sage: cm_j_invariants(QQ)
    [0, 54000, -12288000, 1728, 287496, -3375, 16581375, 8000, -32768, -884736, -884736000, -147197952000, -262537412640768000]
    """
    if K == RationalField():
        return [Integer(x) for x in [0, 54000, -12288000, 1728, \
                               287496, -3375, 16581375, 8000, \
                               -32768,  -884736, -884736000,\
                               -147197952000, -262537412640768000]]
    else:
        raise NotImplementedError, "Enumerate of CM j-invariants over %s not yet implemented or tabulated"%K


def cm_j_invariants_and_orders(K):
    r"""
    Return a list of all j-invariants of CM elliptic curves such that $j \in K$.

    INPUT:
        K -- a number field (but only implemented for K=QQ)

    OUTPUT:
        list -- of CM triple (j,D,f) where j is a CM j-invariant with
        quadratic fundamental discriminant D and conductor f

    EXAMPLE:
        sage: cm_j_invariants_and_orders(QQ)
        [(-3, 3, -12288000), (-3, 2, 54000), (-3, 1, 0), (-4, 2, 287496), (-4, 1, 1728), (-7, 2, 16581375), (-7, 1, -3375), (-8, 1, 8000), (-11, 1, -32768), (-19, 1, -884736), (-43, 1, -884736000), (-67, 1, -147197952000), (-163, 1, -262537412640768000)]
    """
    if K == RationalField():
        T = [ (0,-3, 1), (54000,-3,2), (-12288000, -3,3), (1728,-1, 1), \
               (287496,-1, 2), (-3375,-7,1), (16581375, -7, 2), (8000,-2,1), \
               (-32768, -11, 1),  (-884736, -19,1), (-884736000,-43,1),\
               (-147197952000, -67,1), (-262537412640768000,-163,1)
               ]
        S = []
        for j, D, f in T:
            j = Integer(j)
            D = Integer(D)
            if not (D%4  in [0,1]):
                D *= 4
            f = Integer(f)
            S.append((D,f,j))
        S.sort()
        S.reverse()
        return S
    else:
        raise NotImplementedError, "Enumerate of CM j-invariants over %s not yet implemented or tabulated"%K

