"""
Implementation of the Berlekamp-Massey Algorithm
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

import sage.rings.polynomial_ring as polynomial_ring

def berlekamp_massey(a):
    """
    Use the Berlekamp-Massey algorithm to find the minimal polynomial
    of a linearly generated sequence a.

    The minimal polynomial of a linear recurrence $\{a_r\}$ is by
    definition the unique monic polynomial~$g$, such that if $\{a_r\}$
    satisfies a linear recurrence $a_{j+k} + b_{j-1} a_{j-1+k} +
    \cdots + b_0 a_k=0$ (for all $k\geq 0$), then $g$ divides the
    polynomial $x^j + \sum_{i=0}^{j-1} b_i x^i$.

    INPUT:
        a -- a list of even length of elements of a single field

    OUTPUT:
        Polynomial -- the minimal polynomial of the sequence

    EXAMPLES:

    """

    if not isinstance(a, list):
        raise TypeError, "Argument 1 must be a list."
    if len(a)%2 != 0:
        raise ValueError, "Argument 1 must have an even number of terms."

    M = len(a)//2

    K = a[0].parent()
    R = polynomial_ring.PolynomialRing(K)
    x = R.gen()

    f = {}
    q = {}
    s = {}
    t = {}
    f[-1] = R(a)
    f[0] = x**(2*M)
    s[-1] = 1
    t[0] = 1
    s[0] = 0
    t[-1] = 0
    j = 0
    while f[j].degree() >= M:
        j += 1
        q[j], f[j] = f[j-2].quo_rem(f[j-1])
        # assert q[j]*f[j-1] + f[j] == f[j-2], "poly divide failed."
        s[j] = s[j-2] - q[j]*s[j-1]
        t[j] = t[j-2] - q[j]*t[j-1]
    t = s[j].reverse()
    return ~(t[t.degree()]) * t  # make monic  (~ is inverse in python)

