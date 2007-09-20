#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.misc.misc import prod
from sage.rings.integer_ring import IntegerRing

def int(n, p=None):
    """
    Returns the q-analogue of the integer n.

    If p is unspecified, then it defaults to using
    the generator q for a univariate polynomial
    ring over the integers.

    EXAMPLES:
        sage: q_analogues.int(3)
        q^2 + q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.int(3,p)
        p^2 + p + 1
    """
    if p == None:
        ZZ = IntegerRing()
        p = ZZ['q'].gens()[0]
        #pass
    return sum([p**i for i in range(n)])

def factorial(n, p=None):
    """
    Returns the q-analogue of the n!.

    If p is unspecified, then it defaults to using
    the generator q for a univariate polynomial
    ring over the integers.

    EXAMPLES:
        sage: q_analogues.factorial(3)
        q^3 + 2*q^2 + 2*q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.factorial(3, p)
        p^3 + 2*p^2 + 2*p + 1
    """
    return prod([int(i,p) for i in range(1, n+1)])

def binomial(n,k,p=None):
    """
    Returns the q-binomial coefficient.

    If p is unspecified, then it defaults to using
    the generator q for a univariate polynomial
    ring over the integers.

    EXAMPLES:
        sage: q_analogues.binomial(4,2)
        q^4 + q^3 + 2*q^2 + q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.binomial(4,2,p)
        p^4 + p^3 + 2*p^2 + p + 1
    """

    return factorial(n,p)/(factorial(k,p)*factorial(n-k,p))
