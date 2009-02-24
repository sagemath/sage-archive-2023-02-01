r"""
q-Analogues
"""
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
from sage.rings.all import ZZ
from dyck_word import DyckWords

def q_int(n, p=None):
    """
    Returns the q-analogue of the integer n.

    If p is unspecified, then it defaults to using the generator q for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_int(3)
        q^2 + q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_int(3,p)
        p^2 + p + 1
    """
    if p == None:
        p = ZZ['q'].gens()[0]
        #pass
    return sum([p**i for i in range(n)])

def q_factorial(n, p=None):
    """
    Returns the q-analogue of the n!.

    If p is unspecified, then it defaults to using the generator q for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_factorial(3)
        q^3 + 2*q^2 + 2*q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_factorial(3, p)
        p^3 + 2*p^2 + 2*p + 1
    """
    return prod([q_int(i,p) for i in range(1, n+1)])

def q_binomial(n,k,p=None):
    """
    Returns the q-binomial coefficient.

    If p is unspecified, then it defaults to using the generator q for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_binomial(4,2)
        q^4 + q^3 + 2*q^2 + q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_binomial(4,2,p)
        p^4 + p^3 + 2*p^2 + p + 1
    """

    return q_factorial(n,p)/(q_factorial(k,p)*q_factorial(n-k,p))


def qt_catalan_number(n):
    """
    Returns the q,t-Catalan number.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.qt_catalan_number(1)
        1
        sage: q_analogues.qt_catalan_number(2)
        q + t
        sage: q_analogues.qt_catalan_number(3)
        q^3 + q^2*t + q*t^2 + t^3 + q*t
        sage: q_analogues.qt_catalan_number(4)
        q^6 + q^5*t + q^4*t^2 + q^3*t^3 + q^2*t^4 + q*t^5 + t^6 + q^4*t + q^3*t^2 + q^2*t^3 + q*t^4 + q^3*t + q^2*t^2 + q*t^3
    """
    ZZqt = ZZ['q','t']

    d = {}
    for dw in DyckWords(n):
        tup = (dw.a_statistic(),dw.b_statistic())
        d[tup] = d.get(tup,0)+1
    return ZZqt(d)
