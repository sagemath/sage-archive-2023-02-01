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
    Returns the ``q``-analogue of the integer ``n``.

    If ``p`` is unspecified, then it defaults to using the generator ``q`` for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_int(3)
        q^2 + q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_int(3,p)
        p^2 + p + 1

    The ``q``-analogue of ``n`` is only defined for ``n`` a nonnegative
    integer (trac #11411)::

        sage: q_analogues.q_int(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.
    """
    if n in ZZ and n >= 0:
        if p == None:
            p = ZZ['q'].gens()[0]
        return sum([p**i for i in range(n)])
    else:
        raise ValueError, "Argument (%s) must be a nonnegative integer." %n

def q_factorial(n, p=None):
    """
    Returns the ``q``-analogue of the factorial ``n!``.

    If ``p`` is unspecified, then it defaults to using the generator ``q`` for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_factorial(3)
        q^3 + 2*q^2 + 2*q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_factorial(3, p)
        p^3 + 2*p^2 + 2*p + 1

    The ``q``-analogue of ``n!`` is only defined for ``n`` a nonnegative
    integer (trac #11411)::

        sage: q_analogues.q_factorial(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.
    """
    if n in ZZ and n >= 0:
        return prod([q_int(i, p) for i in range(1, n+1)])
    else:
        raise ValueError, "Argument (%s) must be a nonnegative integer." %n

def q_binomial(n,k,p=None):
    """
    Returns the ``q``-binomial coefficient.

    If ``p`` is unspecified, then it defaults to using the generator ``q`` for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_binomial(4,2)
        q^4 + q^3 + 2*q^2 + q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_binomial(4,2,p)
        p^4 + p^3 + 2*p^2 + p + 1

    The ``q``-analogue of ``binomial(n,k)`` is currently only defined for
    ``n`` a nonnegative integer, it is zero for negative k  (trac #11411)::

        sage: q_analogues.q_binomial(5, -1)
        0
    """
    if not (n in ZZ and k in ZZ):
        raise ValueError, "Argument (%s, %s) must be integers."%(n, k)
    if n < 0:
        raise NotImplementedError
    if 0 <= k and k <= n:
        k=min(k, n-k)
        return (prod(q_int(j, p) for j in range(n-k+1, n+1)) /
                prod(q_int(j, p) for j in range(1, k+1)))
    else:
        return 0

def q_catalan_number(n,p=None):
    """
    Returns the ``q``-Catalan number of index ``n``.

    If ``p`` is unspecified, then it defaults to using the generator ``q`` for
    a univariate polynomial ring over the integers.

    There are several ``q``-Catalan numbers. This procedure
    returns the one which can be written using the ``q``-binomial coefficients.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_catalan_number(4)
        q^12 + q^10 + q^9 + 2*q^8 + q^7 + 2*q^6 + q^5 + 2*q^4 + q^3 + q^2 + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_catalan_number(4,p)
        p^12 + p^10 + p^9 + 2*p^8 + p^7 + 2*p^6 + p^5 + 2*p^4 + p^3 + p^2 + 1

    The ``q``-Catalan number of index ``n`` is only defined for ``n`` a
    nonnegative integer (trac #11411)::

        sage: q_analogues.q_catalan_number(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.
    """
    if n in ZZ and n >= 0:
        return prod(q_int(j, p) for j in range(n+2, 2*n+1)) / prod(q_int(j, p) for j in range(2,n+1))
    else:
        raise ValueError, "Argument (%s) must be a nonnegative integer." %n

def qt_catalan_number(n):
    """
    Returns the ``q,t``-Catalan number of index ``n``.

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

    The ``q,t``-Catalan number of index ``n`` is only defined for ``n`` a
    nonnegative integer (trac #11411)::

        sage: q_analogues.qt_catalan_number(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.
    """
    if n in ZZ and n >= 0:
        ZZqt = ZZ['q','t']
        d = {}
        for dw in DyckWords(n):
            tup = (dw.a_statistic(),dw.b_statistic())
            d[tup] = d.get(tup,0)+1
        return ZZqt(d)
    else:
        raise ValueError, "Argument (%s) must be a nonnegative integer." %n

