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

from sage.misc.cachefunc import cached_function
from sage.misc.misc import prod
from sage.rings.all import ZZ
from dyck_word import DyckWords

from partition import Partition

def q_int(n, p=None):
    r"""
    Returns the `q`-analogue of the integer `n` which is given by

    .. MATH::

            [n]_q =  \begin{cases}
            1+q+\dots+q^{n-1},  & \text{if }n\ge 0, \\
            -q^{-n} [-n]_q,     & \text{if }n\le 0.
            \end{cases}

    Consequently, if `q=1` then `[n]_1=n` and if `q\ne1` then `[n]_q=(q^n-1)/(q-1)`.

    If the argument `p` is not specified then it defaults to the generator `q`
    of the univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_int(3)
        q^2 + q + 1
        sage: q_analogues.q_int(-3)
        (-q^2 - q - 1)/q^3
        sage: p = ZZ['p'].0
        sage: q_analogues.q_int(3,p)
        p^2 + p + 1
        sage: q_analogues.q_int(3/2)
        Traceback (most recent call last):
        ...
        ValueError: 3/2 must be an integer
    """
    if not n in ZZ:
        raise ValueError, '%s must be an integer' % n

    if p == None:
        p = ZZ['q'].gens()[0]
    if n >= 0:
        return sum(p**i for i in range(n))
    else:
        return -p**n*sum(p**i for i in  range(-n))

def q_factorial(n, p=None):
    """
    Returns the `q`-analogue of the factorial `n!`.

    If `p` is unspecified, then it defaults to using the generator `q` for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_factorial(3)
        q^3 + 2*q^2 + 2*q + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_factorial(3, p)
        p^3 + 2*p^2 + 2*p + 1

    The `q`-analogue of `n!` is only defined for `n` a nonnegative
    integer (trac #11411)::

        sage: q_analogues.q_factorial(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.

    """
    if n in ZZ and n >= 0:
        return prod([q_int(i, p) for i in range(1, n+1)])
    else:
        raise ValueError("Argument (%s) must be a nonnegative integer." %n)

def q_binomial(n,k,p=None):
    """
    Returns the `q`-binomial coefficient.

    If `p` is unspecified, then it defaults to using the generator `q` for
    a univariate polynomial ring over the integers.

    The q-binomials are computed as a product of cyclotomic polynomials (cf. [CH2006]_).

    REFERENCES:

    .. [CH2006] William Y.C. Chen and Qing-Hu Hou, "Factors of the Gaussian
       coefficients", Discrete Mathematics 306 (2006), 1446-1449.
       http://dx.doi.org/10.1016/j.disc.2006.03.031

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_binomial(4,2)
        q^4 + q^3 + 2*q^2 + q + 1
        sage: q_analogues.q_binomial(4,5)
        0
        sage: p = ZZ['p'].0
        sage: q_analogues.q_binomial(4,2,p)
        p^4 + p^3 + 2*p^2 + p + 1
        sage: q_analogues.q_binomial(2,3)
        0

    The `q`-analogue of ``binomial(n,k)`` is currently only defined for
    `n` a nonnegative integer, it is zero for negative k  (trac #11411)::

        sage: q_analogues.q_binomial(5, -1)
        0
    """
    if not (n in ZZ and k in ZZ):
        raise ValueError("Argument (%s, %s) must be integers."%(n, k))
    if n < 0:
        raise NotImplementedError
    if 0 <= k and k <= n:
        if p == None:
            R = ZZ['q']
        else:
            R = ZZ[p]
        from sage.functions.all import floor
        return prod(R.cyclotomic_polynomial(d) for d in range(2,n+1) if floor(n/d)!=floor(k/d)+floor((n-k)/d))
    else:
        return 0

def q_catalan_number(n,p=None):
    """
    Returns the `q`-Catalan number of index `n`.

    If `p` is unspecified, then it defaults to using the generator `q` for
    a univariate polynomial ring over the integers.

    There are several `q`-Catalan numbers. This procedure
    returns the one which can be written using the `q`-binomial coefficients.

    EXAMPLES::

        sage: import sage.combinat.q_analogues as q_analogues
        sage: q_analogues.q_catalan_number(4)
        q^12 + q^10 + q^9 + 2*q^8 + q^7 + 2*q^6 + q^5 + 2*q^4 + q^3 + q^2 + 1
        sage: p = ZZ['p'].0
        sage: q_analogues.q_catalan_number(4,p)
        p^12 + p^10 + p^9 + 2*p^8 + p^7 + 2*p^6 + p^5 + 2*p^4 + p^3 + p^2 + 1

    The `q`-Catalan number of index `n` is only defined for `n` a
    nonnegative integer (trac #11411)::

        sage: q_analogues.q_catalan_number(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.
    """
    if n in ZZ and n >= 0:
        return prod(q_int(j, p) for j in range(n+2, 2*n+1)) / prod(q_int(j, p) for j in range(2,n+1))
    else:
        raise ValueError("Argument (%s) must be a nonnegative integer." %n)

def qt_catalan_number(n):
    """
    Returns the ``q,t``-Catalan number of index `n`.

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

    The ``q,t``-Catalan number of index `n` is only defined for `n` a
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
            tup = (dw.area(),dw.bounce())
            d[tup] = d.get(tup,0)+1
        return ZZqt(d)
    else:
        raise ValueError("Argument (%s) must be a nonnegative integer." %n)

@cached_function
def q_jordan(t, q):
    r"""
    INPUT:

    -  `t` -- a partition of an integer

    -  `q` -- an integer or an indeterminate

    OUTPUT:

    If `q` is the power of a prime number, the output is the number of
    complete flags in `F_q^N` (where `N` is the size of `t`) stable
    under a linear nilpotent endomorphism `f` whose Jordan type is
    given by `t`, i.e. such that for all `i`::

    .. MATH::

        \dim (\ker f^i) = t[0] + \cdots + t[i-1]

    If `q` is an indeterminate, the output is a polynomial whose
    values at powers of prime numbers are the previous numbers.

    The result is cached.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_jordan
        sage: [q_jordan(mu,2) for mu in Partitions(5)]
        [9765, 1029, 213, 93, 29, 9, 1]
        sage: [q_jordan(mu,2) for mu in Partitions(6)]
        [615195, 40635, 5643, 2331, 1491, 515, 147, 87, 47, 11, 1]

        sage: q=PolynomialRing(ZZ,'q').gen()
        sage: q_jordan(Partition([3,2,1]),q)
        16*q^4 + 24*q^3 + 14*q^2 + 5*q + 1

    If the partition is trivial (i.e. has only one part), we get
    the `q`-factorial (in this case, the nilpotent endomorphism is
    necessarily `0`)::

        sage: from sage.combinat.q_analogues import q_factorial
        sage: q_jordan(Partition([5]),3) == q_factorial(5,3)
        True
        sage: q_jordan(Partition([11]),5) == q_factorial(11,5)
        True

    TESTS::

        sage: q_jordan(Partition([4,3,1]),1)
        Traceback (most recent call last):
        ...
        ValueError: q must not be equal to 1

    AUTHOR:

    - Xavier Caruso (2012-06-29)
    """

    if q == 1:
        raise ValueError("q must not be equal to 1")

    if len(t) == 0:
        return 1
    tj = 0
    res = 0
    for i in range(len(t)-1,-1,-1):
        ti = t[i]
        if ti > tj:
            tp = t.to_list()
            tp[i] -= 1
            res += q_jordan(Partition(tp),q) * ((q**ti - q**tj) // (q-1))
            tj = ti
    return res
