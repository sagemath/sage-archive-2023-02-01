r"""
q-Numbers
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>,
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
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing

@cached_function
def q_number(n, q=None):
    r"""
    Return the `q`-number analogue of the integer `n`.

    The `q`-number of `n` is given by

    .. MATH::

        [n]_q = \frac{q^n - q^{-n}}{q - q^-1} = q^{n-1} + q^{n-3} + \cdots
        + q^{-n+3} + q^{-n+1}.

    Note that this is not the "usual" `q`-analogue of `n`, but this is
    analogue is useful in quantum groups.

    If ``q`` is unspecified, then it defaults to using the generator `q` for
    a Laurent polynomial ring over the integers. 

    INPUT:

    - ``n`` -- the value `n` defined above

    - ``q`` -- (default: ``None``) the variable `q`; if ``None``, then use a
      default variable in `\ZZ[q, q^{-1}]`

    EXAMPLES::

        sage: from sage.algebras.quantum_groups.q_numbers import q_number
        sage: q_number(2)
        q + q^-1
        sage: q_number(3)
        q^2 + 1 + q^-2
        sage: q_number(5)
        q^4 + q^2 + 1 + q^-2 + q^-4
        sage: q_number(5, 1)
        5

    TESTS::

        sage: from sage.algebras.quantum_groups.q_numbers import q_number
        sage: q_number(1)
        1
        sage: q_number(0)
        0
    """
    if q is None:
        q = LaurentPolynomialRing(ZZ, ['q']).gens()[0]
    if n == 0:
        return 0
    return sum(q**(n-2*i-1) for i in range(n))

def q_factorial(n, q=None):
    """
    Returns the `q`-analogue of the factorial `n!`

    If `q` is unspecified, then it defaults to using the generator `q` for
    a Laurent polynomial ring over the integers.

    INPUT:

    - ``n`` -- the value `n` defined above

    - ``q`` -- (default: ``None``) the variable `q`; if ``None``, then use a
      default variable in `\ZZ[q,q^{-1}]`

    EXAMPLES::

        sage: from sage.algebras.quantum_groups.q_numbers import q_factorial
        sage: q_factorial(3)
        q^3 + 2*q + 2*q^-1 + q^-3
        sage: p = LaurentPolynomialRing(QQ, 'q').gen(0)
        sage: q_factorial(3, p)
        q^3 + 2*q + 2*q^-1 + q^-3
        sage: p = ZZ['p'].gen(0)
        sage: q_factorial(3, p)
        (p^6 + 2*p^4 + 2*p^2 + 1)/p^3

    The `q`-analogue of `n!` is only defined for `n` a nonnegative
    integer (:trac:`11411`)::

        sage: q_factorial(-2)
        Traceback (most recent call last):
        ...
        ValueError: argument (-2) must be a nonnegative integer
    """
    if n in ZZ and n >= 0:
        return prod([q_number(i, q) for i in range(1, n+1)])
    raise ValueError("argument ({}) must be a nonnegative integer".format(n))

def q_binomial(n, k, q=None):
    r"""
    Return the `q`-binomial coefficient.

    This is also known as the Gaussian binomial coefficient, and is defined by

    .. MATH::

        \binom{n}{k}_q = \frac{(1-q^n)(1-q^{n-1}) \cdots (1-q^{n-k+1})}
        {(1-q)(1-q^2)\cdots (1-q^k)}.

    See :wikipedia:`Gaussian binomial coefficient`.

    If `q` is unspecified, then the variable is the generator `q` for
    a Laurent polynomial ring over the integers.

    INPUT:

    - ``n, k`` -- the values, `n` and `k` defined above

    - ``q`` -- (default: ``None``) the variable `q`; if ``None``, then use a
      default variable in `\ZZ[q,q^{-1}]`

    EXAMPLES::

        sage: from sage.algebras.quantum_groups.q_numbers import q_binomial
        sage: q_binomial(2, 1)
        q + q^-1
        sage: q_binomial(2, 0)
        q^2 + 1 + q^-2
        sage: q_binomial(4, 1)
        sage: q_binomial(4, 3)

    TESTS::

        sage: from sage.algebras.quantum_groups.q_numbers import q_binomial
        sage: all(q_binomial(n, k, 1) == binomial(n, k) for n in range(7) for k in range(n+1))
        True
        sage: q_binomial(-2, 1)
        Traceback (most recent call last):
        ...
        ValueError: n must be nonnegative
    """
    # sanity checks
    if not( n in ZZ and k in ZZ ):
        raise ValueError("arguments ({}, {}) must be integers".format(n, k))
    if n < 0:
        raise ValueError('n must be nonnegative')
    if not(0 <= k and k <= n):
        return 0

    k = min(n-k,k) # Pick the smallest k
    denomin = prod([1 - q**i for i in range(1, k+1)])
    numerat = prod([1 - q**i for i in range(n-k+1, n+1)])
    try:
        return numerat // denomin
    except TypeError:
        return numerat / denomin

