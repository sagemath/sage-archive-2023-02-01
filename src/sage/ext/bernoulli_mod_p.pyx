r"""
Bernoulli numbers modulo p

AUTHOR:
    - David Harvey (2006-07-26): initial version
    - William Stein (2006-07-28): some touch up.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Harvey <dmharvey@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cimport arith
import arith
cdef arith.arith_int arith_int
arith_int  = arith.arith_int()

ctypedef long long llong

from sage.libs.ntl import all as ntl
import sage.libs.pari.gen


def bernoulli_mod_p(int p, scale_properly=True):
    r"""
    Computes bernoulli numbers $B_0, B_2, ... B_{p-3}$ modulo $p$.

    Returns a list of integers [B(0), B(2), ... B(p-3)].

    If scale_properly is False, then it will save a little time by skipping
    a final scaling step; the kth entry of the result (for k > 0) will be
        $$ \frac{1-2k}{(2k)!} B_{2k} $$
    instead of $B_{2k}$ itself. This is useful if you only care about whether
    $B_{2k}$ is zero mod $p$ (for computing irregular primes etc).

    ALGORITHM:
        Use the series identity:

            $$ \frac{x^2}{\cosh x - 1} =
                    -2 + \sum_{n \geq 0} \frac{2n-1}{(2n)!} B_{2n} x^{2n} $$

        (see "Irregular Primes to One Million", Buhler et al.)

        Use NTL to invert the power series mod p.

    PERFORMANCE:
        Should be complexity $O(p \log p)$. Can handle p around $10^5$ in
        5 seconds on my desktop machine.

    INPUT:
        p -- integer, a prime
        scale_properly -- bool (Default: True)

    OUTPUT:
        list -- the bernoulli numbers

    EXAMPLES:
    Check the results against PARI's C-library implemention (that
    computes exact rationals) for $p = 37$:

        sage: bernoulli_mod_p(37)
         [1, 31, 16, 15, 16, 4, 17, 32, 22, 31, 15, 15, 17, 12, 29, 2, 0, 2]
        sage: [bernoulli(n) % 37 for n in xrange(0, 36, 2)]
         [1, 31, 16, 15, 16, 4, 17, 32, 22, 31, 15, 15, 17, 12, 29, 2, 0, 2]

    An example with scale_properly = False; note that the second-last entry
    is still 0, corresponding to $B_{32} = 0$ mod 37.

        sage: bernoulli_mod_p(37, False)
         [1, 3, 35, 13, 26, 3, 5, 11, 32, 20, 14, 8, 36, 26, 14, 24, 0, 21]
        sage: [(1-2*k)/(factorial(2*k))*bernoulli(2*k) % 37 for k in xrange(18)]
         [1, 3, 35, 13, 26, 3, 5, 11, 32, 20, 14, 8, 36, 26, 14, 24, 0, 21]

    Boundary case:
        sage: bernoulli_mod_p(3)
         [1]

    AUTHOR: David Harvey
    """
    # todo: check p is >= 3 and not too big (should fit into regular
    # int I guess). Not sure how I should be reporting errors here...

    if p <= 2:
        raise ValueError, "p (=%s) must be a prime >= 3"%p

    if not sage.libs.pari.gen.pari(p).isprime():
        raise ValueError, "p (=%s) must be a prime"%p

    # (The conversion of p to an int is automatically by Pyrex since p
    # is declared to be an int.)

    # todo: The article by Buhler et al. explains how to use "multisectioning"
    # to reduce the problem to inverting a series of length p/8, instead of
    # length p/2 as we do here; it would be faster that way.

    N = (p-1)/2     # length of polynomials

    ntl.set_modulus(ntl.ZZ(p))

    # Compute the series that we want to invert.

    # We normalise the series (multiply by 2) so that the constant
    # term is 1 (otherwise NTL gets angry).

    series = ntl.ZZ_pX()
    series.preallocate_space(N)
    cdef int product, factor1, factor2
    cdef int i
    product = p-2
    for i from N-1 >= i >= 1:
        series[i] = product
        # multiply product by (2*i+1)*(2*i+2) mod p
        factor1 = ((i << 1) + 1)
        factor2 = ((i << 1) + 2)
        product = <int> (((((<llong> product) * factor1) % p) * factor2) % p)
    series[0] = product

    # invert the series
    series = series.invert_and_truncate(N)

    # extract the result
    cdef int value
    if scale_properly:
        output = [1]
        product = -1
        for i from 1 <= i <= N-1:
            # multiply product by (2*i)*(2*i-1) mod p
            factor1 = <llong> ((i << 1) - 1)
            factor2 = <llong> (i << 1)
            product = <int> (((((<llong> product) * factor1) % p) * factor2) % p)
            value = ((<llong> product) * series[i] * arith_int.c_inverse_mod_int(2*i-1, p)) % p
            output.append(value)
    else:
        output = [1]
        for i from 1 <= i <= N-1:
            output.append(series[i])

    return output



def verify_bernoulli_mod_p(data):
    """
    Computes checksum for bernoulli numbers.

    It checks the identity
        $$ \sum_{n=0}^{(p-3)/2} 2^{2n} (2n+1) B_{2n}  \equiv  -2 \pmod p $$

    (see "Irregular Primes to One Million", Buhler et al)

    INPUT:
        data -- list, same format as output of bernoulli_mod_p function

    OUTPUT:
        bool -- True if checksum passed

    EXAMPLES:
        sage: from sage.ext.bernoulli_mod_p import verify_bernoulli_mod_p
        sage: verify_bernoulli_mod_p(bernoulli_mod_p(next_prime(3)))
        True
        sage: verify_bernoulli_mod_p(bernoulli_mod_p(next_prime(1000)))
        True
        sage: verify_bernoulli_mod_p([1, 2, 4, 5, 4])
        True
        sage: verify_bernoulli_mod_p([1, 2, 3, 4, 5])
        False

    This one should test that long longs are working:
        sage: verify_bernoulli_mod_p(bernoulli_mod_p(next_prime(20000)))
        True

    AUTHOR: David Harvey
    """

    cdef int N, p, i, product, sum, value, element
    N = len(data)
    p = N*2 + 1
    product = 1
    sum = 0
    for i from 0 <= i < N:
        element = data[i]
        value = <int> (((((<llong> product) * (2*i+1)) % p) * element) % p)
        sum = (sum + value) % p
        product = (4 * product) % p

    if (sum + 2) % p == 0:
        return True
    else:
        return False


# ============ end of file
