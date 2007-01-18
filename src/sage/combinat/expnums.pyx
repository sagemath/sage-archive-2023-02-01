r"""
Compute Bell and Uppuluri-Carpenter numbers.

ALGORITHM: We use the same integer addition algorithm as GAP.  This an extension of Bell's triangle to the general case of exponential numbers.  The recursion performs $O(n^2)$ additions, but the running time is dominated by the cost of the last integer addition, because the growth of the integer results of partial computations is exponential in $n$.  The algorithm stores $O(n)$ integers, but each is exponential in $n$.

EXAMPLES:
    sage: expnums(10, 1)
    [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]

    sage: expnums(10, -1)
    [1, -1, 0, 1, 1, -2, -9, -9, 50, 267]

    sage: expnums(1, 1)
    [1]
    sage: expnums(0, 1)
    []
    sage: expnums(-1, 0)
    []

AUTHOR: Nick Alexander
"""

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"
include "../ext/gmp.pxi"

from sage.rings.integer cimport Integer

def expnums(int n, int aa):
    r"""Compute the first $n$ exponential numbers around $aa$, starting with the zero-th.

    Returns a list of length $n$.
    """
    if n < 1:
        return []

    cdef Integer z
    if n == 1:
        z = Integer.__new__(Integer)
        mpz_set_si(z.value, 1)
        return [z]

    n = n - 1

    r = []
    z = Integer.__new__(Integer)
    mpz_set_si(z.value, 1)
    r.append(z)

    z = Integer.__new__(Integer)
    mpz_set_si(z.value, aa)
    r.append(z)

    cdef mpz_t *bell
    bell = <mpz_t *>sage_malloc(sizeof(mpz_t) * (n+1))
    if bell == NULL:
        raise MemoryError, "out of memory allocating temporary storage in expnums"
    cdef mpz_t a
    mpz_init_set_si(a, aa)
    mpz_init_set_si(bell[1], aa)
    cdef int i
    cdef int k
    for i from 1 <= i < n:
        mpz_init(bell[i+1])
        mpz_mul(bell[i+1], a, bell[1])
        for k from 0 <= k < i:
            mpz_add(bell[i-k], bell[i-k], bell[i-k+1])

        z = Integer.__new__(Integer)
        mpz_set(z.value, bell[1])
        r.append(z)

    for i from 1 <= i <= n:
        mpz_clear(bell[i])
    sage_free(bell)

    return r

# The following code is from GAP 4.4.9.
###################################################
# InstallGlobalFunction(Bell,function ( n )
#     local   bell, k, i;
#     bell := [ 1 ];
#     for i  in [1..n-1]  do
#         bell[i+1] := bell[1];
#         for k  in [0..i-1]  do
#             bell[i-k] := bell[i-k] + bell[i-k+1];
#         od;
#     od;
#     return bell[1];
# end);
