#  Free for any use.
#  Unfit for any purpose.
#
#  Copyright 2010, Tom Boothby

"""

This is a nearly-straightforward implementation of what
Knuth calls "Algorithm P" in TAOCP 7.2.1.2.  The intent is
to be able to enumerate permutation by "plain changes", or
multiplication by adjacent transpositions, as a generator.
This is useful when a class of objects is inherently
enumerated by permutations, but it is faster to swap items
in a permutation than construct the next object directly
from the next permutation in a list. The backtracking
algorithm in sage/graphs/genus.pyx is an example of this.

The lowest level is implemented as a struct with auxilliary
methods.  This is because Cython does not allow pointers to
class instances, so a list of these objects is inherently
slower than a list of structs.  The author prefers ugly code
to slow code.

For those willing to sacrifice a (very small) amount of
speed, we provide a class that wraps our struct.

"""

include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"
from cpython.list cimport *

##########################################################
#
# The next two functions, reset_swap and next_swap do the
# real work.  They've been implemented separately because
# the target application in sage/graphs/genus.pyx is
# cleaner and faster if it can manage its own memory.
# Both take three arguments:
#
#    int  n: number of elements to permute
#    int *c: an int array of length at least n
#    int *o: an int array of length at least o
#
# The user is advised to call reset_swap before next_swap.
#
##########################################################

cdef void reset_swap(int n, int *c, int *o):
    """
    Reset the plain_swapper to the initial state.
    """

    for i in range(n):
        c[i] = -1
        o[i] = 1

cdef int next_swap(int n, int *c, int *o):
    """
    Here's the traslation of Algorithm P.  We've modified
    it to

       a) work on zero-indexed lists
       b) operate as a generator
       c) yield the swap index rather than the current
          permutation.

    Note, Knuth's descriptions of algorithms tend to encourage
    one to think of finite state machines.  For convenience,
    we have added comments to show what state the machine is
    in at any given point in the algorithm. `plain_swap_reset`
    sets the state to 1, and this function begins and ends in
    state 2.

    Returns the index i such that the next permutation can be
    obtained by swapping P[i] <-> P[i+1]

    """

    cdef int j,s,q,offset

    #state 3
    j = n-1
    if j <= 0:
        return -1
    s = -1

    while 1:
        #state 4
        q = c[j] + o[j]
        if q == j:
            #state 6
            if j == 1:
                return -1
            s = s+1
        elif q >= -1:
            break

        #state 7
        o[j] = -o[j]
        j = j-1

    #state 5
    offset = c[j]
    if q > offset:
        offset = q
    c[j] = q
    return j - offset + s

def permutation_iterator_transposition_list(int n):
    """

    Returns a list of transposition indices to enumerate the
    permutations on `n` letters by adjacent transpositions.
    Assumes zero-based lists.  We artificially limit the
    argument to `n < 12` to avoid overflowing 32-bit pointers.
    While the algorithm works for larger `n`, the user is
    encouraged to avoid filling anything more than 4GB of
    memory with the output of this function.

    EXAMPLES:
        sage: import sage.combinat.permutation_cython
        sage: from sage.combinat.permutation_cython import permutation_iterator_transposition_list
        sage: permutation_iterator_transposition_list(4)
        [2, 1, 0, 2, 0, 1, 2, 0, 2, 1, 0, 2, 0, 1, 2, 0, 2, 1, 0, 2, 0, 1, 2]
        sage: permutation_iterator_transposition_list(200)
        Traceback (most recent call last):
        ...
        ValueError: Cowardly refusing to enumerate the permutations on more than 12 letters.
        sage: permutation_iterator_transposition_list(1)
        []

        sage: # Generate the permutations of [1,2,3,4] fixing 4.
        sage: Q = [1,2,3,4]
        sage: L = [copy(Q)]
        sage: for t in permutation_iterator_transposition_list(3):
        ...      Q[t], Q[t+1] = Q[t+1], Q[t]
        ...      L.append(copy(Q))
        sage: print L
        [[1, 2, 3, 4], [1, 3, 2, 4], [3, 1, 2, 4], [3, 2, 1, 4], [2, 3, 1, 4], [2, 1, 3, 4]]

    """

    cdef int *c, *o, N, m
    cdef list T

    if n <= 1:
        return []
    if n > 12:
        raise ValueError, "Cowardly refusing to enumerate the permutations on more than 12 letters."

    m = n-1
    N = n
    while m > 1:
        N *= m
        m -= 1

    c = <int *>sage_malloc(2*n*sizeof(int))
    if c is NULL:
        raise MemoryError, "Failed to allocate memory in permutation_iterator_transposition_list"
    o = c + n

    try:
        T = PyList_New(N-1)
    except Exception:
        sage_free(c)
        raise MemoryError, "Failed to allocate memory in permutation_iterator_transposition_list"

    reset_swap(n,c,o)
    for m in range(N-1):
        PyList_SET_ITEM(T, m, next_swap(n,c,o))
    sage_free(c)

    return T
