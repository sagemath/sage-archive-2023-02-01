"""
Number of partitions of integer

AUTHOR:
    -- William Stein (2007-07-28): initial version
    -- Jonathon Bober (2007-07-28): wrote the program partitions_c.cc
                  that does all the actual heavy lifting.
"""

import sys

cdef extern from "gmp.h":
    ctypedef void* mpz_t

cdef extern from "partitions_c.h":
    int part(mpz_t answer, unsigned int n)

include "../ext/interrupt.pxi"

from sage.rings.integer cimport Integer

def number_of_partitions(n):
    cdef unsigned int nn = n
    if n != nn:
        raise ValueError, "input must be a nonnegative integer less than %s"%(2*sys.maxint)

    if nn <= 0:
        return Integer(0)
    elif nn == 1:
        return Integer(1)  # part hangs on n=1 as input.

    cdef Integer ans = Integer(0)

    _sig_on
    part(ans.value, nn)
    _sig_off

    return ans
