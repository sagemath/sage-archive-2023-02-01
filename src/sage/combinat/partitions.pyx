"""
Number of partitions of integer

AUTHOR:
    -- William Stein (2007-07-28): initial version
    -- Jonathan Bober (2007-07-28): wrote the program partitions_c.cc
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
    n = Integer(n)
    if n <= 0:
        return Integer(0)
    elif n == 1:
        return Integer(1)  # part hangs on n=1 as input.
    if n >= Integer('4294967296'):
        raise ValueError, "input must be a nonnegative integer less than 4294967296."
    cdef unsigned int nn = n

    cdef Integer ans = Integer(0)

    _sig_on
    part(ans.value, nn)
    _sig_off

    return ans
