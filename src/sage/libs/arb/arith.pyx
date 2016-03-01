"""
Arithmetic functions using the arb library
"""

#*****************************************************************************
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from ..flint.types cimport ulong
from ..flint.fmpq cimport fmpq_t, fmpq_init, fmpq_clear, fmpq_get_mpq
from .bernoulli cimport bernoulli_fmpq_ui
from sage.rings.rational cimport Rational

def bernoulli(n):
    """
    Return the ``n``-th Bernoulli number using ``arb``.

    INPUT:

    - ``n`` -- unsigned integer

    OUTPUT: a rational number

    EXAMPLES::

        sage: from sage.libs.arb.arith import bernoulli
        sage: bernoulli(24)
        -236364091/2730
        sage: bernoulli(0)
        1
        sage: bernoulli(1)
        -1/2

    TESTS::

        sage: bernoulli(-1)
        Traceback (most recent call last):
        ...
        OverflowError: can't convert negative value to mp_limb_t
    """
    cdef ulong i = n
    cdef Rational q = <Rational>Rational.__new__(Rational)
    cdef fmpq_t x
    fmpq_init(x)
    bernoulli_fmpq_ui(x, i)
    fmpq_get_mpq(q.value, x)
    fmpq_clear(x)
    return q
