"""
TESTS:

We verify that Trac #6919 is correctly fixed::

    sage: R.<x> = PolynomialRing(ZZ)
    sage: A = 2^(2^17+2^15)
    sage: a = A * x^31
    sage: b = (A * x) * x^30
    sage: a == b
    True
"""
include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"

def free_flint_stack():
    flint_stack_cleanup()
