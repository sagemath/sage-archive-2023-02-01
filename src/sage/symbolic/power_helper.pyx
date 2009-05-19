import ring

cdef bint __pynac_pow = False

cpdef object try_symbolic_power(object obj, object p):
    """
    While evaluating powers, pynac calls Sage again to check if a fractional
    power is an exact root. This leads to infinite recursions and crashes.
    We use this hack to prevent an infinite recursion.

    TESTS::

        sage: 2^(1/2)
        sqrt(2)
        sage: (2*I)^(1/2)
        sqrt(2*I)

    """
    global __pynac_pow
    if __pynac_pow:
        __pynac_pow = False
        return None
    else:
        try:
            __pynac_pow = True
            return ring.SR(obj)**ring.SR(p)
        finally:
            __pynac_pow = False
