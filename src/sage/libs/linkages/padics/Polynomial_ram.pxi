"""

"""

cdef inline int ccmp(celement a, celement b, long prec, bint reduce_a, bint reduce_b, PowComputer_class prime_pow) except -2:
    """
    Comparison of two elements.

    INPUT:

    - ``a`` -- an ``celement``.
    - ``b`` -- an ``celement``.
    - ``prec`` -- a long, the precision of the comparison.
    - ``reduce_a`` -- a bint, whether ``a`` needs to be reduced.
    - ``reduce_b`` -- a bint, whether ``b`` needs to be reduced.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - If neither ``a`` nor ``b`` needs to be reduced, returns
      -1 (if `a < b`), 0 (if `a == b`) or 1 (if `a > b`)

    - If at least one needs to be reduced, returns
      0 (if ``a == b mod p^prec``) or 1 (otherwise)
    """
    pass

cdef inline bint creduce(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Reduce modulo a power of the maximal ideal.

    INPUT:

    - ``out`` -- an ``celement`` to store the reduction.
    - ``a`` -- the element to be reduced.
    - ``prec`` -- a long, the precision to reduce modulo.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if the reduction is zero; False otherwise.
    """
    pass

cdef inline bint creduce_small(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Reduce modulo a power of the maximal ideal.

    This function assumes that at most one addition/subtraction has
    happened on reduced inputs.  For integral inputs this translates
    to the assumption that `-p^prec < a < 2p^prec`.

    INPUT:

    - ``out`` -- an ``celement`` to store the reduction.
    - ``a`` -- the element to be reduced.
    - ``prec`` -- a long, the precision to reduce modulo.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if the reduction is zero; False otherwise.
    """
    pass

cdef inline long cremove(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Extract the maximum power of the uniformizer dividing this element.

    INPUT:

    - ``out`` -- an ``celement`` to store the unit.
    - ``a`` -- the element whose valuation and unit are desired.
    - ``prec`` -- a long, used if `a = 0`.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - if `a = 0`, returns prec (the value of ``out`` is undefined).
      Otherwise, returns the number of times `p` divides `a`.
    """
    pass

cdef inline long cvaluation(celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Returns the maximum power of the uniformizer dividing this
    element.

    This function differs from :meth:`cremove` in that the unit is
    discarded.

    INPUT:

    - ``a`` -- the element whose valuation is desired.
    - ``prec`` -- a long, used if `a = 0`.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - if `a = 0`, returns prec.  Otherwise, returns the number of
      times p divides a.
    """
    pass

cdef inline bint cisunit(celement a, PowComputer_class prime_pow) except -1:
    """
    Returns whether this element has valuation zero.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a` has valuation 0, and False otherwise.
    """
    pass

cdef inline int cshift(celement out, celement a, long n, long prec, PowComputer_class prime_pow, bint reduce_afterward) except -1:
    """
    Multiplies by a power of the uniformizer.

    INPUT:

    - ``out`` -- an ``celement`` to store the result.  If `n >= 0`
      then out will be set to `a * p^n`.
      If `n < 0`, out will be set to `a // p^-n`.
    - ``a`` -- the element to shift.
    - ``n`` -- long, the amount to shift by.
    - ``prec`` -- long, a precision modulo which to reduce.
    - ``prime_pow`` -- the PowComputer for the ring.
    - ``reduce_afterward`` -- whether to reduce afterward.
    """
    pass

cdef inline int cshift_notrunc(celement out, celement a, long n, long prec, PowComputer_class prime_pow) except -1:
    """
    Multiplies by a power of the uniformizer, assuming that the
    valuation of a is at least -n.

    INPUT:

    - ``out`` -- an ``celement`` to store the result.  If `n >= 0`
      then out will be set to `a * p^n`.
      If `n < 0`, out will be set to `a // p^-n`.
    - ``a`` -- the element to shift.  Assumes that the valuation of a
      is at least -n.
    - ``n`` -- long, the amount to shift by.
    - ``prec`` -- long, a precision modulo which to reduce.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    pass

cdef inline int cinvert(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Inversion

    The result will be reduced modulo p^prec.

    INPUT:

    - ``out`` -- an ``celement`` to store the inverse.
    - ``a`` -- an ``celement``, the element to be inverted.
    - ``prec`` -- a long, the precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    pass

cdef inline int cdivunit(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
    """
    Division.

    The inversion is performed modulo p^prec.  Note that no reduction
    is performed after the product.

    INPUT:

    - ``out`` -- an ``celement`` to store the quotient.
    - ``a`` -- an ``celement``, the first input.
    - ``b`` -- an ``celement``, the second input.
    - ``prec`` -- a long, the precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    pass

cdef inline int cpow(celement out, celement a, mpz_t n, long prec, PowComputer_class prime_pow) except -1:
    """
    Exponentiation.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result.
    - ``a`` -- the base.
    - ``n`` -- an ``mpz_t``, the exponent.
    - ``prec`` -- a long, the working absolute precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    pass

cdef clist(celement a, long prec, bint pos, PowComputer_class prime_pow):
    """
    Returns a list of digits in the series expansion.

    This function is used in printing, and expresses ``a`` as a series
    in the standard uniformizer ``p``.  Note that for extension
    elements, "digits" could be lists themselves.

    INPUT:

    - ``a`` -- an ``celement`` giving the underlying `p`-adic element.
    - ``prec`` -- a precision giving the number of digits desired.
    - ``pos`` -- if True then representatives in 0..(p-1) are used;
                 otherwise the range (-p/2..p/2) is used.
    - ``prime_pow`` -- a PowComputer for the ring.

    OUTPUT:

    - A list of p-adic digits `[a_0, a_1, \ldots]` so that `a = a_0 + a_1*p + \cdots` modulo `p^{prec}`.
    """
    pass

