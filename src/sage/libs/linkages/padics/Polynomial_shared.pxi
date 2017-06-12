"""
This linkage file implements the padics API using Sage Polynomials.

It is included into Polynomial_ram.pxi and Polynomial_unram.pxi,
where functions that depend on the ramification of the defining polynomial are placed.

AUTHORS:

- David Roe (2017-6-11) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.list cimport *
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.all import ZZ, QQ
from sage.ext.stdsage cimport PY_NEW
from copy import copy

cdef inline int cconstruct(celement value, PowComputer_ prime_pow) except -1:
    """
    Construct a new element.

    INPUT:

    - ``unit`` -- an ``celement`` to be initialized.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    if value is None:
        raise TypeError
    value.__init__(prime_pow.poly_ring)

cdef inline int cdestruct(celement value, PowComputer_ prime_pow) except -1:
    """
    Deallocate an element.

    INPUT:

    - ``unit`` -- an ``celement`` to be cleared.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    pass

cdef inline int cneg(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    """
    Negation

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the negation.
    - ``a`` -- an ``celement`` to be negated.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = (<celement?>(-a)).__coeffs

cdef inline int cadd(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    """
    Addition

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the sum.
    - ``a`` -- an ``celement``, the first summand.
    - ``b`` -- an ``celement``, the second summand.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = (<celement?>(a + b)).__coeffs

cdef inline int csub(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    """
    Subtraction.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the difference.
    - ``a`` -- an ``celement``, the first input.
    - ``b`` -- an ``celement``, the second input.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = (<celement?>(a - b)).__coeffs

cdef inline int cmul(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    """
    Multiplication.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the product.
    - ``a`` -- an ``celement``, the first input.
    - ``b`` -- an ``celement``, the second input.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = (<celement?>(a*b)).__coeffs

cdef inline int csetone(celement out, PowComputer_ prime_pow) except -1:
    """
    Sets to 1.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 1.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = [prime_pow.base_ring(1)]
    pass

cdef inline int csetzero(celement out, PowComputer_ prime_pow) except -1:
    """
    Sets to 0.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 0.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = []

cdef inline bint cisone(celement a, PowComputer_ prime_pow) except -1:
    """
    Returns whether this element is equal to 1.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a = 1`, and False otherwise.
    """
    # Can do this more efficiently
    return a == 1

cdef inline bint ciszero(celement a, PowComputer_ prime_pow) except -1:
    """
    Returns whether this element is equal to 0.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a = 0`, and False otherwise.
    """
    # Can do this more efficiently
    return a == 0

cdef inline int ccopy(celement out, celement a, PowComputer_ prime_pow) except -1:
    """
    Copying.

    INPUT:

    - ``out`` -- the ``celement`` to store the result.
    - ``a`` -- the element to copy.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = copy((<celement?>a).__coeffs)

cdef inline cpickle(celement a, PowComputer_ prime_pow):
    """
    Serialization into objects that Sage knows how to pickle.

    INPUT:

    - ``a`` the element to pickle.
    - ``prime_pow`` the PowComputer for the ring.

    OUTPUT:

    - a serializable object storing ``a``.
    """
    return a.__coeffs

cdef inline int cunpickle(celement out, x, PowComputer_ prime_pow) except -1:
    """
    Reconstruction from the output of meth:`cpickle`.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result.
    - ``x`` -- the result of `meth`:cpickle.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    out.__coeffs = x

cdef inline long chash(celement a, long ordp, long prec, PowComputer_ prime_pow) except -1:
    """
    Hashing.

    INPUT:

    - ``a`` -- an ``celement`` storing the underlying element to hash.
    - ``ordp`` -- a long storing the valuation.
    - ``prec`` -- a long storing the precision.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    pass

# The element is filled in for zero in the output of clist if necessary.
# It could be [] for some other linkages.
_list_zero = Integer(0)

cdef int cconv(celement out, x, long prec, long valshift, PowComputer_ prime_pow) except -2:
    """
    Conversion from other Sage types.

    INPUT:

    - ``out`` -- an ``celement`` to store the output.
    - ``x`` -- a Sage element that can be converted to a `p`-adic element.
    - ``prec`` -- a long, giving the precision desired: absolute if
                  `valshift = 0`, relative if `valshift > 0`.
    - ``valshift`` -- the power of the uniformizer to divide by before
      storing the result in ``out``.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    if out is None:
        raise TypeError
    # This needs to be improved
    out.__coeffs = (<celement?>(prime_pow.poly_ring(x) % prime_pow.modulus)).__coeffs

cdef inline long cconv_mpz_t(celement out, mpz_t x, long prec, bint absolute, PowComputer_ prime_pow) except -2:
    """
    A fast pathway for conversion of integers that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``celement`` to store the output.
    - ``x`` -- an ``mpz_t`` giving the integer to be converted.
    - ``prec`` -- a long, giving the precision desired: absolute or
                  relative depending on the ``absolute`` input.
    - ``absolute`` -- if False then extracts the valuation and returns
                      it, storing the unit in ``out``; if True then
                      just reduces ``x`` modulo the precision.
    - ``prime_pow`` -- a PowComputer for the ring.

    OUTPUT:

    - If ``absolute`` is False then returns the valuation that was
      extracted (``maxordp`` when `x = 0`).
    """
    cdef Integer n = PY_NEW(Integer)
    mpz_set(n.value, x)
    out.__coeffs = [prime_pow.base_ring(n)]
    out.__normalize()

cdef inline int cconv_mpz_t_out(mpz_t out, celement x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    """
    Converts the underlying `p`-adic element into an integer if
    possible.

    - ``out`` -- stores the resulting integer as an integer between 0
                 and `p^{prec + valshift}`.
    - ``x`` -- an ``celement`` giving the underlying `p`-adic element.
    - ``valshift`` -- a long giving the power of `p` to shift `x` by.
    -` ``prec`` -- a long, the precision of ``x``: currently not used.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    cdef Integer n
    if len(x.__coeffs) > 1:
        raise ValueError("Cannot convert to integer")
    elif len(x.__coeffs) == 0:
        mpz_set_ui(out, 0)
    else:
        n = ZZ(x.__coeffs[0])
        mpz_set(out, n.value)

cdef inline long cconv_mpq_t(celement out, mpq_t x, long prec, bint absolute, PowComputer_ prime_pow) except? -10000:
    """
    A fast pathway for conversion of rationals that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``celement`` to store the output.
    - ``x`` -- an ``mpq_t`` giving the rational to be converted.
    - ``prec`` -- a long, giving the precision desired: absolute or
                  relative depending on the ``absolute`` input.
    - ``absolute`` -- if False then extracts the valuation and returns
                      it, storing the unit in ``out``; if True then
                      just reduces ``x`` modulo the precision.
    - ``prime_pow`` -- a PowComputer for the ring.

    OUTPUT:

    - If ``absolute`` is False then returns the valuation that was
      extracted (``maxordp`` when `x = 0`).
    """
    cdef Rational c = Rational.__new__(Rational)
    mpq_set(c.value, x)
    out.__coeffs = [prime_pow.base_ring(c)]
    out.__normalize()

cdef inline int cconv_mpq_t_out(mpq_t out, celement x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    """
    Converts the underlying `p`-adic element into a rational.

    - ``out`` -- gives a rational approximating the input.  Currently
                 uses rational reconstruction but may change in the
                 future to use a more naive method.
    - ``x`` -- an ``celement`` giving the underlying `p`-adic element.
    - ``valshift`` -- a long giving the power of `p` to shift `x` by.
    -` ``prec`` -- a long, the precision of ``x``, used in rational
                   reconstruction.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    cdef Rational c
    if len(x.__coeffs) > 1:
        raise ValueError("Cannot convert to integer")
    elif len(x.__coeffs) == 0:
        mpq_set_ui(out, 0, 1)
    else:
        c = QQ(x.__coeffs[0])
        mpq_set(out, c.value)
