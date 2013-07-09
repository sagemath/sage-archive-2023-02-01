"""
"""
# COPYRIGHT

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"

include "sage/flint/fmpz_poly.pxi"

from sage.rings.padics.pow_computer_flint cimport PowComputer_flint_unram

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
import sage.rings.finite_rings.integer_mod

cdef inline int cconstruct(celement value, PowComputer_class prime_pow) except -1:
    """
    Construct a new element.

    INPUT:

    - ``unit`` -- an ``celement`` to be initialized.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    fmpz_poly_init(value)

cdef inline int cdestruct(celement value, PowComputer_class prime_pow) except -1:
    """
    Deallocate an element.

    INPUT:

    - ``unit`` -- an ``celement`` to be cleared.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    fmpz_poly_clear(value)

cdef inline int ccmp(celement a, celement b, long prec, bint reduce_a, bint reduce_b, PowComputer_class prime_pow_) except -2:
    """
    Comparison of two elements.

    INPUT:

    - ``a`` -- an ``celement``.
    - ``b`` -- an ``celement``.
    - ``prec`` -- a long, the precision of the comparison.
    - ``reduce_a`` -- a bint, whether ``a`` needs to be reduced.
    - ``reduce_b`` -- a bint, whether ``b`` needs to be reduced.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUPUT:

    - If neither ``a`` nor ``b`` needs to be reduced, returns
      -1 (if `a < b`), 0 (if `a == b`) or 1 (if `a > b`)

    - If at least one needs to be reduced, returns
      0 (if ``a == b mod p^prec``) or 1 (otherwise)
    """
    cdef PowComputer_flint_unram prime_pow = <PowComputer_flint_unram>prime_pow_
    cdef long cmp

    if reduce_a or reduce_b:
        csub(prime_pow.tmp_poly, a, b, prec, prime_pow)
        creduce(prime_pow.tmp_poly, prime_pow.tmp_poly, prec, prime_pow)
        return not ciszero(prime_pow.tmp_poly, prime_pow)

    cdef long da = fmpz_poly_degree(a)
    cdef long db = fmpz_poly_degree(b)
    if da < db: return -1
    elif da > db: return 1
    for i from 0 <= i <= da:
        fmpz_poly_get_coeff_fmpz(prime_pow.ftmp, a, i)
        fmpz_poly_get_coeff_fmpz(prime_pow.ftmp2, b, i)
        cmp = fmpz_cmp(prime_pow.ftmp, prime_pow.ftmp2)
        if cmp < 0: return -1
        elif cmp > 0: return 1
    return 0

cdef inline int cneg(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Negation

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the negation.
    - ``a`` -- an ``celement`` to be negated.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    fmpz_poly_neg(out, a)

cdef inline int cadd(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
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
    fmpz_poly_add(out, a, b)

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

    fmpz_poly_scalar_mod_fmpz(out, a, (<PowComputer_flint_unram>prime_pow).pow_fmpz_t_tmp(prec)[0])
    return ciszero(out, prime_pow)

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
    return creduce(out, a, prec, prime_pow)

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
    if ciszero(a, prime_pow):
        return prec
    cdef long ret = cvaluation(a, prec, prime_pow)
    fmpz_poly_scalar_divexact_fmpz(out, a, (<PowComputer_flint_unram>prime_pow).pow_fmpz_t_tmp(ret)[0])
    return ret

cdef inline long cvaluation(celement a, long prec, PowComputer_class prime_pow_) except -1:
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
    cdef PowComputer_flint_unram prime_pow = <PowComputer_flint_unram>prime_pow_
    if ciszero(a, prime_pow):
        return prec
    cdef long ret = maxordp
    cdef long val
    fmpz_set_si(prime_pow.ftmp, prec)
    for i from 0 <= i <= fmpz_poly_degree(a):
        fmpz_poly_get_coeff_fmpz(prime_pow.ftmp2, a, i)
        val = fmpz_remove(prime_pow.ftmp2, prime_pow.ftmp2, prime_pow.ftmp)
        if val < ret: ret = val
    return ret

cdef inline bint cisunit(celement a, PowComputer_class prime_pow_) except -1:
    """
    Returns whether this element has valuation zero.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a` has valuation 0, and False otherwise.
    """
    cdef PowComputer_flint_unram prime_pow = <PowComputer_flint_unram>prime_pow_
    fmpz_poly_scalar_mod_fmpz(prime_pow.tmp_poly, a, prime_pow.fprime)
    return ciszero(prime_pow.tmp_poly, prime_pow)

cdef inline int cshift(celement out, celement a, long n, long prec, PowComputer_class prime_pow_, bint reduce_afterward) except -1:
    """
    Mulitplies by a power of the uniformizer.

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
    cdef PowComputer_flint_unram prime_pow = <PowComputer_flint_unram>prime_pow_
    if n > 0:
        fmpz_poly_scalar_mul_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(n)[0])
    elif n < 0:
        sig_on()
        fmpz_poly_scalar_fdiv_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(-n)[0])
        sig_off()
    else:
        fmpz_poly_set(out, a)
    if reduce_afterward:
        creduce(out, out, prec, prime_pow)

cdef inline int cshift_notrunc(celement out, celement a, long n, long prec, PowComputer_class prime_pow_) except -1:
    """
    Mulitplies by a power of the uniformizer, assuming that the
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
    cdef PowComputer_flint_unram prime_pow = <PowComputer_flint_unram>prime_pow_
    if n > 0:
        fmpz_poly_scalar_mul_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(n)[0])
    elif n < 0:
        sig_on()
        fmpz_poly_scalar_divexact_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(-n)[0]) # ??
        sig_off()
    else:
        fmpz_poly_set(out, a)

cdef inline int csub(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
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
    fmpz_poly_sub(out, a, b)

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
    raise NotImplementedError

cdef inline int cmul(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
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
    fmpz_poly_mul(out, a, b)

cdef inline int cdivunit(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
    """
    Division.

    The inversion is perfomed modulo p^prec.  Note that no reduction
    is performed after the product.

    INPUT:

    - ``out`` -- an ``celement`` to store the quotient.
    - ``a`` -- an ``celement``, the first input.
    - ``b`` -- an ``celement``, the second input.
    - ``prec`` -- a long, the precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    raise NotImplementedError

cdef inline int csetone(celement out, PowComputer_class prime_pow) except -1:
    """
    Sets to 1.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 1.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    fmpz_poly_set_ui(out, 1)

cdef inline int csetzero(celement out, PowComputer_class prime_pow) except -1:
    """
    Sets to 0.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 0.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    fmpz_poly_set_ui(out, 0)

cdef inline bint cisone(celement a, PowComputer_class prime_pow_) except -1:
    """
    Returns whether this element is equal to 1.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a = 1`, and False otherwise.
    """
    return fmpz_poly_is_one(a)

cdef inline bint ciszero(celement a, PowComputer_class prime_pow) except -1:
    """
    Returns whether this element is equal to 0.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a = 0`, and False otherwise.
    """
    return fmpz_poly_is_zero(a)

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
    raise NotImplementedError

cdef inline int ccopy(celement out, celement a, PowComputer_class prime_pow) except -1:
    """
    Copying.

    INPUT:

    - ``out`` -- the ``celement`` to store the result.
    - ``a`` -- the element to copy.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    fmpz_poly_set(out, a)

cdef inline cpickle(celement a, PowComputer_class prime_pow):
    """
    Serialization into objects that Sage knows how to pickle.

    INPUT:

    - ``a`` the element to pickle.
    - ``prime_pow`` the PowComputer for the ring.

    OUTPUT:

    - a serializable object storing ``a``.
    """
    raise NotImplementedError

cdef inline int cunpickle(celement out, x, PowComputer_class prime_pow) except -1:
    """
    Reconstruction from the output of meth:`cpickle`.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result.
    - ``x`` -- the result of `meth`:cpickle.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    raise NotImplementedError

cdef inline long chash(celement a, long ordp, long prec, PowComputer_class prime_pow) except -1:
    """
    Hashing.

    INPUT:

    - ``a`` -- an ``celement`` storing the underlying element to hash.
    - ``ordp`` -- a long storing the valuation.
    - ``prec`` -- a long storing the precision.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    raise NotImplementedError

cdef clist(celement a, long prec, bint pos, PowComputer_class prime_pow_):
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
    cdef PowComputer_flint_unram prime_pow = <PowComputer_flint_unram>prime_pow_
    ret = []
    cdef Integer digit
    for i from 0 <= i <= fmpz_poly_degree(a):
        fmpz_poly_get_coeff_fmpz(prime_pow.ftmp, a, i)
        j = 0
        while True:
            if fmpz_is_zero(prime_pow.ftmp): break
            if len(ret) <= j: ret.append([])
            while len(ret[j]) <= i: ret[j].append([])
            digit = PY_NEW(Integer)
            fmpz_fdiv_qr(prime_pow.ftmp, prime_pow.ftmp2, prime_pow.ftmp, prime_pow.fprime)
            fmpz_get_mpz(digit.value, prime_pow.ftmp2)
            ret[j][i] =digit
            j += 1
    return ret

# The element is filled in for zero in the output of clist if necessary.
_list_zero = []

cdef int cteichmuller(celement out, celement value, long prec, PowComputer_class prime_pow) except -1:
    """
    Teichmuller lifting.

    INPUT:

    - ``out`` -- an ``celement`` which is set to a `q-1` root of unity
                 congruent to `value` mod `\pi`; or 0 if `a \equiv 0
                 \pmod{\pi}`.
    - ``value`` -- an ``celement``, the element mod `\pi` to lift.
    - ``prec`` -- a long, the precision to which to lift.
    - ``prime_pow`` -- the Powcomputer of the ring.
    """
    raise NotImplementedError

cdef int cconv(celement out, x, long prec, long valshift, PowComputer_class prime_pow) except -2:
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
    raise NotImplementedError

cdef inline long cconv_mpz_t(celement out, mpz_t x, long prec, bint absolute, PowComputer_class prime_pow) except -2:
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
    fmpz_poly_set_mpz(out, x)
    if absolute:
        creduce(out, out, prec, prime_pow)
    elif mpz_sgn(x) == 0:
        return maxordp
    else:
        return cremove(out, out, prec, prime_pow)

cdef inline int cconv_mpz_t_out(mpz_t out, celement x, long valshift, long prec, PowComputer_class prime_pow) except -1:
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
    raise NotImplementedError

cdef inline long cconv_mpq_t(celement out, mpq_t x, long prec, bint absolute, PowComputer_class prime_pow) except? -10000:
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
    raise NotImplementedError

cdef inline int cconv_mpq_t_out(mpq_t out, celement x, long valshift, long prec, PowComputer_class prime_pow) except -1:
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
    raise NotImplementedError
