"""

"""

from sage.rings.integer cimport Integer
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport *

include "sage/libs/linkages/padics/Polynomial_shared.pxi"

cdef inline int ccmp(celement a, celement b, long prec, bint reduce_a, bint reduce_b, PowComputer_ prime_pow) except -2:
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
    # need to account for reduce_a and reduce_b
    if a == b:
        return 0
    elif a < b:
        return -1
    else:
        return 1

cdef inline bint creduce(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
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
    out.__coeffs = (<celement?>(a % prime_pow.modulus)).__coeffs
    cdef long coeff_prec = prec / prime_pow.e + 1
    cdef long break_pt = prec % prime_pow.e
    if break_pt > len(out.__coeffs):
        break_pt = len(out.__coeffs)
    for i in range(break_pt):
        out.__coeffs[i] %= prime_pow.base_ring.uniformizer_pow(coeff_prec)
    coeff_prec -= 1
    for i in range(break_pt, len(out.__coeffs)):
        out.__coeffs[i] %= prime_pow.base_ring.uniformizer_pow(coeff_prec)
    out.__normalize()
    return out == 0

cdef inline bint creduce_small(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
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

cdef inline long cremove(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Extract the maximum power of the uniformizer dividing this element.

    INPUT:

    - ``out`` -- an ``celement`` to store the unit.
    - ``a`` -- the element whose valuation and unit are desired.
    - ``prec`` -- a long, used if `a = 0`.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - if `a = 0`, returns prec (the value of ``out`` is undefined).
      Otherwise, returns the number of times `\pi` divides `a`.
    """
    if a == 0:
        return prec
    cdef long v = cvaluation(a, prec, prime_pow)
    cshift(out, a, -v, prec, prime_pow, True)
    return v

cdef inline long cvaluation(celement a, long prec, PowComputer_ prime_pow) except -1:
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
    if a == 0:
        return prec
    return min(c.valuation()*prime_pow.e + i for i, c in enumerate(a.__coeffs))

cdef inline bint cisunit(celement a, PowComputer_ prime_pow) except -1:
    """
    Returns whether this element has valuation zero.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a` has valuation 0, and False otherwise.
    """
    return (cvaluation(a, 1, prime_pow) == 0)

cdef inline int cshift(celement out, celement a, long n, long prec, PowComputer_ prime_pow, bint reduce_afterward) except -1:
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
    cdef long q, r
    modulus = prime_pow.modulus
    cdef celement ans
    if n > 0:
        ans = a * prime_pow.poly_ring.gen()**n
        ans %= modulus
    elif n == 0:
        ans = a
    else:
        q = n / prime_pow.e
        r = n % prime_pow.e
        # Multiply by (p/x^e)^n
        if q:
            ans = a * prime_pow.pxe ** q # should do this modulo prime_pow.modulus, rather than reducing afterward
            ans %= modulus
        else:
            ans = a
        if r:
            K = modulus.base_ring().fraction_field()
            Qpmodulus = modulus.change_ring(K)
            R = Qpmodulus.parent()
            # split modulus in half:
            # modulus = p*c - x^r*d, where c and d are integral polynomials, and c has unit const term.
            # Then p/x^r = d/c
            c = ans[:r] / K.uniformizer()
            d = -R(ans.list()[r:])
            _, _, ci = Qpmodulus.xgcd(c)
            ans *= d * ci
            ans = ans / K.uniformizer()
            ans = ans.change_ring(modulus.base_ring())
            ans = ans % modulus
    if reduce_afterward:
        creduce(out, ans, prec, prime_pow)
    else:
        out.__coeffs = ans.__coeffs

cdef inline int cshift_notrunc(celement out, celement a, long n, long prec, PowComputer_ prime_pow) except -1:
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
    cshift(out, a, n, prec, prime_pow, True)

cdef inline int cinvert(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    """
    Inversion

    The result will be reduced modulo p^prec.

    INPUT:

    - ``out`` -- an ``celement`` to store the inverse.
    - ``a`` -- an ``celement``, the element to be inverted.
    - ``prec`` -- a long, the precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    K = prime_pow.base_ring.fraction_field()
    Qpmodulus = prime_pow.modulus.change_ring(K)
    cdef celement inv
    g, _, inv = Qpmodulus.xgcd(a)
    if g != 1:
        raise ArithmeticError("Not invertible")
    inv = inv.change_ring(prime_pow.base_ring) % prime_pow.modulus # The % may be unnecessary
    # Need to reduce modulo prec
    out.__coeffs = inv.__coeffs

cdef inline int cdivunit(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
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
    cinvert(out, b, prec, prime_pow)
    cmul(out, out, a, prec, prime_pow)

cdef inline int cpow(celement out, celement a, mpz_t n, long prec, PowComputer_ prime_pow) except -1:
    """
    Exponentiation.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result.
    - ``a`` -- the base.
    - ``n`` -- an ``mpz_t``, the exponent.
    - ``prec`` -- a long, the working absolute precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    # We do this the stupid way for now.
    cdef Integer zn = PY_NEW(Integer)
    mpz_set(zn.value, n)
    cdef celement ans = a**zn
    ans %= prime_pow.modulus
    out.__coeffs = ans.__coeffs

cdef clist(celement a, long prec, bint pos, PowComputer_ prime_pow):
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

    - A list of p-adic digits `[a_0, a_1, \ldots]` so that `a = a_0 + a_1*pi + \cdots` modulo `pi^{prec}`.
    """
    raise NotImplementedError

cdef int cteichmuller(celement out, celement value, long prec, PowComputer_ prime_pow) except -1:
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
    if value[0].valuation() > 0:
        out.__coeffs = []
    else:
        out.__coeffs = [value[0].parent().teichmuller(value[0])]

