# -*- coding: utf-8 -*-
r"""
This linkage file implements the padics API for ramified extensions using Sage
Polynomials.

It contains the bits that are specific for ramified extensions. Everything that
is independent of ramification is in Polynomial_shared.pxi.

.. NOTE::

    There are no doctests in this file since the functions here can not be
    called directly from Python. Testing of this function is necessarily
    indirect and mostly done through arithmetic black-box tests that are part
    of the test suites of the `p`-adic parents.

AUTHORS:

- David Roe, Julian Rüth (2017-06-11): initial version

"""
#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#                     2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.integer cimport Integer
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport *

include "sage/libs/linkages/padics/Polynomial_shared.pxi"

cdef inline bint creduce(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Reduce ``a`` modulo a power of the maximal ideal.

    INPUT:

    - ``out`` -- a ``celement`` to store the reduction

    - ``a`` -- the ``celement`` to be reduced

    - ``prec`` -- a ``long``, the precision to reduce modulo

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    OUTPUT:

    ``True`` if the reduction is zero, ``False`` otherwise

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
    r"""
    Reduce ``a`` modulo a power of the maximal ideal.

    Similar to ``creduce`` but this function assumes that at most one
    addition/subtraction has happened on reduced inputs.  For integral inputs
    this translates to the assumption that `-p^\mathrm{prec} < a < 2p^\mathrm{prec}`.

    INPUT:

    - ``out`` -- a ``celement`` to store the reduction

    - ``a`` -- the ``celement`` to be reduced

    - ``prec`` -- a ``long``, the precision to reduce modulo

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    OUTPUT:

    ``True`` if the reduction is zero, ``False`` otherwise

    """
    return creduce(out, a, prec, prime_pow)

cdef inline long cvaluation(celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Return the maximum power of the uniformizer dividing ``a``.

    This function differs from :meth:`cremove` in that the unit is discarded.

    INPUT:

    - ``a`` -- the element whose valuation is desired

    - ``prec`` -- a ``long``; if the valuation of ``a`` exceeds ``prec``, this
      function returns ``prec``. In particular, ``prec`` is returned if ``a``
      is zero.

    - ``prec`` -- a ``long``, the return value if ``a`` is zero

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    OUTPUT:

    The number of times the uniformizer divides ``a``, or ``prec`` if that is
    higher.

    """
    long ret = prec

    for i,c in enumerate(a.__coeffs):
        if i >= prec:
            break
        ret = min(ret, c.valuation()*prime_pow.e + i)

    return ret

cdef inline int cshift(celement out, celement a, long n, long prec, PowComputer_ prime_pow, bint reduce_afterward) except -1:
    r"""
    Multiply ``a`` with an ``n``-th power of the uniformizer.

    This function shifts the `\pi`-adic expansion of ``a`` by ``n``, i.e., it
    multiplies ``a`` by the `n`-th power of the uniformizer and drops any terms
    with negative powers of the uniformizer in the `\pi`-adic expansion.

    INPUT:

    - ``out`` -- a ``celement`` to store the result

    - ``a`` -- the ``celement`` to shift

    - ``n`` -- a ``long``, the amount to shift by

    - ``prec`` -- a ``long``, a precision modulo which to reduce if
      ``reduce_afterward`` is set

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    - ``reduce_afterward`` -- whether to :meth:`creduce` ``out`` before
      returning

    """
    cdef long q, r

    if n > 0:
        a *= prime_pow.uniformizer_pow(n)
    elif n < 0:
        q = -n / prime_pow.e # ≥ 0
        r = -n % prime_pow.e # ≥ 0
        # As 0 > n = -q*e - r, π^n = p^-q * (p/π^e)^q * π^-r
        if q:
            # Multiply with p^-q: this kills the digits in the π-adic expansion
            # that belong to terms below π^(q*e).
            a = a.map_coefficients(lambda c: c>>q)
            # Multiply with (p/π^e)^q.
            a *= prime_pow.pxe_pow(q)
        if r:
            # Multiply with (p/x^r)
            a *= prime_pow.px_pow(r)
            a %= prime_pow.modulus()
            # Divide by p
            a = a.map_coefficients(lambda c: c>>1)

    if reduce_afterward:
        creduce(out, a, prec, prime_pow)
    else:
        out._coeffs = a.__coeffs

cdef inline int cshift_notrunc(celement out, celement a, long n, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Multiply ``a`` with an ``n``-th power of the uniformizer.

    This method is identical to :meth:`cshift` but assumes that the valuation
    of ``a`` is at least ``-n``.

    INPUT:

    - ``out`` -- a ``celement`` to store the result

    - ``a`` -- the ``celement`` to shift

    - ``n`` -- a ``long``, the amount to shift by

    - ``prec`` -- a ``long``, a precision modulo which to reduce if
      ``reduce_afterward`` is set

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    """
    cshift(out, a, n, prec, prime_pow, True)

cdef inline int cinvert(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Compute the inverse of ``a``.

    INPUT:

    - ``out`` -- a ``celement`` to store the inverse

    - ``a`` -- a ``celement``, the element to be inverted

    - ``prec`` -- a ``long``, ``out`` is reduced to this precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    """
    out.__coeffs = prime_pow.invert(a, prec).__coeffs
    creduce(out, out, prec, prime_pow)

cdef inline int cdivunit(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Divide ``a`` by ``b``.

    This function computes ``a*(1/b)`` where the inverse of ``b`` is determined
    to precision ``prec``. No reduction is performed after the product.

    INPUT:

    - ``out`` -- a ``celement`` to store the quotient

    - ``a`` -- a ``celement``, the dividend

    - ``b`` -- a ``celement``, the divisor, an element of valuation zero

    - ``prec`` -- a ``long``, the precision to which the inverse of ``b`` is
      determined

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    """
    cinvert(out, b, prec, prime_pow)
    cmul(out, out, a, prec, prime_pow)

cdef inline int cpow(celement out, celement a, mpz_t n, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Raise ``a`` to the ``n``-th power.

    INPUT:

    - ``out`` -- a ``celement`` in which to store the result

    - ``a`` -- a ``celement``, the base

    - ``n`` -- an ``mpz_t``, the exponent

    - ``prec`` -- a ``long``, the working absolute precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    """
    cdef Integer zn = PY_NEW(Integer)
    mpz_set(zn.value, n)

    csetone(out, prime_pow)
    for digit in zn.binary():
        if digit == '1':
            cmul(out, out, a, prec, prime_pow)
        cmul(out, out, out, prec, prime_pow)
        # we should probably not creduce that frequently to increase the performance
        creduce(out, out, a, prec, prime_pow)

# The element is filled in for zero in the output of clist if necessary.
# This WON'T work if the absolute inertia degree is 1.
_list_zero = []

cdef clist(celement a, long prec, bint pos, PowComputer_ prime_pow):
    r"""
    Return a list of digits in the series expansion.

    This function is used in printing, and expresses ``a`` as a series
    in the standard uniformizer `\pi`.  Note that for extension
    elements, the coefficients in the series could be lists themselves.

    INPUT:

    - ``a`` -- a ``celement`` for which to determine the list expansion

    - ``prec`` -- a ``long``, the number of digits desired

    - ``pos`` -- if ``True``, then representatives in the range 0…p-1 are used;
      otherwise, in the range `-\lfloor p/2\rlfloor … \lceil p/2\rceil`

    - ``prime_pow`` -- a ``PowComputer`` for the ring

    OUTPUT:

    A list of p-adic digits `[a_0, a_1, \ldots]` so that `a = a_0 + a_1\pi +
    \cdots` modulo `\pi^\mathrm{prec}`

    """
    # This is not very efficient, but there's no clear better way.
    # We assume this is only called on two-step extensions (for more general
    # extensions, convert to the absolute field).
    R = a.base_ring()
    if R.degree() == 1:
        raise NotImplementedError("Absolute extensions using Sage polynomials not completely supported")
    if R.base_ring().degree() != 1:
        raise TypeError("clist only allowed on towers of height 2")
    cdef long j
    ans = []
    p = a.base_ring().prime()
    p2 = (p-1)//2
    ccopy(prime_pow.poly_clist, a, prime_pow)
    for j in range(prec):
        # the following is specific to the ramified over unramified case.
        const_term = prime_pow.poly_clist[0]
        if const_term._is_exact_zero():
            term = []
        else:
            flint_rep = const_term._flint_rep_abs()[0]
            term = [c % p for c in flint_rep.list()]
            while term and not term[-1]:
                del term[-1]
            if not pos:
                term = [c - p if c > p2 else c for c in term]
        ans.append(term)
        if j != prec-1:
            cshift(prime_pow.poly_clist, prime_pow.poly_clist, -1, prec-j, prime_pow, False)
        if not prime_pow.poly_clist:
            break
    return ans

cdef int cteichmuller(celement out, celement value, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Compute a Teichmüller representative congruent to ``value``.

    INPUT:

    - ``out`` -- a ``celement`` which is set to a `q-1`-th root of unity
      congruent to ``value`` modulo `\pi`; or 0 if `a \equiv 0 \pmod{\pi}`.

    - ``value`` -- n ``celement``, the element mod `\pi` to lift

    - ``prec`` -- a ``long``, the precision to which to lift

    - ``prime_pow`` -- the ``PowComputer`` of the ring

    """
    if value[0].valuation() > 0:
        out.__coeffs = []
    else:
        out.__coeffs = [value[0].parent().teichmuller(value[0])]
