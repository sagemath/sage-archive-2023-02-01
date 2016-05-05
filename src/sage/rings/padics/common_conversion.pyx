r"""
The functions in this file are used in creating new p-adic elements.

When creating a p-adic element, the user can specify that the absolute
precision be bounded and/or that the relative precision be bounded.
Moreover, different p-adic parents impose their own bounds on the
relative or absolute precision of their elements.  The precision
determines to what power of `p` the defining data will be reduced, but
the valuation of the resulting element needs to be determined before
the element is created.  Moreover, some defining data can impose their
own precision bounds on the result.

AUTHORS:

- David Roe (2012-03-01)
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.int cimport *
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.all cimport *
from sage.libs.gmp.rational_reconstruction cimport mpq_rational_reconstruction
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
import sage.rings.finite_rings.integer_mod
from sage.libs.pari.types cimport *
from sage.libs.pari.gen cimport gen as pari_gen
from sage.rings.infinity import infinity

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1
# The following Integer is used so that the functions here don't need to initialize an mpz_t.
cdef Integer temp = PY_NEW(Integer)

cdef long get_ordp(x, PowComputer_class prime_pow) except? -10000:
    """
    This function determines the valuation of the `p`-adic element
    that will result from the given data ``x``.

    Note that the valuation can differ depending on the ring: if the
    new `p`-adic element is being created in a ring with ramification
    then the valuation will be larger.

    Some preprocessing is done in the initialization methods of the element, so
    the type of x is restricted.

    Also note that for some kinds of elements conversion to and from
    Integers and Rationals is done in a custom morphism rather than
    through this function.

    INPUT:

    - ``x`` -- data defining a new p-adic element: a Python int, an
      Integer, Rational, an element of Zp or Qp with the same prime, a
      PARI p-adic element, a list, a tuple, or an IntegerMod.

    - a PowComputer associated to a `p`-adic ring, which determines
      the prime and the ramification degree.

    OUTPUT:

    - a long, giving the valuation of the resulting `p`-adic element.
      If the input is zero, returns ``maxordp``
    """
    cdef long k, n, p, curterm, shift, f, e = prime_pow.e
    cdef Integer value
    cdef GEN pari_tmp
    if isinstance(x, int):
        if x == 0:
            return maxordp
        try:
            n = PyInt_AsLong(x)
        except OverflowError:
            return get_ordp(Integer(x), prime_pow)
        else:
            if mpz_fits_slong_p(prime_pow.prime.value) == 0:
                # x is not divisible by p
                return 0
            p = mpz_get_si(prime_pow.prime.value)
            k = 0
            while n % p == 0:
                k += 1
                n = n / p
    elif isinstance(x, Integer):
        if mpz_sgn((<Integer>x).value) == 0:
            return maxordp
        k = mpz_remove(temp.value, (<Integer>x).value, prime_pow.prime.value)
    elif isinstance(x, Rational):
        if mpq_sgn((<Rational>x).value) == 0:
            return maxordp
        k = mpz_remove(temp.value, mpq_numref((<Rational>x).value), prime_pow.prime.value)
        if k == 0:
            k = -mpz_remove(temp.value, mpq_denref((<Rational>x).value), prime_pow.prime.value)
    elif isinstance(x, (list,tuple)):
        f = prime_pow.f
        if (e == 1 and len(x) > f) or (e != 1 and len(x) > e):
            # could reduce modulo the defining polynomial but that isn't currently supported
            raise ValueError("List too long")
        k = maxordp
        shift = 0
        for a in x:
            if isinstance(a, (list,tuple)):
                if e == 1 or f == 1:
                    raise ValueError("nested lists not allowed for unramified and eisenstein extensions")
                for b in a:
                    if isinstance(b, (list,tuple)):
                        raise ValueError("list nesting too deep")
                    curterm = get_ordp(b, prime_pow)
                    k = min(k, curterm + shift)
            else:
                curterm = get_ordp(a, prime_pow)
                k = min(k, curterm + shift)
            if e != 1: shift += 1
        # We don't want to multiply by e again.
        return k
    elif isinstance(x, pAdicGenericElement):
        k = (<pAdicGenericElement>x).valuation_c()
        if not (<pAdicGenericElement>x)._is_base_elt(prime_pow.prime):
            k //= x.parent().ramification_index()
    elif isinstance(x, pari_gen):
        pari_tmp = (<pari_gen>x).g
        if typ(pari_tmp) == t_PADIC:
            k = valp(pari_tmp)
        else: # t_INT and t_FRAC were converted before this function
            raise TypeError("unsupported coercion from pari: only p-adics, integers and rationals allowed")
    elif sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
        value = <Integer>x.lift()
        if mpz_sgn(value.value) == 0:
            return maxordp
        k = mpz_remove(temp.value, value.value, prime_pow.prime.value)
    else:
        from sage.structure.element import parent
        raise NotImplementedError("Can not determine p-adic valuation of an element of %s"%parent(x))
    # Should check for overflow
    return k * e

cdef long get_preccap(x, PowComputer_class prime_pow) except? -10000:
    """
    This function determines the maximum absolute precision possible
    for an element created from the given data ``x``.

    Note that the valuation can differ depending on the ring: if the
    new `p`-adic element is being created in a ring with ramification
    then the precision bound will be larger.

    Some preprocessing is done in the initialization methods of the element, so
    the type of x is restricted.

    Also note that for some kinds of elements conversion to and from
    Integers and Rationals is done in a custom morphism rather than
    through this function.

    INPUT:

    - ``x`` -- data defining a new p-adic element: an Integer,
      Rational, an element of Zp or Qp with the same prime, a PARI
      p-adic element, a list, a tuple, or an IntegerMod.
    - ``prime_pow`` -- the PowComputer for the ring into which ``x``
      is being converted.  This is used to determine the prime and the
      ramification degree.

    OUTPUT:

    - a long, giving the absolute precision modulo which the input is
      defined.  If the input is exact, returns ``maxordp``
    """
    cdef long k, shift, e = prime_pow.e
    cdef Integer prec
    cdef GEN pari_tmp
    if isinstance(x, int) or isinstance(x, Integer) or isinstance(x, Rational):
        return maxordp
    elif isinstance(x, (list,tuple)):
        k = maxordp
        shift = 0
        for a in x:
            if isinstance(a, (list,tuple)):
                for b in a:
                    curterm = get_preccap(b, prime_pow)
                    k = min(k, curterm + shift)
            else:
                curterm = get_preccap(a, prime_pow)
                k = min(k, curterm + shift)
            if e != 1: shift += 1
        # We don't want to multiply by e again.
        return k
    elif isinstance(x, pAdicGenericElement):
        if (<pAdicGenericElement>x)._is_exact_zero():
            return maxordp
        prec = <Integer>x.precision_absolute()
        k = mpz_get_si(prec.value)
        if not (<pAdicGenericElement>x)._is_base_elt(prime_pow.prime):
            # since x lives in a subfield, the ramification index of x's parent will divide e.
            return k * (e // x.parent().ramification_index())
    elif isinstance(x, pari_gen):
        pari_tmp = (<pari_gen>x).g
        # since get_ordp has been called typ(x.g) == t_PADIC
        k = valp(pari_tmp) + precp(pari_tmp)
    elif sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
        k = mpz_remove(temp.value, (<Integer>x.modulus()).value, prime_pow.prime.value)
        if mpz_cmp_ui(temp.value, 1) != 0:
            raise TypeError("cannot coerce from the given integer mod ring (not a power of the same prime)")
    else:
        from sage.structure.element import parent
        raise NotImplementedError("Can not determine p-adic precision of an element of %s"%parent(x))
    return k * e

cdef long comb_prec(iprec, long prec) except? -10000:
    """
    This function returns the minimum of iprec and prec.

    INPUT:

    - ``iprec`` -- either infinity, an Integer, a Python int or
      something that can be converted to an Integer.

    - ``prec`` -- a long.
    """
    if iprec is infinity: return prec
    cdef Integer intprec
    if isinstance(iprec, Integer):
        intprec = <Integer>iprec
        if mpz_cmp_si(intprec.value, prec) >= 0:
            return prec
        if mpz_fits_slong_p(intprec.value) == 0:
            raise OverflowError("precision overflow")
        return mpz_get_si(intprec.value)
    if isinstance(iprec, int):
        return min(PyInt_AS_LONG(iprec), prec)
    return comb_prec(Integer(iprec), prec)

cdef int _process_args_and_kwds(long *aprec, long *rprec, args, kwds, bint absolute, PowComputer_class prime_pow) except -1:
    """
    This function obtains values for absprec and relprec from a
    combination of positional and keyword arguments.

    When creating a p-adic element, the user can pass in two arguments: ``absprec`` and ``relprec``.
    In implementing morphisms to speed up conversion from Integers and Rationals,
    we need to determine ``absprec`` and ``relprec`` from the ``args`` and ``kwds`` arguments of
    ``_call_with_args``.  This function collects the code to do so.

    INPUT:

    - ``args`` -- a tuple of positional arguments (at most two)

    - ``kwds`` -- a dictionary of keyword arguments (only
      ``'relprec'`` and ``'absprec'`` are used)

    - ``absolute`` -- (boolean) True if the precision cap of the ring
      is a cap on absolute precision, False if a cap on relative
      precision.

    - ``prime_pow`` -- a
      :class:`sage.rings.padics.pow_computer.PowComputer_class`
      instance

    OUTPUT:

    - ``aprec`` -- (first argument) the maximum absolute precision of
      the resulting element.

    - ``rprec`` -- (second argument) the maximum relative precision of
      the resulting element.

    - error status
    """
    if "empty" in kwds:
        # For backward compatibility
        aprec[0] = 0
        rprec[0] = 0
        return 0
    if len(args) > 2:
        raise TypeError("too many positional arguments")
    if len(args) == 2:
        if "relprec" in kwds:
            raise TypeError("_call_with_args() got multiple values for keyword argument 'relprec'")
        relprec = args[1]
    else:
        relprec = kwds.get("relprec",infinity)
    if len(args) >= 1:
        if "absprec" in kwds:
            raise TypeError("_call_with_args() got multiple values for keyword argument 'absprec'")
        absprec = args[0]
    else:
        absprec = kwds.get("absprec",infinity)
    if absolute:
        aprec[0] = comb_prec(absprec, prime_pow.prec_cap)
        rprec[0] = comb_prec(relprec, maxordp)
    else:
        rprec[0] = comb_prec(relprec, prime_pow.prec_cap)
        aprec[0] = comb_prec(absprec, maxordp)

cdef inline long cconv_mpq_t_shared(mpz_t out, mpq_t x, long prec, bint absolute, PowComputer_class prime_pow) except? -10000:
    """
    A fast pathway for conversion of rationals that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the output.
    - ``x`` -- an ``mpq_t`` giving the integer to be converted.
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
    cdef long numval, denval
    cdef bint success
    if prec <= 0:
        raise ValueError
    if absolute:
        success = mpz_invert(out, mpq_denref(x), prime_pow.pow_mpz_t_tmp(prec))
        if not success:
            raise ValueError("p divides denominator")
        mpz_mul(out, out, mpq_numref(x))
        mpz_mod(out, out, prime_pow.pow_mpz_t_tmp(prec))
    elif mpq_sgn(x) == 0:
        mpz_set_ui(out, 0)
        return maxordp
    else:
        denval = mpz_remove(out, mpq_denref(x), prime_pow.prime.value)
        mpz_invert(out, out, prime_pow.pow_mpz_t_tmp(prec))
        if denval == 0:
            numval = mpz_remove(temp.value, mpq_numref(x), prime_pow.prime.value)
            mpz_mul(out, out, temp.value)
        else:
            numval = 0
            mpz_mul(out, out, mpq_numref(x))
        mpz_mod(out, out, prime_pow.pow_mpz_t_tmp(prec))
        return numval - denval

cdef inline int cconv_mpq_t_out_shared(mpq_t out, mpz_t x, long valshift, long prec, PowComputer_class prime_pow) except -1:
    """
    Converts the underlying `p`-adic element into a rational

    - ``out`` -- gives a rational approximating the input.  Currently uses rational reconstruction but
                 may change in the future to use a more naive method
    - ``x`` -- an ``mpz_t`` giving the underlying `p`-adic element
    - ``valshift`` -- a long giving the power of `p` to shift `x` by
    -` ``prec`` -- a long, the precision of ``x``, used in rational reconstruction
    - ``prime_pow`` -- a PowComputer for the ring
    """
    mpq_rational_reconstruction(out, x, prime_pow.pow_mpz_t_tmp(prec))

    # if valshift is nonzero then we starte with x as a p-adic unit,
    # so there will be no powers of p in the numerator or denominator
    # and the following operations yield reduced rationals.
    if valshift > 0:
        mpz_mul(mpq_numref(out), mpq_numref(out), prime_pow.pow_mpz_t_tmp(valshift))
    elif valshift < 0:
        mpz_mul(mpq_denref(out), mpq_denref(out), prime_pow.pow_mpz_t_tmp(-valshift))

cdef inline int cconv_shared(mpz_t out, x, long prec, long valshift, PowComputer_class prime_pow) except -2:
    """
    Conversion from other Sage types.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the output.

    - ``x`` -- a Sage element that can be converted to a `p`-adic element.

    - ``prec`` -- a long, giving the precision desired: absolute if
                  `valshift = 0`, relative if `valshift != 0`.

    - ``valshift`` -- the power of the uniformizer to divide by before
      storing the result in ``out``.

    - ``prime_pow`` -- a PowComputer for the ring.

    """
    if PyInt_Check(x):
        x = Integer(x)
    elif isinstance(x, pari_gen):
        x = x.sage()
    if isinstance(x, pAdicGenericElement) or sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
        x = x.lift()
    if isinstance(x, Integer):
        if valshift > 0:
            mpz_divexact(out, (<Integer>x).value, prime_pow.pow_mpz_t_tmp(valshift))
            mpz_mod(out, out, prime_pow.pow_mpz_t_tmp(prec))
        elif valshift < 0:
            raise RuntimeError("Integer should not have negative valuation")
        else:
            mpz_mod(out, (<Integer>x).value, prime_pow.pow_mpz_t_tmp(prec))
    elif isinstance(x, Rational):
        if valshift == 0:
            mpz_invert(out, mpq_denref((<Rational>x).value), prime_pow.pow_mpz_t_tmp(prec))
            mpz_mul(out, out, mpq_numref((<Rational>x).value))
        elif valshift < 0:
            mpz_divexact(out, mpq_denref((<Rational>x).value), prime_pow.pow_mpz_t_tmp(-valshift))
            mpz_invert(out, out, prime_pow.pow_mpz_t_tmp(prec))
            mpz_mul(out, out, mpq_numref((<Rational>x).value))
        else:
            mpz_invert(out, mpq_denref((<Rational>x).value), prime_pow.pow_mpz_t_tmp(prec))
            mpz_divexact(temp.value, mpq_numref((<Rational>x).value), prime_pow.pow_mpz_t_tmp(valshift))
            mpz_mul(out, out, temp.value)
        mpz_mod(out, out, prime_pow.pow_mpz_t_tmp(prec))
    elif isinstance(x, list):
        if valshift == 0:
            if len(x) == 0:
                cconv_shared(out, Integer(0), prec, valshift, prime_pow)
            elif len(x) == 1:
                cconv_shared(out, x[0], prec, valshift, prime_pow)
            else:
                raise NotImplementedError("conversion not implemented from non-prime residue field")
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError("No conversion defined for %s which is a %s in %s"%(x,type(x),x.parent() if hasattr(x,"parent") else "no parent"))

cdef inline long cconv_mpz_t_shared(mpz_t out, mpz_t x, long prec, bint absolute, PowComputer_class prime_pow) except -2:
    """
    A fast pathway for conversion of integers that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the output.
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
    cdef long val
    if absolute:
        mpz_mod(out, x, prime_pow.pow_mpz_t_tmp(prec))
    elif mpz_sgn(x) == 0:
        mpz_set_ui(out, 0)
        return maxordp
    else:
        val = mpz_remove(out, x, prime_pow.prime.value)
        mpz_mod(out, out, prime_pow.pow_mpz_t_tmp(prec))
        return val

cdef inline int cconv_mpz_t_out_shared(mpz_t out, mpz_t x, long valshift, long prec, PowComputer_class prime_pow) except -1:
    """
    Converts the underlying `p`-adic element into an integer if
    possible.

    - ``out`` -- stores the resulting integer as an integer between 0
      and `p^{prec + valshift}`.
    - ``x`` -- an ``mpz_t`` giving the underlying `p`-adic element.
    - ``valshift`` -- a long giving the power of `p` to shift `x` by.
    -` ``prec`` -- a long, the precision of ``x``: currently not used.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    if valshift == 0:
        mpz_set(out, x)
    elif valshift < 0:
        raise ValueError("negative valuation")
    else:
        mpz_mul(out, x, prime_pow.pow_mpz_t_tmp(valshift))
