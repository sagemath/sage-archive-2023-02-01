"""
This file implements elements of eisenstein and unramified extensions of Zp and Qp with capped relative precision.
For the parent class see padic_extension_leaves.pyx.

The underlying implementation is through NTL's ZZ_pX class.  Each element contains the following data:
  ordp (long)    -- A power of the uniformizer to scale the unit by.  For unramified extensions this uniformizer is p,
                     for eisenstein extensions it is not.  A value equal to the maximum value of a long indicates that
                     the element is an exact zero.
  relprec (long) -- A signed integer giving the precision to which this element is defined.  For nonzero relprecs, the absolute value
                     gives the power of the uniformizer modulo which the unit is defined.  A positive value indicates
                     that the element is normalized (ie unit is actually a unit: in the case of eisenstein extensions
                     the constant term is not divisible by p, in the case of unramified extensions that there is at
                     least one coefficient that is not divisible by p).  A negative value indicates that the element
                     may or may not be normalized.  A zero value indicates that the element is zero to some precision.
                     If so, ordp gives the absolute precision of the element.  If ordp is the maximum value for a long,
                     then the element is an exact zero.
  unit (ZZ_pX_c) -- An ntl ZZ_pX storing the unit part.  The varible x is the uniformizer in the case of eisenstein extensions.
                     If the element is not normalized, the "unit"  may or may not actually
                     be a unit.  This ZZ_pX is created with global ntl modulus determined by the absolute value of
                     relprec.  If relprec is 0, unit IS NOT INITIALIZED, or destructed if normalized and found to be
                     zero.  Otherwise, let r be relprec and e be the ramification index over Qp or Zp.  Then the modulus
                     of unit is given by p^ceil(r/e).  Note that all kinds of problems arise if you try to mix moduli.
                     ZZ_pX_conv_modulus gives a semi-safe way to convert between different moduli without having
                     to pass through ZZX (see sage/libs/ntl/decl.pxi and c_lib/src/ntl_wrap.cpp)
  prime_pow (some subclass of PowComputer_ZZ_pX) -- a class, identical among all elements with the same parent, holding
                     common data.
    prime_pow.deg -- The degree of the extension
    prime_pow.e   -- The ramification index
    prime_pow.f   -- The inertia degree
    prime_pow.prec_cap -- the unramified precision cap.  For eisenstein extensions this is the smallest power of p that is zero.
    prime_pow.ram_prec_cap -- the ramified precision cap.  For eisenstein extensions this will be the smallest power of x that is
                     indistinugishable from zero.
    prime_pow.pow_ZZ_tmp, prime_pow.pow_mpz_t_tmp, prime_pow.pow_Integer -- functions for accessing powers of p.
                     The first two return pointers.  See sage/rings/padics/pow_computer_ext for examples and important warnings.
    prime_pow.get_context, prime_pow.get_context_capdiv, prime_pow.get_top_context -- obtain an ntl_ZZ_pContext_class corresponding to p^n.
                     The capdiv version divides by prime_pow.e as appropriate.  top_context corresponds to prec_cap.
    prime_pow.restore_context, prime_pow.restore_context_capdiv, prime_pow.restore_top_context -- restores the given context.
    prime_pow.get_modulus, get_modulus_capdiv, get_top_modulus -- Returns a ZZ_pX_Modulus_c* pointing to a polynomial modulus defined modulo
                     p^n (appropriately divided by prime_pow.e in the capdiv case).

EXAMPLES:
An eisenstein extension:
    sage: R = Zp(5,5)
    sage: S.<x> = R[]
    sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
    sage: W.<w> = R.ext(f); W
    Eisenstein Extension of 5-adic Ring with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (3*5^2 + O(5^7))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^8))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
    sage: z = (1+w)^5; z
    1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
    sage: y = z >> 1; y
    w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
    sage: y.valuation()
    4
    sage: y.precision_relative()
    20
    sage: y.precision_absolute()
    24
    sage: z - (y << 1)
    1 + O(w^25)
    sage: (1/w)^12+w
    w^-12 + w + O(w^13)
    sage: (1/w).parent()
    Eisenstein Extension of 5-adic Field with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (3*5^2 + O(5^7))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^8))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))

An unramified extension:
    sage: g = x^3 + 3*x + 3
    sage: A.<a> = R.ext(g)
    sage: z = (1+a)^5; z
    (2*a^2 + 4*a) + (3*a^2 + 3*a + 1)*5 + (4*a^2 + 3*a + 4)*5^2 + (4*a^2 + 4*a + 4)*5^3 + (4*a^2 + 4*a + 4)*5^4 + O(5^5)
    sage: z - 1 - 5*a - 10*a^2 - 10*a^3 - 5*a^4 - a^5
    O(5^5)
    sage: y = z >> 1; y
    (3*a^2 + 3*a + 1) + (4*a^2 + 3*a + 4)*5 + (4*a^2 + 4*a + 4)*5^2 + (4*a^2 + 4*a + 4)*5^3 + O(5^4)
    sage: 1/a
    (3*a^2 + 4) + (a^2 + 4)*5 + (3*a^2 + 4)*5^2 + (a^2 + 4)*5^3 + (3*a^2 + 4)*5^4 + O(5^5)

Different printing modes:
    sage: R = Zp(5, print_mode='digits'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; repr(z)
    '...4110403113210310442221311242000111011201102002023303214332011214403232013144001400444441030421100001'
    sage: R = Zp(5, print_mode='bars'); S.<x> = R[]; g = x^3 + 3*x + 3; A.<a> = R.ext(g)
    sage: z = (1+a)^5; repr(z)
    '...[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 3, 4]|[1, 3, 3]|[0, 4, 2]'
    sage: R = Zp(5, print_mode='terse'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; z
    6 + 95367431640505*w + 25*w^2 + 95367431640560*w^3 + 5*w^4 + O(w^100)
    sage: R = Zp(5, print_mode='val-unit'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: y = (1+w)^5 - 1; y
    w^5 * (2090041 + 19073486126901*w + 1258902*w^2 + 674*w^3 + 16785*w^4) + O(w^100)

You can get at the underlying ntl unit:
    sage: z._ntl_rep()
    [6 95367431640505 25 95367431640560 5]
    sage: y._ntl_rep()
    [2090041 19073486126901 1258902 674 16785]
    sage: y._ntl_rep_abs()
    ([5 95367431640505 25 95367431640560 5], 0)

NOTES:
    If you get an error 'internal error: can't grow this _ntl_gbigint,' it indicates that moduli are being mixed
    inappropriately somewhere.
    For example, when calling a function with a ZZ_pX_c as an argument, it copies.  If the modulus is not set
    to the modulus of the ZZ_pX_c, you can get errors.

AUTHORS:
    -- David Roe  (2008-01-01) initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/stdsage.pxi"
include "../../ext/interrupt.pxi"

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.libs.pari.gen import gen as pari_gen
from sage.rings.all import is_IntegerMod
from sage.rings.padics.padic_ext_element cimport pAdicExtElement
from sage.rings.padics.precision_error import PrecisionError

from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX_small_Eis
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX_big_Eis

from sage.rings.real_double cimport RealDoubleElement

cdef object infinity
from sage.rings.infinity import infinity

cdef long maxordp = (1 << (sizeof(long) * 8 - 2)) -1

cdef class pAdicZZpXCRElement(pAdicZZpXElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, empty = False):
        """
        Creates an element of a capped relative precision, unramified or eisenstein extension of Zp or Qp.

        INPUT:
        parent -- either an EisensteinRingCappedRelative or UnramifiedRingCappedRelative
        x -- an integer, rational, p-adic element, polynomial, list, integer_mod, pari int/frac/poly_t/pol_mod, an ntl_ZZ_pX, an ntl_ZZ, an ntl_ZZ_p or an ntl_ZZX
        absprec -- an upper bound on the absolute precision of the element created
        relprec -- an upper bound on the relative precision of the element created
        empty -- whether to return after initializing to zero (without setting the valuation).

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1+w)^5; z # indirect doctest
        1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        """
        pAdicZZpXElement.__init__(self, parent)
        self.relprec = 0
        if empty:
            return
        cdef long aprec, rprec, ctx_prec
        if relprec is not infinity and not PY_TYPE_CHECK(relprec, Integer):
            relprec = Integer(relprec)
        if (relprec is infinity) or (relprec > parent.precision_cap()):
            rprec = self.prime_pow.ram_prec_cap
        else:
            rprec = mpz_get_si((<Integer>relprec).value)
        if rprec < 0:
            rprec = 0
        if absprec is not infinity:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                absprec = infinity
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
        cdef mpz_t tmp
        cdef ZZ_c tmp_z
        cdef Py_ssize_t i
        cdef Integer tmp_Int
        if PY_TYPE_CHECK(x, pAdicGenericElement):
            if self.prime_pow.in_field == 0 and x.valuation() < 0:
                raise ValueError, "element has negative valuation"
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes."
        if PY_TYPE_CHECK(x, pAdicBaseGenericElement):
            mpz_init(tmp)
            (<pAdicBaseGenericElement>x)._set_to_mpz(tmp)
            if absprec is infinity:
                self._set_from_mpz_rel(tmp, rprec)
            else:
                self._set_from_mpz_both(tmp, aprec, rprec)
            mpz_clear(tmp)
            return
        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                x = x.lift()
            if x.type() == 't_INT':
                x = Integer(x)
            elif x.type() == 't_FRAC':
                x = Rational(x)
            elif x.type() == 't_POLMOD' or x.type == 't_POL':
                # This code doesn't check to see if the primes are the same.
                L = []
                x = x.lift().lift()
                for i from 0 <= i <= x.poldegree():
                    L.append(Integer(x.polcoeff(i)))
                x = L
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers, rationals, polynomials and pol_mods allowed"
        elif is_IntegerMod(x):
            mpz_init(tmp)
            ctx_prec = mpz_remove(tmp, (<Integer>x.modulus()).value, self.prime_pow.prime.value)
            if mpz_cmp_ui(tmp, 1):
                mpz_clear(tmp)
                x = x.lift()
                if absprec is infinity or ctx_prec < aprec:
                    aprec = ctx_prec
                    absprec = Integer(1) # absprec just has to be non-infinite: everything else uses aprec
            else:
                mpz_clear(tmp)
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"
        elif PY_TYPE_CHECK(x, ntl_ZZ_p):
            ctx_prec = ZZ_remove(tmp_z, (<ntl_ZZ>x.modulus()).x, self.prime_pow.pow_ZZ_tmp(1)[0])
            if ZZ_IsOne(tmp_z):
                x = x.lift()
                tmp_Int = PY_NEW(Integer)
                ZZ_to_mpz(&tmp_Int.value, &(<ntl_ZZ>x).x)
                x = tmp_Int
                if absprec is infinity or ctx_prec < aprec:
                    aprec = ctx_prec
                    absprec = Integer(1) # absprec just has to be non-infinite: everything else uses aprec
            else:
                raise TypeError, "cannot coerce the given ntl_ZZ_p (modulus not a power of the same prime)"
        elif PY_TYPE_CHECK(x, ntl_ZZ):
            tmp_Int = PY_NEW(Integer)
            ZZ_to_mpz(&tmp_Int.value, &(<ntl_ZZ>x).x)
            x = tmp_Int
        elif isinstance(x, (int, long)):
            x = Integer(x)
        cdef pAdicZZpXCRElement _x
        if PY_TYPE_CHECK(x, Integer):
            if absprec is infinity:
                self._set_from_mpz_rel((<Integer>x).value, rprec)
            else:
                self._set_from_mpz_both((<Integer>x).value, aprec, rprec)
        elif PY_TYPE_CHECK(x, Rational):
            if absprec is infinity:
                self._set_from_mpq_rel((<Rational>x).value, rprec)
            else:
                self._set_from_mpq_both((<Rational>x).value, aprec, rprec)
        elif PY_TYPE_CHECK(x, ntl_ZZ_pX):
            if absprec is infinity:
                self._set_from_ZZ_pX_rel(&(<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c, rprec)
            else:
                self._set_from_ZZ_pX_both(&(<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c, aprec, rprec)
        elif PY_TYPE_CHECK(x, ntl_ZZX):
            if absprec is infinity:
                self._set_from_ZZX_rel((<ntl_ZZX>x).x, rprec)
            else:
                self._set_from_ZZX_both((<ntl_ZZX>x).x, aprec, rprec)
        elif PY_TYPE_CHECK(x, pAdicExtElement):
            if x.parent() is parent:
                _x = <pAdicZZpXCRElement>x
                self._set(&_x.unit, _x.ordp, _x.relprec)
            elif x.parent().fraction_field() is parent:
                if PY_TYPE_CHECK(x, pAdicZZpXCRElement):
                    _x = <pAdicZZpXCRElement>x
                    if _x.relprec < 0:
                        _x._normalize()
                    if _x._is_exact_zero():
                        self._set_exact_zero()
                    elif _x._is_inexact_zero():
                        self._set_inexact_zero(_x.ordp)
                    else:
                        if _x.relprec < rprec:
                            rprec = _x.relprec
                        self._set(&_x.unit, _x.ordp, rprec)
                else:
                    # x is a pAdicZZpXCAElement
                    xordp = x.valuation()
                    xprec = x.precision_absolute()
                    if xordp == xprec:
                        self._set_inexact_zero(mpz_get_si((<Integer>xordp).value))
                    else:
                        poly = x._ntl_rep_abs()[0]
                        if absprec is infinity:
                            self._set_from_ZZ_pX_rel(&(<ntl_ZZ_pX>poly).x,(<ntl_ZZ_pX>poly).c, rprec)
                        else:
                            self._set_from_ZZ_pX_both(&(<ntl_ZZ_pX>poly).x,(<ntl_ZZ_pX>poly).c, aprec, rprec)
            elif x.parent() is parent.fraction_field():
                _x = <pAdicZZpXCRElement>x
                if _x.relprec < 0:
                    _x._normalize()
                if _x._is_exact_zero():
                    self._set_exact_zero()
                elif _x._is_inexact_zero():
                    self._set_inexact_zero(_x.ordp)
                else:
                    if _x.relprec < rprec:
                        rprec = _x.relprec
                    self._set(&_x.unit, _x.ordp, rprec)
            else:
                raise NotImplementedError, "Conversion from different p-adic extensions not yet supported"
        else:
            try:
                x = list(x)
            except TypeError:
                try:
                    x = x.list()
                except AttributeError:
                    raise TypeError, "cannot convert x to a p-adic element"
            if absprec is infinity:
                self._set_from_list_rel(x, rprec)
            else:
                self._set_from_list_both(x, aprec, rprec)

    cdef void _set_inexact_zero(self, long absprec):
        """
        Sets self to be zero with valuation absprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(0,6); z # indirect doctest
        O(w^6)
        sage: z.valuation()
        6
        sage: z.precision_absolute()
        6
        sage: z.precision_relative()
        0
        """
        self.ordp = absprec
        self.relprec = 0

    cdef void _set_exact_zero(self):
        """
        Sets self to be an exact zero.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = R(0); z # indirect doctest
        0
        sage: z.valuation()
        +Infinity
        sage: z.precision_absolute()
        +Infinity
        sage: z.precision_relative()
        0
        """
        self.ordp = maxordp
        self.relprec = 0

    cpdef bint _is_exact_zero(self):
        """
        Tests if self is an exact zero.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(0)
        sage: z._is_exact_zero()
        True
        sage: z = W(0,6)
        sage: z._is_exact_zero()
        False
        """
        if self.ordp == maxordp:
            return 1
        else:
            return 0

    cpdef bint _is_inexact_zero(self):
        """
        Tests if self is an inexact zero.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(0)
        sage: z._is_inexact_zero()
        False
        sage: z = W(0,6)
        sage: z._is_inexact_zero()
        True
        """
        self._normalize()
        if self.relprec == 0:
            return not self._is_exact_zero()
        else:
            return False

    cdef void _set(self, ZZ_pX_c* unit, long ordp, long relprec):
        """
        Sets unit, ordp and relprec directly.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: F = W.fraction_field()
        sage: z = F(1+w); z # indirect doctest
        1 + w + O(w^25)
        """
        self.ordp = ordp
        self._set_prec_rel(relprec)
        self.prime_pow.restore_context_capdiv(self.relprec)
        ZZ_pX_conv_modulus(self.unit, unit[0], self.prime_pow.get_context_capdiv(relprec).x)

    cdef int _set_from_mpz_rel(self, mpz_t x, long relprec) except -1:
        """
        Sets self from an mpz_t with relative precision bounded by relprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: W(70, relprec = 8) # indirect doctest
        4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
        """
        if mpz_sgn(x) == 0:
            self._set_exact_zero()
            return 0
        cdef mpz_t tmp_m
        cdef ZZ_c tmp_z
        cdef long shift
        mpz_init(tmp_m)
        _sig_on
        shift = mpz_remove(tmp_m, x, self.prime_pow.prime.value)
        _sig_off
        self._set_prec_rel(relprec)
        mpz_to_ZZ(&tmp_z, &tmp_m)
        mpz_clear(tmp_m)
        self.prime_pow.restore_context_capdiv(relprec)
        ZZ_pX_SetCoeff(self.unit, 0, ZZ_to_ZZ_p(tmp_z))
        self.ordp = 0
        self._pshift_self(shift)

    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, long relprec) except -1:
        """
        Sets self from an mpz_t with relative precision bounded by relprec
        and absolute precision bounded by absprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: W(70, 8) # indirect doctest
        4*w^5 + 3*w^7 + O(w^8)
        """
        if mpz_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return 0
        cdef mpz_t tmp_m
        cdef ZZ_c tmp_z
        cdef long shift
        mpz_init(tmp_m)
        _sig_on
        shift = mpz_remove(tmp_m, x, self.prime_pow.prime.value)
        _sig_off
        self.ordp = shift * self.prime_pow.e
        if self._set_prec_both(absprec, relprec) == 1:
            # This indicates that self._set_inexact_zero was called
            mpz_clear(tmp_m)
            return 0
        mpz_to_ZZ(&tmp_z, &tmp_m)
        mpz_clear(tmp_m)
        ZZ_pX_SetCoeff(self.unit, 0, ZZ_to_ZZ_p(tmp_z))
        self.ordp = 0
        self._pshift_self(shift)

    cdef int _set_from_mpq_rel(self, mpq_t x, long relprec) except -1:
        """
        Sets self from an mpq_t with relative precision bounded by relprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(70/3, relprec = 9); z # indirect doctest
        3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + O(w^14)
        sage: z * 3
        4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + O(w^14)
        sage: W(70)
        4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23 + 3*w^25 + w^27 + O(w^30)
        sage: F = W.fraction_field()
        sage: y = F(3/700); y
        w^-10 + w^-8 + 4*w^-6 + w^-3 + 4*w^-2 + 3*w^-1 + 3 + 4*w + w^3 + 4*w^4 + w^5 + 4*w^6 + 2*w^7 + 3*w^8 + 4*w^9 + 3*w^10 + 4*w^11 + w^12 + O(w^15)
        sage: y * 700
        3 + O(w^25)
        """
        if mpq_sgn(x) == 0:
            self._set_exact_zero()
            return 0
        cdef mpz_t num_unit, den_unit
        self._set_from_mpq_part1(num_unit, den_unit, x)
        self._set_prec_rel(relprec)
        self._set_from_mpq_part2(num_unit, den_unit)

    cdef int _set_from_mpq_both(self, mpq_t x, long absprec, long relprec) except -1:
        """
        Sets self from an mpq_t with relative precision bounded by relprec and absolute precision bounded by absprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(70/3, 14); z # indirect doctest
        3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + O(w^14)
        sage: z * 3
        4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + O(w^14)
        sage: W(70)
        4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23 + 3*w^25 + w^27 + O(w^30)
        sage: F = W.fraction_field()
        sage: y = F(3/700,-2); y
        w^-10 + w^-8 + 4*w^-6 + w^-3 + O(w^-2)
        sage: y * 700
        3 + O(w^8)
        """
        if mpq_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return 0
        cdef mpz_t num_unit, den_unit
        self._set_from_mpq_part1(num_unit, den_unit, x)
        if self._set_prec_both(absprec, relprec) == 1:
            # indicates an inexact zero
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            return 0
        self._set_from_mpq_part2(num_unit, den_unit)

    cdef int _set_from_mpq_part1(self, mpz_t num_unit, mpz_t den_unit, mpq_t x) except -1:
        """
        Sets num_unit to be the unit of the numerator, den_unit to be the unit of the denominator and sets self.ordp correctly.

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(7000/3, 23); z # indirect doctest
        2*w^15 + 2*w^17 + 3*w^19 + w^22 + O(w^23)
        """
        cdef long num_ordp, den_ordp
        _sig_on
        mpz_init(num_unit)
        mpz_init(den_unit)
        num_ordp = mpz_remove(num_unit, mpq_numref(x), self.prime_pow.prime.value)
        den_ordp = mpz_remove(den_unit, mpq_denref(x), self.prime_pow.prime.value)
        _sig_off
        self.ordp = (num_ordp - den_ordp) * self.prime_pow.e
        if self.ordp < 0 and self.prime_pow.in_field == 0:
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            raise ValueError, "p divides the denominator"

    cdef int _set_from_mpq_part2(self, mpz_t num_unit, mpz_t den_unit) except -1:
        """
        Given that self.ordp and self.relprec have been set, takes num_unit and den_unit and sets self.unit.

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: W(QQ(0), 23) # indirect doctest
        O(w^23)
        sage: W(QQ(0))
        0
        """
        cdef ZZ_c num_zz, den_zz
        cdef ZZ_p_c tmp_zp
        cdef long val = self.ordp / self.prime_pow.e
        cdef mpz_t tmp_m
        mpz_init(tmp_m)
        mpz_set(tmp_m, num_unit)
        mpz_to_ZZ(&num_zz, &tmp_m)
        mpz_set(tmp_m, den_unit)
        mpz_to_ZZ(&den_zz, &tmp_m)
        self.prime_pow.restore_context_capdiv(self.relprec)
        ZZ_p_div(tmp_zp, ZZ_to_ZZ_p(num_zz), ZZ_to_ZZ_p(den_zz))
        ZZ_pX_SetCoeff(self.unit, 0, tmp_zp)
        self.ordp = 0
        self._pshift_self(val)

    cdef int _set_from_ZZX_rel(self, ZZX_c poly, long relprec) except -1:
        """
        Sets self from a ZZX with relative precision bounded by relprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(ntl.ZZX([4,1,16]), relprec = 14); z # indirect doctest
        4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + 4*w^13 + O(w^14)
        sage: z._ntl_rep()
        [4 1 16]
        sage: z = W(ntl.ZZX([5^40,5^42,3*5^41]), relprec = 14); z
        w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + O(w^214)
        sage: W(5)^40 + w*W(5)^42 + w^2 * W(3) * W(5)^41
        w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + 2*w^215 + w^217 + 2*w^218 + w^220 + w^221 + w^222 + 3*w^224 + O(w^225)
        """
        if ZZX_IsZero(poly):
            self._set_exact_zero()
            return 0
        if ZZX_deg(poly) >= self.prime_pow.deg:
            raise NotImplementedError
        # the -1 in the next line signals that there is no absprec specified
        if self._set_from_ZZX_part1(poly, -1, relprec) == -2:
            # indicates _set_inexact_zero was called
            return 0
        # context was restored in _set_from_ZZX_part1
        ZZX_to_ZZ_pX(self.unit, poly)
        self._internal_lshift(-self.ordp)

    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, long relprec) except -1:
        """
        Sets self from a ZZX with relative precision bounded by relprec and absolute precision bounded by absprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(ntl.ZZX([4,1,16]), 12); z # indirect doctest
        4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + O(w^12)
        sage: z._ntl_rep()
        [4 1 16]
        sage: z = W(ntl.ZZX([5^40,5^42,3*5^41]), 212); z
        w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + O(w^212)
        """
        if ZZX_IsZero(poly) or absprec <= 0:
            self._set_inexact_zero(absprec)
            return 0
        if ZZX_deg(poly) >= self.prime_pow.deg:
            raise NotImplementedError
        if self._set_from_ZZX_part1(poly, absprec, relprec) == -2:
            # indicates _set_inexact_zero was called
            return 0
        # context was restored in _set_from_ZZX_part1
        ZZX_to_ZZ_pX(self.unit, poly)
        self._internal_lshift(-self.ordp)

    cdef int _set_from_ZZX_part1(self, ZZX_c poly, long absprec, long relprec) except -1:
        """
        Sets self.ordp from poly and restores the context.  poly must have degree less than self.prime_pow.deg

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(ntl.ZZX([4,1,16]), 12); z # indirect doctest
        4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + O(w^12)
        """
        cdef long i = 0
        cdef long deg = ZZX_deg(poly)
        cdef long mini = -1
        cdef long minval
        cdef long curval
        cdef ZZ_c tmp_z
        while mini == -1:
            if not ZZ_IsZero(ZZX_coeff(poly,i)):
                minval = ZZ_remove(tmp_z, ZZX_coeff(poly, i), self.prime_pow.pow_ZZ_tmp(1)[0])
                mini = i
            i += 1
        while i <= deg:
            if not ZZ_IsZero(ZZX_coeff(poly,i)):
                curval = ZZ_remove(tmp_z, ZZX_coeff(poly, i), self.prime_pow.pow_ZZ_tmp(1)[0])
                if curval < minval:
                    minval = curval
                    mini = i
            i += 1
        if self.prime_pow.e == 1:
            self.ordp = minval
        else:
            self.ordp = minval * self.prime_pow.e + mini
        if absprec == -1: # indicates that _set_from_ZZX_rel is calling
            self._set_prec_rel(relprec)
        elif self._set_prec_both(absprec, relprec):
            # indicates self._set_inexact_zero was called
            return -2
        # We restore the context appropriately so that part2 works.
        self.prime_pow.restore_context_capdiv(self.relprec + self.ordp)

    cdef int _set_from_ZZ_pX_rel(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long relprec) except -1:
        """
        Sets self from a ZZ_pX with relative precision bounded by relprec.
        If ctx is None and poly is 0 this function will raise an error (a ZZ_pX cannot represent something with infinite absolute precision).

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(ntl.ZZ_pX([4,1,16],5^2)); z # indirect doctest
        4 + w + w^2 + 3*w^7 + w^9 + O(w^10)
        sage: z._ntl_rep()
        [4 1 16]
        sage: z = W(ntl.ZZ_pX([5^40,5^42,3*5^41], 5^44)); z
        w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + 2*w^215 + w^217 + 2*w^218 + O(w^220)
        """
        cdef long ctx_prec = -1
        if ctx is not None:
            ctx_prec = self._check_ZZ_pContext(ctx) * self.prime_pow.e
        if ZZ_pX_IsZero(poly[0]):
            if ctx_prec == -1:
                raise ValueError, "must specify either a context or an absolute precision bound"
            else:
                self._set_inexact_zero(ctx_prec)
            return 0
        self._set_from_ZZ_pX_part1(poly)
        if ctx_prec == -1:
            self._set_prec_rel(relprec)
        else:
            self._set_prec_both(ctx_prec, relprec)
        self._set_from_ZZ_pX_part2(poly)

    cdef int _set_from_ZZ_pX_both(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec, long relprec) except -1:
        """
        Sets self from a ZZ_pX with relative precision bounded by relprec and absolute precision bounded by absprec.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(ntl.ZZ_pX([4,1,16],5^2), absprec = 8, relprec = 12); z # indirect doctest
        4 + w + w^2 + 3*w^7 + O(w^8)
        sage: z._ntl_rep()
        [4 1 16]
        sage: z = W(ntl.ZZ_pX([5^40,5^42,3*5^41], 5^50), 220); z
        w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + 2*w^215 + w^217 + 2*w^218 + O(w^220)
        """
        cdef long ctx_prec
        if ctx is not None:
            ctx_prec = self._check_ZZ_pContext(ctx)
            if ctx_prec * self.prime_pow.e < absprec:
                absprec = ctx_prec * self.prime_pow.e
        if ZZ_pX_IsZero(poly[0]):
            self._set_inexact_zero(absprec)
            return 0
        self._set_from_ZZ_pX_part1(poly)
        self._set_prec_both(absprec, relprec)
        self._set_from_ZZ_pX_part2(poly)

    cdef int _set_from_ZZ_pX_part1(self, ZZ_pX_c* poly) except -1:
        """
        Sets self.ordp based on poly.  poly must not be 0.

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(ntl.ZZ_pX([4,1,16],5^2), absprec = 8, relprec = 12); z # indirect doctest
        4 + w + w^2 + 3*w^7 + O(w^8)
        """
        cdef long val, index
        ZZ_pX_min_val_coeff(val, index, poly[0], self.prime_pow.pow_ZZ_tmp(1)[0])
        if self.prime_pow.e == 1:
            self.ordp = val
        else:
            self.ordp = val * self.prime_pow.e + index

    cdef int _set_from_ZZ_pX_part2(self, ZZ_pX_c* poly) except -1:
        """
        Assuming that self.ordp and self.relprec have been set, sets self.unit to poly and then normalizes.

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(ntl.ZZ_pX([4,1,16],5^2), absprec = 8, relprec = 12); z # indirect doctest
        4 + w + w^2 + 3*w^7 + O(w^8)
        """
        self.prime_pow.restore_context_capdiv(self.ordp + self.relprec)
        ZZ_pX_conv_modulus(self.unit, poly[0], self.prime_pow.get_context_capdiv(self.ordp + self.relprec).x)
        self._internal_lshift(-self.ordp)

    cdef bint _set_prec_rel(self, long relprec):
        """
        Safely sets the relative precision of self to be the absolute value of relprec.

        Returns True iff self.relprec was reset.

        Note that this will wipe out anything in self.unit.  Be careful
        resetting self.unit directly: if you set it to a different modulus, NTL may
        have problems.  The safest way to reset self.unit to a different modulus is:
        self.prime_pow.restore_context_capdiv(self.relprec)
        cdef ZZ_pX_c tmp = self.unit
        self._set_prec_rel(new_rel_prec)
        ZZ_pX_conv_modulus(self.unit, tmp, self.prime_pow.get_context_capdiv(self.relprec).x)

        You may be able to just set self.relprec and ZZ_pX_conv_modulus if you're decreasing
        precision.  I'm not sure.

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: W(70, relprec = 8) # indirect doctest
        4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
        """
        if self.relprec == relprec:
            return False
        if self.relprec != 0:
            ZZ_pX_destruct(&self.unit)
        self.prime_pow.restore_context_capdiv(relprec)
        if relprec != 0:
            ZZ_pX_construct(&self.unit)
        self.relprec = relprec
        return True

    cdef bint _set_prec_both(self, long absprec, long relprec):
        """
        Assuming self.ordp is set, sets the relative precision of self to the minimum of abs(relprec) and absprec-self.ordp.

        If relprec is negative, will set self.relprec to be negative (indicating unnormalized unit)

        Returns True iff self.relprec = 0, ie self was set to an inexact zero

        Note that this will wipe out anything in self.unit.  Be careful
        resetting self.unit directly: if you set it to a different modulus, NTL may
        have problems.  The safest way to reset self.unit to a different modulus is:
        self.prime_pow.restore_context_capdiv(self.relprec)
        cdef ZZ_pX_c tmp = self.unit
        self._set_prec_rel(new_rel_prec)
        ZZ_pX_conv_modulus(self.unit, tmp, self.prime_pow.get_context_capdiv(self.relprec).x)

        You may be able to just set self.relprec and ZZ_pX_conv_modulus if you're decreasing
        precision.  I'm not sure.

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: W(70, 8) # indirect doctest
        4*w^5 + 3*w^7 + O(w^8)
        """
        if self.relprec != 0:
            ZZ_pX_destruct(&self.unit)
        self.relprec = absprec - self.ordp
        cdef long arelprec
        if relprec < 0:
            arelprec = -relprec
        else:
            arelprec = relprec
        if self.relprec <= 0:
            self._set_inexact_zero(absprec)
        else:
            if arelprec < self.relprec:
                self.relprec = arelprec
            self.prime_pow.restore_context_capdiv(self.relprec)
            ZZ_pX_construct(&self.unit)
            if relprec < 0:
                self.relprec = -self.relprec
        return self.relprec == 0

    cdef int _normalize(self) except -1:
        """
        Normalizes self, adjusting self.ordp, self.relprec, and self.unit so that self.unit actually represents a unit.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1+w)^5
        sage: y = z - 1
        sage: y._ntl_rep_unnormalized()
        [5 3005 25 3060 5]
        sage: y # indirect doctest
        w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        sage: y._ntl_rep_unnormalized()
        [41 26 152 49 535]
        """
        cdef long minval, mini, shift
        if self.relprec < 0:
            if ZZ_pX_IsZero(self.unit):
                self.ordp -= self.relprec # note that self.relprec < 0
                self.relprec = 0
                ZZ_pX_destruct(&self.unit)
            else:
                ZZ_pX_min_val_coeff(minval, mini, self.unit, self.prime_pow.pow_ZZ_tmp(1)[0])
                if self.prime_pow.e == 1:
                    shift = minval
                else:
                    shift = minval * self.prime_pow.e + mini
                if shift >= -self.relprec:
                    self.ordp -= self.relprec # note that self.relprec < 0
                    self.relprec = 0
                    ZZ_pX_destruct(&self.unit)
                elif shift > 0:
                    self.relprec = -self.relprec - shift
                    self.ordp += shift
                    self._internal_lshift(-shift)
                else:
                    self.relprec = -self.relprec

    cdef int _internal_lshift(self, long shift) except -1:
        """
        Multiplies self.unit by x^shift.

        Note that self.relprec must be set before calling this function, and self.unit must be defined to precision self.relprec + shift
        This function does not alter self.ordp even though it WILL change the valuation of self.unit
        Also note that if you call this function you should usually manually set self.relprec = -self.relprec since this function will usually unnormalize self.

        TESTS:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1+w)^5
        sage: y = z - 1
        sage: y # indirect doctest
        w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        """
        cdef ZZ_pX_c tmpP
        cdef ZZ_pX_Modulus_c mod
        if self.prime_pow.e == 1:
            if shift > 0:
                ZZ_pX_left_pshift(self.unit, self.unit, self.prime_pow.pow_ZZ_tmp(shift)[0], self.prime_pow.get_context(self.relprec).x)
            else:
                ZZ_pX_right_pshift(self.unit, self.unit, self.prime_pow.pow_ZZ_tmp(-shift)[0], self.prime_pow.get_context(self.relprec).x)
        else:
            if shift > 0:
                self.prime_pow.restore_context_capdiv(self.relprec)
                mod = self.prime_pow.get_modulus_capdiv(self.relprec)[0]
                ZZ_pX_PowerXMod_long_pre(tmpP, shift, mod)
                ZZ_pX_MulMod_pre(self.unit, self.unit, tmpP, mod)
            elif shift < 0:
                self.prime_pow.eis_shift_capdiv(&self.unit, &self.unit, -shift, self.relprec)

    cdef int _pshift_self(self, long shift) except -1:
        """
        Multiplies self by p^shift.

        This function assumes that self.relprec, self.ordp and self.unit are already set (in the case self.prime_pow.e != 1),
        and is more reasonable to call externally than _internal_lshift

        EXAMPLES:
        sage: R = Qp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(3/25, relprec = 6); z
        3*w^-10 + 3*w^-8 + 2*w^-6 + O(w^-4)
        sage: z * 25
        3 + O(w^6)
        """
        cdef ZZ_pX_c high_shifter, high_shifter2
        cdef ZZ_pX_Modulus_c modulus, modulus_up
        cdef ntl_ZZ_pContext_class c
        cdef PowComputer_ZZ_pX_small_Eis sm
        cdef PowComputer_ZZ_pX_big_Eis big
        cdef ntl_ZZ_pX printer
        cdef ZZ_pX_c* high_array
        cdef long i, high_length
        if self.prime_pow.e == 1:
            self.ordp += shift
        else:
            self.ordp += shift * self.prime_pow.e
            if shift < 0:
                shift = -shift
                c = self.prime_pow.get_context_capdiv(self.relprec)
                c.restore_c()
                modulus = self.prime_pow.get_modulus_capdiv(self.relprec)[0]
                if PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_big_Eis):
                    high_array = (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_shifter
                elif PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_small_Eis):
                    high_array = (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_shifter
                else:
                    raise RuntimeError, "unrecognized PowComputer type"
                ZZ_pX_conv_modulus(high_shifter, high_array[0], c.x)
                ZZ_pX_InvMod_newton_ram(high_shifter, high_shifter, modulus, c.x)
                ZZ_pX_PowerMod_long_pre(high_shifter, high_shifter, shift, modulus)
                ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus)

                #modulus_up = self.prime_pow.get_modulus_capdiv(self.relprec + self.prime_pow.e)[0]
                #c = self.prime_pow.get_context_capdiv(self.relprec + self.prime_pow.e)
                #c.restore_c()
                #ZZ_pX_SetX(high_shifter)
                #ZZ_pX_LeftShift(high_shifter, high_shifter, self.prime_pow.e)
                #ZZ_pX_sub(high_shifter, high_shifter, modulus_up.val())
                #c = self.prime_pow.get_context_capdiv(self.relprec)
                #c.restore_c()
                #ZZ_pX_right_pshift(high_shifter2, high_shifter, self.prime_pow.pow_ZZ_tmp(1)[0], c.x)
                #printer = ntl_ZZ_pX([],c)
                #printer.x = high_shifter2
                #print "high_shifter2 = %s"%(printer)
                #print "shift = %s"%shift
                #printer.x = modulus.val()
                #print "modulus = %s"%(printer)
                #print c
                #print "before PowerMod1"
                #ZZ_pX_PowerMod_long_pre(high_shifter2, high_shifter2, shift, modulus)
                #print "after PowerMod"
                #ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter2, modulus)
            elif shift > 0:
                i = 0
                c = self.prime_pow.get_context_capdiv(self.relprec)
                c.restore_c()
                modulus = self.prime_pow.get_modulus_capdiv(self.relprec)[0]
                if PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_big_Eis):
                    high_array = (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_shifter
                    high_length = (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_length
                elif PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_small_Eis):
                    high_array = (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_shifter
                    high_length = (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_length
                else:
                    raise RuntimeError, "unrecognized PowComputer type"
                if shift >= self.prime_pow.prec_cap:
                    # high_shifter = p^(2^(high_length - 1))/x^(e*2^(high_length - 1))
                    ZZ_pX_conv_modulus(high_shifter, high_array[high_length-1], c.x)
                    # if shift = r + s * 2^(high_length - 1)
                    # then high_shifter = p^(s*2^(high_length - 1))/x^(e*s*2^(high_length - 1))
                    ZZ_pX_PowerMod_long_pre(high_shifter, high_shifter, (shift / (1 << (high_length - 1))), modulus)
                    ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus)
                    # Now we only need to multiply self.unit by p^r/x^(e*r) where r < 2^(high_length - 1), which is tractible.
                    shift = shift % (1 << (high_length - 1))
                while shift > 0:
                    if shift & 1:
                        ZZ_pX_conv_modulus(high_shifter, high_array[i], c.x)
                        ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus)
                    shift = shift >> 1
                    i += 1

    def __dealloc__(self):
        """
        Deallocates self.unit if needed.

        EXAMPLES:
        sage: R = Qp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = W(3/25, relprec = 6); z
        3*w^-10 + 3*w^-8 + 2*w^-6 + O(w^-4)
        sage: del z #indirect doctest
        """
        if self.relprec != 0:
            ZZ_pX_destruct(&self.unit)

    cdef pAdicZZpXCRElement _new_c(self, long relprec):
        """
        Returns a new element with the same parent as self and relative precision relprec

        Note that if relprec is non-positive, the convention is that relprec = 0 indicates an exact or inexact zero,
        relprec < 0 indicates an unnormalized element.

        EXAMPLES:
        sage: R = Qp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: w^5 + 1 # indirect doctest
        1 + w^5 + O(w^25)
        """
        cdef pAdicZZpXCRElement ans = PY_NEW(pAdicZZpXCRElement)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        if relprec > 0:
            self.prime_pow.restore_context_capdiv(relprec)
            ans.relprec = relprec
            ZZ_pX_construct(&ans.unit)
        elif relprec == 0:
            ans._set_exact_zero()
        else:
            self.prime_pow.restore_context_capdiv(-relprec)
            ans.relprec = relprec
            ZZ_pX_construct(&ans.unit)
        return ans

    def __reduce__(self):
        """
        Pickles self.

        EXAMPLES:
        sage: R = Qp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5 - 1
        sage: loads(dumps(z)) == z
        True
        """
        self.prime_pow.restore_context_capdiv(self.relprec)
        cdef ntl_ZZ_pX holder = PY_NEW(ntl_ZZ_pX)
        holder.c = self.prime_pow.get_context_capdiv(self.relprec)
        holder.x = self.unit
        cdef Integer relprec, ordp
        relprec = PY_NEW(Integer)
        ordp = PY_NEW(Integer)
        mpz_set_si(relprec.value, self.relprec)
        mpz_set_si(ordp.value, self.ordp)
        return make_ZZpXCRElement, (self.parent(), holder, ordp, relprec, 0)

    cdef int _cmp_units(left, pAdicGenericElement right) except -2:
        """
        For units left and right, returns 0 if they are equal up to the
        lesser of the two precisions, or 1 if they are not.

        EXAMPLES:
        sage: R = Qp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: w == 1 # indirect doctest
        False
        sage: y = 1 + w + O(w^7)
        sage: z = 1 + w + w^10 + O(w^13)
        sage: y == z
        True
        """
        # This function needs improvement.  In particular, there are a lot of
        # speed improvements to be had, and it should be changed so that it
        # returns 1 only half the time (and -1 the other half) when left and
        # right are not equal.
        cdef pAdicZZpXCRElement diff = <pAdicZZpXCRElement> (left - right)
        diff._normalize()
        if diff.relprec == 0:
            return 0
        # for now, just return 1
        return 1

    def __invert__(self):
        """
        Returns the inverse of self.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5
        sage: y = ~z; y # indirect doctest
        1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14 + 3*w^16 + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24 + O(w^25)
        sage: y.parent()
        Eisenstein Extension of 5-adic Field with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (3*5^2 + O(5^7))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^8))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
        sage: z = z - 1
        sage: ~z
        w^-5 + 4*w^-4 + 4*w^-3 + 4*w^-2 + 2*w^-1 + 1 + w + 4*w^2 + 4*w^3 + 4*w^4 + w^5 + w^6 + w^7 + 4*w^8 + 4*w^9 + 2*w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^14 + O(w^15)
        """
        return self._invert_c_impl()

    cdef RingElement _invert_c_impl(self):
        """
        Returns the inverse of self.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5
        sage: y = ~z; y # indirect doctest
        1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14 + 3*w^16 + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24 + O(w^25)
        sage: y.parent()
        Eisenstein Extension of 5-adic Field with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (3*5^2 + O(5^7))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^8))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
        sage: z = z - 1
        sage: ~z
        w^-5 + 4*w^-4 + 4*w^-3 + 4*w^-2 + 2*w^-1 + 1 + w + 4*w^2 + 4*w^3 + 4*w^4 + w^5 + w^6 + w^7 + 4*w^8 + 4*w^9 + 2*w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^14 + O(w^15)
        """
        if self._is_exact_zero():
            raise ZeroDivisionError, "cannot divide by zero"
        if self._is_inexact_zero(): # this calls _normalize
            raise PrecisionError, "cannot divide by something indistinguishable from zero"
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        if not ans.prime_pow.in_field:
            ans._parent = self._parent.fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        ans.ordp = -self.ordp
        _sig_on
        if self.prime_pow.e == 1:
            ZZ_pX_InvMod_newton_unram(ans.unit, self.unit, self.prime_pow.get_modulus(ans.relprec)[0], self.prime_pow.get_context(ans.relprec).x, self.prime_pow.get_context(1).x)
        else:
            ZZ_pX_InvMod_newton_ram(ans.unit, self.unit, self.prime_pow.get_modulus_capdiv(ans.relprec)[0], self.prime_pow.get_context_capdiv(ans.relprec).x)
        _sig_off
        return ans

    cdef pAdicZZpXCRElement _lshift_c(self, long n):
        """
        Multiplies self by the uniformizer raised to the power n.  If n is negative, right shifts by -n.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5
        sage: z
        1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        sage: z << 17 # indirect doctest
        w^17 + w^22 + w^23 + 2*w^24 + 4*w^25 + 3*w^27 + w^29 + 4*w^30 + 4*w^31 + 4*w^32 + 4*w^33 + 4*w^34 + 4*w^37 + w^38 + 4*w^41 + O(w^42)
        sage: z << (-1)
        w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
        """
        if not self.prime_pow.in_field and n < -self.ordp:
            return self._rshift_c(-n)
        if n >= maxordp or n <= -maxordp:
            raise ValueError, "overflow in valuation"
        cdef pAdicZZpXCRElement ans
        if self._is_exact_zero() or n == 0:
            return self
        elif self._is_inexact_zero():
            ans = self._new_c(0)
        else:
            ans = self._new_c(self.relprec)
            ans.unit = self.unit
        ans.ordp = self.ordp + n
        if ans.ordp >= maxordp or n <= -maxordp:
            raise ValueError, "overflow in valuation"
        return ans

    def __lshift__(pAdicZZpXCRElement self, shift):
        """
        Multiplies self by the uniformizer raised to the power n.  If n is negative, right shifts by -n.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5
        sage: z
        1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        sage: z << 17 # indirect doctest
        w^17 + w^22 + w^23 + 2*w^24 + 4*w^25 + 3*w^27 + w^29 + 4*w^30 + 4*w^31 + 4*w^32 + 4*w^33 + 4*w^34 + 4*w^37 + w^38 + 4*w^41 + O(w^42)
        sage: z << (-1)
        w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
        """
        cdef pAdicZZpXCRElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            if self._is_exact_zero():
                return self
            if self.prime_pow.in_field or mpz_sgn((<Integer>shift).value) > 0:
                raise ValueError, "Shift does not fit in long"
            else:
                ans = self._new_c(0)
                ans.ordp = 0
                return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicZZpXCRElement _rshift_c(self, long n):
        """
        Divides self by the uniformizer raised to the power n.  If parent is not a field, throws away
        the non-positive part of the series expansion.
        If n is negative, left shifts by -n.

        EXAMPLES:
        sage: R = Zp(5,5,print_mode='digits')
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5
        sage: for m in range(26): repr(z >> m) # indirect doctest
        '...4001400444441030421100001'
        '...400140044444103042110000'
        '...40014004444410304211000'
        '...4001400444441030421100'
        '...400140044444103042110'
        '...40014004444410304211'
        '...4001400444441030421'
        '...400140044444103042'
        '...40014004444410304'
        '...4001400444441030'
        '...400140044444103'
        '...40014004444410'
        '...4001400444441'
        '...400140044444'
        '...40014004444'
        '...4001400444'
        '...400140044'
        '...40014004'
        '...4001400'
        '...400140'
        '...40014'
        '...4001'
        '...400'
        '...40'
        '...4'
        '...'
        sage: repr(z >> (-4))
        '...40014004444410304211000010000'
        """
        if self.prime_pow.in_field or n <= self.ordp:
            return self._lshift_c(-n)
        if self._is_exact_zero() or n == 0:
            return self
        cdef long arelprec
        if self.relprec < 0:
            arelprec = -self.relprec
        else:
            arelprec = self.relprec
        cdef pAdicZZpXCRElement ans
        if arelprec > n - self.ordp:
            ans = self._new_c(arelprec - (n - self.ordp))
            if self.prime_pow.e == 1:
                ZZ_pX_right_pshift(ans.unit, self.unit, self.prime_pow.pow_ZZ_tmp(n - self.ordp)[0], self.prime_pow.get_context(ans.relprec).x)
            else:
                self.prime_pow.eis_shift_capdiv(&ans.unit, &self.unit, n - self.ordp, ans.relprec)
        else:
            ans = self._new_c(0)
        ans.ordp = 0
        ans.relprec = -ans.relprec
        return ans

    def __rshift__(pAdicZZpXCRElement self, shift):
        """
        Divides self by the uniformizer raised to the power n.  If parent is not a field, throws away the non-positive part of the series expansion.
        If n is negative, left shifts by -n.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5
        sage: z
        1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        sage: z >> (6) # indirect doctest
        1 + 2*w + 4*w^2 + 3*w^4 + w^6 + 4*w^7 + 4*w^8 + 4*w^9 + 4*w^10 + 4*w^11 + 4*w^14 + w^15 + 4*w^18 + O(w^19)
        sage: z >> (-4)
        w^4 + w^9 + w^10 + 2*w^11 + 4*w^12 + 3*w^14 + w^16 + 4*w^17 + 4*w^18 + 4*w^19 + 4*w^20 + 4*w^21 + 4*w^24 + w^25 + 4*w^28 + O(w^29)
        sage: F = W.fraction_field()
        sage: z = F(z)
        sage: z >> 7
        w^-7 + w^-2 + w^-1 + 2 + 4*w + 3*w^3 + w^5 + 4*w^6 + 4*w^7 + 4*w^8 + 4*w^9 + 4*w^10 + 4*w^13 + w^14 + 4*w^17 + O(w^18)
        """
        cdef pAdicZZpXCRElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            if self._is_exact_zero():
                return self
            if self.prime_pow.in_field or mpz_sgn((<Integer>shift).value) < 0:
                raise ValueError, "Shift does not fit in long"
            else:
                ans = self._new_c(0)
                ans.ordp = 0
                return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cdef ModuleElement _neg_c_impl(self):
        """
        Returns -self.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5; z
        1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        sage: -z # indirect doctest
        4 + 3*w^5 + 4*w^6 + w^7 + w^8 + w^9 + w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^15 + 3*w^16 + w^17 + 2*w^18 + 3*w^19 + 2*w^21 + 4*w^23 + 4*w^24 + O(w^25)
        sage: y = z + (-z); y
        O(w^25)
        sage: -y
        O(w^25)
        sage: -W(0)
        0
        """
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        ans.ordp = self.ordp
        if self.relprec != 0:
            self.prime_pow.restore_context_capdiv(self.relprec)
            ZZ_pX_negate(ans.unit, self.unit)
        return ans

    def __pow__(pAdicZZpXCRElement self, _right, m): # m ignored
        r"""
        Computes self^right.

        Note: when right is divisible by p then one can get more precision than expected.
        Lemma 2.1 (Constructing Class Fields over Local Fields, Sebastian Pauli):
        Let $\alpha$ be in $\mathcal{O}_K$.  Let $p = -\pi_K^{e_K} \epsilon$ be the factorization
        of $p$ where $\epsilon$ is a unit.  Then the $p$-th power of $1 + \alpha \pi_K^{\lambda}$ satisifes
                                             / 1 + \alpha^p \pi_K^{p \lambda}                      mod \mathfrak{p}_K^{p \lambda + 1}   if 1 \le \lambda < \frac{e_K}{p-1}
        (1 + \alpha \pi^{\lambda})^p \equiv {  1 + (\alpha^p - \epsilon \alpha) \pi_K^{p \lambda}  mod \mathfrak{p}_K^{p \lambda + 1}   if \lambda = \frac{e_K}{p-1}
                                             \ 1 - \epsilon \alpha \pi_K^{\lambda + e}             mod \mathfrak{p}_K^{\lambda + e + 1} if \lambda > \frac{e_K}{p-1}

        So if right is divisible by $p^k$ we can multiply the relative precision by $p$ until we exceed $e/(p-1)$, then add $e$ until we have done a total of $k$ things:
        the precision of the result can therefore be greater than the precision of self.

        There is also the issue of $p$-adic exponents, and determining how the precision of the exponent affects the precision of the result.
        In computing $(a + O(\pi^k))^{b + O(p^m)}$, one needs that the reduction of $a$ mod $\pi$ is in the prime field $F_p$ (so that the $p^m$ power of the Teichmuller
        part is constant as $m$ increases).  Given this restriction, we can factor out the Teichmuller part and use the above lemma to find the first spot where
        $(1 + \alpha \pi^{\lambda})^(p^m)$ differs from 1.  We compare this with the precision bound given by computing $(a + O(\pi^k))^b$ and take the lesser of the two.

        In order to do this we need to compute the valuation of (self / self.parent().teichmuller(self)) - 1.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: (1 + w)^5 # indirect doctest
        1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        sage: (1 + w + O(w^19))^5
        1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + O(w^24)
        sage: (1 + O(w))^5
        1 + O(w^5)
        sage: (1 + w + O(w^3))^25
        1 + w^10 + w^11 + 4*w^12 + O(w^13)
        sage: (3 + 2*w + w^2 + O(w^6))^(15 + O(125))
        2 + 4*w^6 + w^7 + 3*w^8 + 3*w^9 + 4*w^10 + O(w^11)
        sage: (3 + 2*w + w^2 + O(w^6))^(15 + O(25))
        2 + 4*w^6 + w^7 + 3*w^8 + 3*w^9 + O(w^10)
        sage: (3 + w^2 + O(w^6))^(15+O(25))
        2 + w^5 + 4*w^7 + w^9 + 3*w^10 + O(w^11)
        sage: R = Zp(2, 10)
        sage: S.<x> = R[]
        sage: f = x^34 + 18*x^5 - 72*x^3 + 2
        sage: W.<w> = R.ext(f)
        sage: (1+w+O(w^2))^8
        1 + w^8 + O(w^16)
        sage: (1+w+O(w^2))^16
        1 + w^16 + O(w^32)
        sage: (1+w+O(w^2))^32
        1 + w^32 + w^50 + w^55 + w^60 + O(w^64)
        sage: (1+w+O(w^2))^64
        1 + w^64 + w^66 + w^71 + w^76 + w^81 + w^84 + w^86 + w^91 + w^94 + w^96 + O(w^98)
        """
        self._normalize()
        cdef Integer right
        cdef bint padic_exp
        cdef long exp_prec
        cdef long exp_val
        cdef long relprec
        cdef long threshold # e / (p-1)
        cdef long prime_long
        cdef mpz_t tmp, tmp2
        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            threshold = 0
        else:
            threshold = self.prime_pow.e / (mpz_get_si(self.prime_pow.prime.value) - 1)
        cdef Integer base_level
        cdef pAdicZZpXCRElement ans
        cdef long i
        if self._is_exact_zero():
            # Return 0 except for 0^0 error or type error on the exponent.
            if PY_TYPE_CHECK(_right, Integer) or PY_TYPE_CHECK(_right, Rational) or (PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime)  or isinstance(_right, (int, long)):
                if _right == 0:
                    raise ArithmeticError, "0^0 is undefined"
                return self
            else:
                raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        elif self._is_inexact_zero():
            # If an integer exponent, return an inexact zero of valuation right * self.ordp.  Otherwise raise an error.
            if isinstance(_right, (int, long)):
                _right = Integer(_right)
            if PY_TYPE_CHECK(_right, Integer):
                ans = self._new_c(0)
                mpz_init_set_si(tmp, self.ordp)
                mpz_mul(tmp, tmp, (<Integer>_right).value)
                if mpz_cmp_si(tmp, maxordp) >= 0 or mpz_cmp_si(tmp, -maxordp) <= 0:
                    raise ValueError, "valuation overflow"
                ans.ordp = mpz_get_si(tmp)
                mpz_clear(tmp)
                return ans
            elif PY_TYPE_CHECK(_right, Rational) or (PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime):
                raise ValueError, "Need more precision"
            else:
                raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        if isinstance(_right, (int, long)):
            _right = Integer(_right)
        if PY_TYPE_CHECK(_right, Integer):
            right = <Integer> _right
            if right == 0:
                # return 1 to maximum precision
                ans = self._new_c(self.prime_pow.ram_prec_cap)
                ans.ordp = 0
                ZZ_pX_SetCoeff_long(ans.unit, 0, 1)
                return ans
            padic_exp = False
            exp_val = _right.valuation(self.prime_pow.prime) ##
        elif PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime:
            if self.ordp != 0:
                raise ValueError, "in order to raise to a p-adic exponent, base must be a unit"
            right = Integer(_right)
            padic_exp = True
            exp_prec = _right.precision_absolute() ##
            exp_val = _right.valuation() ##
            if exp_val < 0:
                raise NotImplementedError, "negative valuation exponents not yet supported"
            # checks to see if the residue of self.unit is in the prime field.
            if self.prime_pow.e == 1:
                for i from 1 <= i <= ZZ_pX_deg(self.unit):
                    if not ZZ_divide_test(ZZ_p_rep(ZZ_pX_coeff(self.unit, i)), self.prime_pow.pow_ZZ_tmp(1)[0]):
                        raise ValueError, "in order to raise to a p-adic exponent, base must reduce to an element of F_p mod the uniformizer"
            # compute the "level"
            teich_part = self.parent().teichmuller(self)
            base_level = (self / teich_part - 1).valuation() ##
        elif PY_TYPE_CHECK(_right, Rational):
            raise NotImplementedError
        else:
            raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        # Now we compute the increased relprec due to the exponent having positive p-adic valuation
        if exp_val > 0:
            mpz_init_set_si(tmp, self.relprec)
            while mpz_cmp_si(tmp, threshold) <= 0 and exp_val > 0:
                mpz_mul(tmp, tmp, self.prime_pow.prime.value)
                exp_val -= 1
            if exp_val > 0:
                mpz_init_set_si(tmp2, self.prime_pow.e)
                mpz_addmul_ui(tmp, tmp2, exp_val)
                mpz_clear(tmp2)
            if mpz_cmp_si(tmp, self.prime_pow.ram_prec_cap) > 0:
                relprec = self.prime_pow.ram_prec_cap
            else:
                relprec = mpz_get_si(tmp)
            mpz_clear(tmp)
        else:
            relprec = self.relprec
        # Now we compute the limit on relprec due to a non-infinite precision on the exponent.
        if padic_exp:
            if exp_prec > 0:
                # I can freely change base_level, so I use it in place of tmp above.
                while mpz_cmp_si(base_level.value, threshold) <= 0 and exp_prec > 0:
                    mpz_mul(base_level.value, base_level.value, self.prime_pow.prime.value)
                    exp_prec -= 1
                if exp_prec > 0:
                    mpz_init_set_si(tmp2, self.prime_pow.e)
                    mpz_addmul_ui(base_level.value, tmp2, exp_prec)
                    mpz_clear(tmp2)
                if mpz_cmp_si(base_level.value, relprec) < 0:
                    relprec = mpz_get_si(base_level.value)
            else:
                ans = self._new_c(0)
                ans.ordp = 0
                return ans
        ans = self._new_c(relprec)
        self.prime_pow.restore_context_capdiv(relprec)
        if self.ordp == 0:
            ans.ordp = 0
        else:
            mpz_init_set(tmp, right.value)
            mpz_mul_si(tmp, tmp, self.ordp)
            if mpz_cmp_si(tmp, maxordp) >= 0 or mpz_cmp_si(tmp, -maxordp) <= 0:
                raise ValueError, "valuation overflow"
            ans.ordp = mpz_get_si(tmp)
            mpz_clear(tmp)
        cdef ntl_ZZ rZZ = PY_NEW(ntl_ZZ)
        mpz_to_ZZ(&rZZ.x, &right.value)
        _sig_on
        ZZ_pX_PowerMod_pre(ans.unit, self.unit, rZZ.x, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
        _sig_off
        return ans

    cdef ModuleElement _add_c_impl(self, ModuleElement _right):
        """
        Computes the sum of self and right.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)) - 69 # indirect doctest
        1 + O(w^13)
        sage: -69 + (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13))
        1 + O(w^13)
        sage: y = (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13))
        sage: y - 70
        O(w^13)
        sage: y + 0
        4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
        """
        cdef pAdicZZpXCRElement right = <pAdicZZpXCRElement>_right
        cdef pAdicZZpXCRElement ans
        cdef long tmpL
        cdef ZZ_pX_c tmpP
        if self.relprec == 0:
            if self.ordp >= right.ordp + right.relprec: # or self._is_exact_zero()
                return right
            elif self.ordp <= right.ordp:
                ans = self._new_c(0)
                ans.ordp = self.ordp
            else:
                ans = self._new_c(self.ordp - right.ordp)
                ZZ_pX_conv_modulus(ans.unit, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                if right.relprec < 0:
                    ans.relprec = -ans.relprec
                ans.ordp = right.ordp
            return ans
        if right.relprec == 0:
            if right.ordp >= self.ordp + self.relprec: # or right._is_exact_zero()
                return self
            elif right.ordp <= self.ordp:
                ans = self._new_c(0)
                ans.ordp = right.ordp
            else:
                ans = self._new_c(right.ordp - self.ordp)
                ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                if self.relprec < 0:
                    ans.relprec = -ans.relprec
                ans.ordp = self.ordp
            return ans
        cdef long srprec = self.relprec
        if srprec < 0:
            srprec = -srprec
        cdef long rrprec = right.relprec
        if rrprec < 0:
            rrprec = -rrprec
        #print "self.ordp = %s\nright.ordp = %s"%(self.ordp, right.ordp)
        #print "self = %s\nright = %s"%(self, right)
        if self.ordp == right.ordp:
            # The relative precision of the sum is the minimum of the relative precisions in this case, possibly decreasing if we got cancellation
            # Since the valuations are the same, we could just add the units, if they had the same modulus.
            # But they don't necessarily, so we may have to conv_modulus
            if srprec == rrprec:
                ans = self._new_c(-srprec) # -srprec indicates that ans is not normalized
                self.prime_pow.restore_context_capdiv(srprec)
                ZZ_pX_add(ans.unit, self.unit, right.unit)
            elif srprec < rrprec:
                ans = self._new_c(-srprec)
                ZZ_pX_conv_modulus(ans.unit, right.unit, self.prime_pow.get_context_capdiv(srprec).x)
                self.prime_pow.restore_context_capdiv(srprec)
                # conv_modulus should have restored the context, so we don't need to again.
                ZZ_pX_add(ans.unit, ans.unit, self.unit)
            else:
                ans = self._new_c(-rrprec)
                ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(rrprec).x)
                self.prime_pow.restore_context_capdiv(rrprec)
                # conv_modulus should have restored the context, so we don't need to again.
                ZZ_pX_add(ans.unit, ans.unit, right.unit)
            ans.ordp = self.ordp
        elif self.ordp < right.ordp:
            tmpL = right.ordp - self.ordp
            if tmpL >= srprec:
                return self
            if srprec <= tmpL + rrprec:
                ans = self._new_c(-srprec)
            else:
                ans = self._new_c(-tmpL - rrprec)
            ans.ordp = self.ordp
            ZZ_pX_conv_modulus(ans.unit, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ans._internal_lshift(tmpL)
            if srprec <= tmpL + rrprec:
                ZZ_pX_add(ans.unit, ans.unit, self.unit)
            else:
                ZZ_pX_conv_modulus(tmpP, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                ZZ_pX_add(ans.unit, ans.unit, tmpP)
            # if self is normalized, then the valuations are actually different so the sum will be normalized.
            if self.relprec > 0:
                ans.relprec = -ans.relprec
        else:
            tmpL = self.ordp - right.ordp
            if tmpL >= rrprec:
                return right
            if rrprec <= tmpL + srprec:
                ans = self._new_c(-rrprec)
            else:
                ans = self._new_c(-tmpL - srprec)
            ans.ordp = right.ordp
            ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ans._internal_lshift(tmpL)
            if rrprec <= tmpL + srprec:
                ZZ_pX_add(ans.unit, ans.unit, right.unit)
            else:
                ZZ_pX_conv_modulus(tmpP, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                ZZ_pX_add(ans.unit, ans.unit, tmpP)
            # if right is normalized, then the valuations are actually different so the sum will be normalized.
            if right.relprec > 0:
                ans.relprec = -ans.relprec
        return ans

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Returns the difference of self and right.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(329)
        sage: b = W(111)
        sage: a - b
        3 + 3*w^5 + w^7 + 2*w^9 + 3*w^10 + 4*w^11 + 2*w^13 + 2*w^14 + w^15 + 4*w^16 + 2*w^18 + 3*w^19 + 2*w^20 + 3*w^21 + w^22 + w^24 + O(w^25)
        sage: W(218)
        3 + 3*w^5 + w^7 + 2*w^9 + 3*w^10 + 4*w^11 + 2*w^13 + 2*w^14 + w^15 + 4*w^16 + 2*w^18 + 3*w^19 + 2*w^20 + 3*w^21 + w^22 + w^24 + O(w^25)
        sage: a - O(w^14)
        4 + 3*w^10 + 2*w^12 + O(w^14)
        sage: a - 0
        4 + 3*w^10 + 2*w^12 + w^14 + 2*w^15 + w^16 + 3*w^17 + 3*w^18 + w^19 + 2*w^21 + 4*w^22 + w^23 + 4*w^24 + O(w^25)
        sage: O(w^14) - a
        1 + 4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + O(w^14)
        """
        # For now, a simple implementation
        return self + (-right)

    cdef RingElement _mul_c_impl(self, RingElement _right):
        """
        Returns the product of self and right.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(329)
        sage: b = W(111)
        sage: a*b
        4 + 3*w^5 + w^7 + 2*w^9 + 4*w^11 + 3*w^12 + 2*w^13 + w^14 + 2*w^15 + 3*w^16 + 4*w^17 + 4*w^18 + 2*w^19 + 2*w^21 + 4*w^22 + 2*w^23 + w^24 + O(w^25)
        sage: a * 0
        0
        sage: a * O(w^14)
        O(w^14)
        """
        cdef pAdicZZpXCRElement right = <pAdicZZpXCRElement>_right
        cdef ZZ_pX_c modulus_corrected
        cdef ntl_ZZ_pContext_class ctx
        cdef pAdicZZpXCRElement ans
        if self._is_exact_zero():
            return self
        if right._is_exact_zero():
            return right
        self._normalize()
        right._normalize()
        if self.relprec <= right.relprec:
            ans = self._new_c(self.relprec)
        else:
            ans = self._new_c(right.relprec)
        ans.ordp = self.ordp + right.ordp
        if ans.ordp > maxordp or ans.ordp < -maxordp:
            raise ValueError, "valuation overflow"
        if ans.relprec == 0:
            return ans
        if self.relprec == right.relprec:
            self.prime_pow.restore_context_capdiv(ans.relprec)
            _sig_on
            ZZ_pX_MulMod_pre(ans.unit, self.unit, right.unit, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
            _sig_off
        elif self.relprec < right.relprec:
            _sig_on
            ZZ_pX_conv_modulus(modulus_corrected, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ZZ_pX_MulMod_pre(ans.unit, self.unit, modulus_corrected, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
            _sig_off
        else:
            _sig_on
            ZZ_pX_conv_modulus(modulus_corrected, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ZZ_pX_MulMod_pre(ans.unit, right.unit, modulus_corrected, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
            _sig_off
        return ans

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Returns the quotient of self by right.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: W(14) / W(125)
        4*w^-15 + w^-13 + 3*w^-11 + 2*w^-10 + 3*w^-9 + 4*w^-8 + 4*w^-7 + 3*w^-6 + 2*w^-5 + 4*w^-4 + 3*w^-3 + 2*w^-2 + 4*w^-1 + 2 + w^2 + w^4 + 4*w^5 + w^6 + w^7 + 3*w^9 + O(w^10)
        """
        # for now, a simple implementation
        return self * (~right)

    def copy(self):
        """
        Returns a copy of self.

        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: b = W(45, 17); b
        4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16 + O(w^17)
        sage: c = b.copy(); c
        4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16 + O(w^17)
        sage: c is b
        False
        """
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        ans.ordp = self.ordp
        ans.unit = self.unit
        return ans

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def is_zero(self, absprec = None):
        """
        Returns whether the valuation of self is at least absprec.  If absprec is None,
        returns if self is indistinugishable from zero.

        If self is an inexact zero of valuation less than absprec, raises a PrecisionError.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: O(w^189).is_zero()
        True
        sage: W(0).is_zero()
        True
        sage: a = W(675)
        sage: a.is_zero()
        False
        sage: a.is_zero(7)
        True
        sage: a.is_zero(21)
        False
        """
        cdef bint ans
        cdef long aprec
        self._normalize()
        if self._is_exact_zero():
            ans = True
        elif absprec is None:
            ans = (self.relprec == 0)
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0:
                    ans = True
                elif self.relprec == 0:
                    raise PrecisionError, "Not enough precision to determine if element is zero"
                else:
                    ans = False
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
                if self.relprec == 0 and aprec > self.ordp:
                    raise PrecisionError, "Not enough precision to determine if element is zero"
                else:
                    ans = (self.ordp >= aprec)
        return ans

    cpdef ntl_ZZ_pX _ntl_rep_unnormalized(self):
        """
        Returns an ntl_ZZ_pX holding the current unit part of self.

        self is not normalized before this, so the polynomial returned
        may not actually be a unit.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(566); b = W(209)
        sage: c = a + b; c._ntl_rep_unnormalized()
        [775]
        sage: c
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24 + O(w^25)
        sage: c._ntl_rep_unnormalized()
        [106 60 114 35 112]
        """
        if self.relprec == 0:
            raise ValueError, "self == 0"
        self.prime_pow.restore_context_capdiv(self.relprec)
        cdef ntl_ZZ_pX ans = PY_NEW(ntl_ZZ_pX)
        ans.c = self.prime_pow.get_context_capdiv(self.relprec)
        ans.x = self.unit
        return ans

    cpdef ntl_ZZ_pX _ntl_rep(self):
        """
        Returns an ntl_ZZ_pX that holds the unit part of self.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(566); b = W(209)
        sage: c = a + b; c._ntl_rep()
        [106 60 114 35 112]
        sage: c
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24 + O(w^25)
        sage: c._ntl_rep()
        [106 60 114 35 112]
        """
        self._normalize()
        return self._ntl_rep_unnormalized()

    cpdef _ntl_rep_abs(self):
        """
        Returns a pair (f, k) where f is an ntl_ZZ_pX and k is a non-positive integer
        such that self = f(self.parent.gen())*p^k

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(566); b = W(209)
        sage: c = a + b; c._ntl_rep_abs()
        ([775], 0)
        sage: c
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24 + O(w^25)
        sage: c._ntl_rep_abs()
        ([775], 0)
        sage: (~c)._ntl_rep_abs()
        ([121], -2)
        sage: ~c
        w^-10 + w^-8 + 4*w^-6 + 4*w^-5 + 3*w^-3 + 4*w^-2 + 3*w^-1 + 4 + 4*w + 2*w^2 + 4*w^3 + 3*w^4 + O(w^5)
        sage: ~c * 25
        1 + 4*w^5 + 3*w^7 + w^9 + 4*w^10 + 2*w^11 + 3*w^12 + w^13 + 4*w^14 + O(w^15)
        sage: W(121)
        1 + 4*w^5 + 3*w^7 + w^9 + 4*w^10 + 2*w^11 + 3*w^12 + w^13 + 4*w^14 + 2*w^16 + 3*w^17 + 3*w^18 + 4*w^19 + 4*w^20 + 3*w^21 + w^22 + w^23 + 4*w^24 + O(w^25)
        """
        self._normalize()
        if self.ordp == 0:
            return self._ntl_rep(), Integer(0)
        cdef ntl_ZZ_pContext_class ctx
        cdef long little_shift, ppow
        if self.ordp > 0:
            ctx = self.prime_pow.get_context_capdiv(self.ordp + self.relprec)
        else:
            little_shift = ((-self.ordp) % self.prime_pow.e)
            if little_shift != 0:
                little_shift = self.prime_pow.e - little_shift
            ctx = self.prime_pow.get_context_capdiv(self.relprec + little_shift)
        ctx.restore_c()
        cdef pAdicZZpXCRElement dummy = PY_NEW(pAdicZZpXCRElement)
        cdef ntl_ZZ_pX ans = PY_NEW(ntl_ZZ_pX)
        cdef Integer ans_k = PY_NEW(Integer)
        dummy.unit = self.unit
        dummy.prime_pow = self.prime_pow
        if self.ordp > 0:
            dummy.relprec = self.ordp + self.relprec
            dummy._internal_lshift(self.ordp)
            ans.x = dummy.unit
        else:
            ppow = (self.ordp - little_shift) / self.prime_pow.e
            mpz_set_si(ans_k.value, ppow)
            dummy.ordp = 0 # _pshift_self wants ordp set
            dummy.relprec = self.relprec + little_shift
            # self = x^(self.prime_pow.e * ppow) * x^(little_shift) * self.unit
            # so we want to _internal_lshift dummy.unit by little_shift
            dummy._internal_lshift(little_shift)
            # and then write
            # self = p^(ppow) * (x^e/p)^(ppow) * dummy.unit
            # so we need to multiply dummy.unit by (p/x^e)^(-ppow) in the Eisenstein case
            # which we can do by _pshift_self
            dummy._pshift_self(-ppow)
            ans.x = dummy.unit
        ans.c = ctx
        return ans, ans_k

    cdef ZZ_p_c _const_term(self):
        """
        Returns the constant term of self.unit.

        Note: this may be divisible by p if self is not normalized.
        """
        return ZZ_pX_ConstTerm(self.unit)

    def is_equal_to(self, right, absprec = None):
        """
        Returns if self is equal to right modulo self.uniformizer()^absprec.

        If absprec is None, returns if self is equal to right modulo the lower of their two precisions.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(47); b = W(47 + 25)
        sage: a.is_equal_to(b)
        False
        sage: a.is_equal_to(b, 7)
        True
        """
        # Should be sped up later
        return (self - right).is_zero(absprec)

    def lift(self):
        raise NotImplementedError

    cpdef pAdicZZpXCRElement lift_to_precision(self, absprec):
        """
        Returns a pAdicZZpXCRElement congruent to self but with absolute precision
        at least absprec.  If setting absprec that high would violate the precision cap,
        sets the precision to maximum possible.  If self is an inexact zero and
        absprec is greater than the maximum allowed valuation, raises an error.

        Note that the new 'digits' will not necessarily be zero.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(345, 17); a
        4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15 + O(w^17)
        sage: b = a.lift_to_precision(19); b
        4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15 + w^17 + 2*w^18 + O(w^19)
        sage: c = a.lift_to_precision(24); c
        4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15 + w^17 + 2*w^18 + 4*w^19 + 4*w^20 + 2*w^21 + 4*w^23 + O(w^24)
        sage: a._ntl_rep()
        [19 35 118 60 121]
        sage: b._ntl_rep()
        [19 35 118 60 121]
        sage: c._ntl_rep()
        [19 35 118 60 121]
        """
        cdef pAdicZZpXCRElement ans
        cdef long aprec, rprec
        self._normalize()
        if self._is_exact_zero():
            return self
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_fits_slong_p((<Integer>absprec).value) == 0:
            if mpz_sgn((<Integer>absprec).value) < 0 or self.relprec == self.prime_pow.ram_prec_cap:
                return self
            else:
                if self.relprec == 0:
                    raise ValueError, "absprec larger than maximum allowable valuation"
                ans = self._new_c(self.prime_pow.ram_prec_cap)
                ans.ordp = self.ordp
                ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_top_context().x)
                return ans
        aprec = mpz_get_si((<Integer>absprec).value)
        if aprec <= self.ordp + self.relprec:
            return self
        if self.relprec == 0:
            if self.ordp >= aprec:
                return self
            elif aprec >= maxordp:
                raise ValueError, "absprec larger than maximum allowable valuation"
            else:
                ans = self._new_c(0)
                ans._set_inexact_zero(aprec)
                return ans
        # Now we're done handling all the special cases.
        rprec = aprec - self.ordp
        if rprec > self.prime_pow.ram_prec_cap:
            rprec = self.prime_pow.ram_prec_cap
        ans = self._new_c(rprec)
        ans.ordp = self.ordp
        ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(rprec).x)
        return ans

    def list(self, lift_mode = 'simple'):
        """
        Returns a list giving a series representation of self.

        If lift_mode == 'simple' or 'smallest', the returned list will consist
        of integers (in the eisenstein case) or a list of lists of integers (in the unramified case).
        self can be reconstructed as a sum of elements of the list times powers of the uniformiser (in the eisenstein case),
        or as a sum of powers of the p times polynomials in the generator (in the unramified case).
        If lift_mode == 'simple', all integers will be in the range [0,p-1], if 'smallest' they will be in the range [(1-p)/2, p/2].

        If lift_mode == 'teichmuller', returns a list of pAdicZZpXCRElements, all of which are Teichmuller representatives
        and such that self is the sum of that list times powers of the uniformizer.

        Note that zeros are truncated from the returned list, so you must use the valuation function to fully reconstruct self.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: y = W(775, 19); y
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
        sage: y.list()
        [1, 0, 4, 0, 2, 1, 2, 4, 1]
        sage: y.list('smallest')
        [1, 0, -1, 0, 2, 1, 2, 0, 1]
        sage: w^10 - w^12 + 2*w^14 + w^15 + 2*w^16 + w^18 + O(w^19)
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
        sage: g = x^3 + 3*x + 3
        sage: A.<a> = R.ext(g)
        sage: y = 75 + 45*a + 1200*a^2; y
        4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
        sage: y.list()
        [[0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1]]
        sage: y.list('smallest')
        [[0, -1], [-2, 2, -2], [1], [0, 0, 2]]
        sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
        4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
        """
        if lift_mode == 'simple':
            return self.ext_p_list(1)
        elif lift_mode == 'smallest':
            return self.ext_p_list(0)
        elif lift_mode == 'teichmuller':
            return self.teichmuller_list()
        else:
            raise ValueError, "lift mode must be one of 'simple', 'smallest' or 'teichmuller'"

    def log(self, branch = None, same_ring = False):
        r"""
        Compute the p-adic logarithm of any unit.
        (See below for normalization.)

        INPUT:
        branch -- pAdicZZpXFMElement (default None).  The log of the uniformizer.
                  If None, then an error is raised if self is not a unit.
        same_ring -- bool or pAdicGeneric (default True).  When e > p it is
                  possible (even common) that the image of the log map
                  is not contained in the ring of integers.  If same_ring
                  is True, then this function will return a value in self.parent()
                  or raise an error if the answer would have negative valuation.
                  If same_ring is False, then this function will raise an error
                  (this behavior is for consistency with other p-adic types).
                  If same_ring is a p-adic field into which this fixed mod ring
                  can be successfully cast, then self is cast into that field
                  and the log is taken there.  Note that this casting will
                  assume that self has the full precision possible.

        OUTPUT:
        The p-adic log of self.

        Let K be the parent of self, pi be a uniformizer of K and w be
        a generator for the group of roots of unity in K.  The usual
        power series for log with values in the additive
        group of K only converges for 1-units (units congruent to
        1 modulo pi).  However, there is a unique extension of log to a
        homomorphism defined on all the units.  If u = a*v is a unit
        with v = 1 (mod p), then we define log(u) = log(v).  This is
        the correct extension because the units U of K split as a
        product U = V x <w>, where V is the subgroup of 1-units.
        The <w> factor is torsion, so must go to 0 under any
        homomorphism to the torsion free group $(K, +)$.

        Notes -- What some other systems do with regard to non-1-units:
           PARI:  Seems to define log the same way as we do.
           MAGMA: Gives an error when unit is not a 1-unit.

        In addition, if branch is specified, then the log map
        will work on non-units:

           log(pi^k * u) = k * branch + log(u)

        Algorithm:
           Input: Some unit u.
           1. Check that the input is really a unit
              (i.e., valuation 0), or that branch is specified.
           2. Let $1-x = u^{q-1}$, which is a 1-unit, where q is the
              order of the residue field of K.
           3. Use the series expansion
              $$
                \log(1-x) = F(x) = -x - 1/2*x^2 - 1/3*x^3 - 1/4*x^4 - 1/5*x^5 - ...
              $$
              to compute the logarithm log(u**(q-1)).
              Add on terms until x^k is zero modulo the precision cap, and then
              determine if there are further terms that contribute to the sum
              (those where k is slightly above the precision cap but divisible by p).
           4. Then $$\log(u) = log(u^{q-1})/(q-1) = F(1-u^{q-1})/(q-1).$$

        EXAMPLES:
        First, the Eisenstein case.
        sage: R = ZpFM(5,5)
        sage: S.<x> = R[]
        sage: f = x^4 + 15*x^2 + 625*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = 1 + w^2 + 4*w^7; z
        1 + w^2 + 4*w^7 + O(w^20)
        sage: z.log()
        4*w^2 + 3*w^4 + w^6 + w^7 + w^8 + 4*w^9 + 3*w^10 + w^12 + w^13 + 3*w^14 + w^15 + 4*w^16 + 4*w^17 + 3*w^18 + 3*w^19 + O(w^20)

        Check that log is multiplicative:
        sage: y = 1 + 3*w^4 + w^5
        sage: y.log() + z.log() - (y*z).log()
        O(w^20)

        Now an unramified example.
        sage: g = x^3 + 3*x + 3
        sage: A.<a> = R.ext(g)
        sage: b = 1 + 5*(1 + a^2) + 5^3*(3 + 2*a)
        sage: b.log()
        (4*a^2 + 4)*5 + (a^2 + a + 2)*5^2 + (a^2 + 2*a + 4)*5^3 + (a^2 + 2*a + 2)*5^4 + O(5^5)

        Check that log is multiplicative:
        sage: c = 3 + 5^2*(2 + 4*a)
        sage: b.log() + c.log() - (b*c).log()
        O(5^5)

        AUTHORS:
            -- David Roe: initial version

        TODO:
            -- Currently implemented as $O(N^2)$. This can be improved to
            soft-$O(N)$ using algorithm described by Dan Bernstein:
            http://cr.yp.to/lineartime/multapps-20041007.pdf
        """
        if same_ring is False:
            if self.prime_pow.in_field == 0:
                return self.parent().fraction_field()(self).log(branch=branch, same_ring=True)
        elif same_ring is True:
            if self.prime_pow.in_field == 0 and mpz_cmp_ui(self.prime_pow.prime.value, self.prime_pow.e) < 0:
                raise ValueError, "result is not integral.  Use the option same_ring = False"
        else:
            return same_ring(self).log(branch = branch, same_ring=True)
        if self._is_exact_zero():
            raise ValueError, "log of zero is not defined"
        elif self._is_inexact_zero():
            raise PrecisionError, "need more precision to compute log"
        cdef ZZ_c p, q, j, ZZ_top, leftover, gap, ppow, tester, ZZ_tmp, ZZ_tmp2
        cdef ZZ_pX_c res
        cdef mpz_t tmp_m
        cdef long val = self.valuation_c()
        cdef long mini
        cdef long top
        cdef long extra_prec
        cdef long p_long
        cdef long to_shift, p_shift
        cdef long ans_ordp, ans_aprec, ans_rprec
        cdef pAdicZZpXCRElement ans
        cdef pAdicZZpXCRElement branch_add, big_oh
        cdef ZZ_pX_c y, x, xpow, to_add, to_mul
        cdef Integer Integer_val, Integer_e
        cdef ntl_ZZ to_list
        cdef bint is_one, branched, high_ramification
        cdef RealDoubleElement RDF_e, RDF_val, RDF_p, RDF_rprec, pow_switch
        cdef Integer log_e_val_floor, threshold, log_e_rprec_floor
        cdef Integer vpn

        cdef ntl_ZZ printer_ZZ
        cdef ntl_ZZ_pX printer_ZZ_pX

        self._normalize()
        self.prime_pow.restore_context_capdiv(self.relprec)
        if val != 0:
            if branch is None:
                raise ValueError, "not a unit: specify a branch of the log map"
            branched = True
            branch_add = self._new_c(self.prime_pow.ram_prec_cap)
            mpz_init(tmp_m)
            mpz_set_si(tmp_m, val)
            branch_add._set_from_mpz_rel(tmp_m, self.prime_pow.ram_prec_cap)
            mpz_clear(tmp_m)
            branch_add = <pAdicZZpXCRElement>(self.parent()(branch) * branch_add)
        else:
            branched = False
        # We may have to increase the relative precision above self.relprec.
        # So we use y for the moment, and change to x once we've determined the
        # relative precision of the answer.
        y = self.unit
        p = self.prime_pow.pow_ZZ_tmp(1)[0]
        q = self.prime_pow.pow_ZZ_tmp(self.prime_pow.f)[0]
        Integer_e = PY_NEW(Integer)
        mpz_set_si(Integer_e.value, self.prime_pow.e)
        if self.prime_pow.e == 1:
            # This function was written before residue, so we check
            # the residue the hard way.
            self.prime_pow.restore_context(1)
            ZZ_pX_conv_modulus(res, self.unit, self.prime_pow.get_context(1).x)
            if ZZ_pX_IsOne(res):
                is_one = True
                # It's already a 1-unit, so just use the series
                # (base case of "induction")
                # Set x = (1 - self).unit_part(), val = (1-self).valuation()
                self.prime_pow.restore_context(self.relprec)
            else:
                is_one = False
                # We raise to the q-1 power.  Note that this means our running time has a linear dependence on
                # the residue extension degree.
                # Set x = (1 - self^(q-1)).unit_part(), val = (1 - self^(q-1)).valuation()
                self.prime_pow.restore_context(self.relprec)
                ZZ_add_long(j, q, -1)
                ZZ_pX_PowerMod_pre(y, y, j, self.prime_pow.get_modulus(self.relprec)[0])
            ZZ_pX_sub_long(y, 1, y)
            ZZ_pX_min_val_coeff(val, mini, y, p)
            if mini == -1 or val >= self.relprec:
                #self == 1
                big_oh = self._new_c(0)
                big_oh._set_inexact_zero(self.relprec)
                if branched:
                    return branch_add + big_oh
                else:
                    return big_oh
            high_ramification = False
            ans_ordp = val
            Integer_val = PY_NEW(Integer)
            mpz_set_si(Integer_val.value, val)
            self.prime_pow.restore_context(self.relprec - val)
            ZZ_pX_right_pshift(x, y, self.prime_pow.pow_ZZ_tmp(val)[0], self.prime_pow.get_context(self.relprec - val).x)
        else:
            ZZ_rem(gap, ZZ_p_rep(ZZ_pX_ConstTerm(x)), self.prime_pow.pow_ZZ_tmp(1)[0])
            self.prime_pow.restore_context_capdiv(self.relprec)
            if ZZ_IsOne(gap):
                is_one = True
                # It's already a 1-unit, so just use the series
                # (base case of "induction")
                # Set x = (1 - self).unit_part(), val = (1-self).valuation()
            else:
                is_one = False
                # We raise to the p-1 power.
                # Set x = (1 - self^(p-1)).unit_part(), val = (1 - self^(p-1)).valuation()
                ZZ_add_long(j, p, -1)
                ZZ_pX_PowerMod_pre(y, y, j, self.prime_pow.get_modulus_capdiv(self.relprec)[0])
            ZZ_pX_sub_long(y, 1, y)
            ZZ_pX_min_val_coeff(val, mini, y, p)
            if mini == -1 or val * self.prime_pow.e + mini >= self.relprec:
                #self == 1
                big_oh = self._new_c(0)
                big_oh._set_inexact_zero(self.relprec)
                if branched:
                    return branch_add + big_oh
                else:
                    return big_oh
            val = val * self.prime_pow.e + mini
            Integer_val = PY_NEW(Integer)
            mpz_set_si(Integer_val.value, val)
            # We're looking for the minimum of the valuation of x^k/k.  Such minimums can only occur
            # at powers of p, so consider x^(p^n)/p^n, which has valuation val * p^n - e * n. Differentiating
            # and solving for n gives n = (log(e) - log(val) - log(log(p))) / log(p).

            #print "a"
            RDF_e = RealDoubleElement(float(self.prime_pow.e))
            RDF_val = RealDoubleElement(float(val))
            RDF_p = RealDoubleElement(float(self.prime_pow.prime))
            log_e_val_floor = (RDF_e / (RDF_val * RDF_p.log())).log(RDF_p).floor()
            if log_e_val_floor < 0:
                high_ramification = False
                ans_ordp = val
            else:
                # The high ramification code isn't working yet.  So we raise a NotImplementedError.
                raise NotImplementedError
                vpn = Integer_val * self.prime_pow.pow_Integer_Integer(log_e_val_floor)
                # (val * p^n - e * n) - (val * p^(n+1) - e * (n+1)) = vpn * (1 - p) + e
                if vpn * (1 - self.prime_pow.prime) + Integer_e <= 0: # initial lower
                    high_ramification = (log_e_val_floor > 0)
                    vpn = vpn - Integer_e * log_e_val_floor
                else: # final lower
                    high_ramification = True
                    vpn = vpn * self.prime_pow.prime - Integer_e * (log_e_val_floor + Integer(1))
                if mpz_cmp_si(vpn.value, -maxordp) <= 0: # I don't think this should ever happen, but...
                    raise ValueError, "valuation overflow"
                ans_ordp = mpz_get_si(vpn.value)
            # We now figure out the precision of the answer
            if high_ramification:
                #print "b"
                RDF_rprec = RealDoubleElement(float(self.relprec))
                log_e_rprec_floor = (RDF_e / (RDF_rprec * RDF_p.log())).log(RDF_p).floor()
                if log_e_rprec_floor < 0:
                    log_e_rprec_floor = Integer(0)
                if mpz_cmp_ui(self.prime_pow.prime.value, maxordp) >= 0:
                    # e / (p - 1) < 1 so we never have to multiply by p.
                    # Thus the minimum absolute precision of x^k/k is just self.relprec.
                    ans_aprec = self.relprec
                else:
                    #print "c"
                    p_long = mpz_get_si(self.prime_pow.prime.value)
                    if self.prime_pow.e < p_long - 1:
                        ans_aprec = self.relprec
                    else:
                        # We're looking for the minimum of the absolute precision of x^(p^k)/p^k.
                        # The precision behavior of exponentiation changes when the relative precision
                        # hits e / (p-1).
                        # Below the threshold, our minimum is the minimum of (self.relprec * p^k - e*k)
                        # Above, this is always increasing (subtracting e due to the extra p in the
                        # denominator cancels with adding e because of the extra p in the exponent)
                        threshold = (RealDoubleElement(float(self.prime_pow.e)) / (RealDoubleElement(float(p_long - 1) * float(self.relprec - val)))).log(RDF_p).ceiling()
                        if log_e_rprec_floor < threshold:
                            #print "d"
                            #print log_e_rprec_floor
                            vpn = self.precision_relative() * self.prime_pow.pow_Integer_Integer(log_e_rprec_floor)
                            if vpn * (1 - self.prime_pow.prime) + Integer_e <= 0: # initial lower
                                vpn = vpn - Integer_e * log_e_rprec_floor
                            else: # final lower
                                vpn = vpn * self.prime_pow.prime - Integer_e * (log_e_rprec_floor + Integer(1))
                            ans_aprec = mpz_get_si(vpn.value)
                        else:
                            #print "e"
                            # if e is divisible by p-1 we may have a floating point error.  So we check.
                            if threshold > 0 and self.prime_pow.e % (p_long - 1) == 0:
                                vpn = self.prime_pow.pow_Integer_Integer(threshold - Integer(1))
                                mpz_mul_si(vpn.value, vpn.value, self.relprec - val)
                                mpz_mul_si(vpn.value, vpn.value, p_long - 1)
                                if mpz_cmp_si(vpn.value, self.prime_pow.e) > 0:
                                    threshold = threshold - Integer(1)
                            # After threshold, our function strictly increases.  So we need to check the values at threshold
                            # and threshold - 1
                            if threshold > 0:
                                vpn = self.precision_relative() * self.prime_pow.pow_Integer_Integer(threshold - Integer(1))
                                if vpn * (1 - self.prime_pow.prime + Integer_e) <= 0:
                                    #print "f"
                                    vpn = vpn - Integer_e * (threshold - Integer(1))
                                else:
                                    #print "g"
                                    vpn = vpn * self.prime_pow.prime - Integer_e * threshold
                            else:
                                #print "h"
                                vpn = self.precision_relative()
                            ans_aprec = mpz_get_si(vpn.value)
            else:
                ans_aprec = self.relprec
            ans_rprec = ans_aprec - ans_ordp
            if ans_rprec > self.prime_pow.ram_prec_cap:
                ans_rprec = self.prime_pow.ram_prec_cap
                ans_aprec = ans_ordp + ans_rprec
            self.prime_pow.restore_context_capdiv(ans_rprec)
            self.prime_pow.eis_shift_capdiv(&x, &y, val, ans_rprec)

        if ans_ordp < 0 and self.prime_pow.in_field == 0:
            raise ValueError, "Answer has negative valuation"
        ans = self._new_c(-ans_rprec)
        ans.ordp = ans_ordp

        # Need extra precision to take into account powers of p
        # in the denominators of the series. (Indeed, it's a
        # not-entirely-trivial fact that if x is given mod p^n, that
        # log(x) is well-defined mod p^n !) Specifically:
        # we are only guaranteed that $x^j/j$ is zero mod $p^n$ if
        # j >= floor(log_p(j)) + n.
        # But we only actually need to do this extra computation
        # if there is some j with j - j.valuation(p) < n.

        #print "ans_ordp = %s"%ans_ordp
        #print "ans_aprec = %s"%ans_aprec
        #print "ans_rprec = %s"%ans_rprec

        top = (ans_aprec - 1) / val + 1
        if top < 1:
            top = 1
        L = []
        cdef long check_exp = 1
        cdef ZZ_c check, ZZ_aprec, one
        cdef long k
        #printer_ZZ = PY_NEW(ntl_ZZ)
        if mpz_cmp_si(self.prime_pow.prime.value, maxordp) < 0:
            ZZ_conv_from_long(ZZ_top, top)
            ZZ_add_long(ZZ_top, ZZ_top, -1) # we want divisions by ppow to round up, so we subtract 1 here and add after dividing.
            ZZ_conv_from_long(ZZ_aprec, ans_aprec)
            ZZ_conv_from_long(one, 1)
            to_list = PY_NEW(ntl_ZZ)
            ppow = p
            while True:
                ZZ_div(tester, ZZ_top, ppow)
                #printer_ZZ.x = tester
                #print "pre tester = %s"%printer
                ZZ_add_long(tester, tester, 1)
                if ZZ_divide_test(tester, p):
                    # we want tester to skip multiples of p.
                    ZZ_add_long(tester, tester, 1)
                while True:
                    if ZZ_divide_test(tester, p) == 0:
                        ZZ_mul(to_list.x, tester, ppow)
                        ZZ_mul_long(check, to_list.x, val)
                        ZZ_add_long(check, check, -self.prime_pow.e * check_exp)
                        if ZZ_compare(check, ZZ_aprec) < 0:
                            #print "appending %s"%(to_list)
                            L.append(to_list)
                            to_list = PY_NEW(ntl_ZZ)
                        else:
                            break
                    ZZ_add_long(tester, tester, 1)
                #to_list.x = ppow
                #print "ppow = %s"%to_list
                #to_list.x = tester
                #print "tester = %s"%to_list
                if ZZ_compare(tester, one) == 0 and mpz_cmp_si(log_e_val_floor.value, check_exp) < 0:
                    break
                check_exp += 1
                ZZ_mul(ppow, ppow, p)

        xpow = x
        ZZ_conv_from_long(j, 1)

        #printer_ZZ = PY_NEW(ntl_ZZ)
        #printer_ZZ_pX = ntl_ZZ_pX([], self.prime_pow.get_context_capdiv(ans_rprec))
        while ZZ_compare(j, ZZ_top) <= 0:
            ZZ_conv_to_long(to_shift, j)
            to_shift *= val
            if ZZ_divide_test(j, p):
                p_shift = ZZ_remove(leftover, j, p)
            else:
                p_shift = 0
                leftover = j
            #print "a"
            #printer_ZZ.x = leftover
            #print "leftover = %s"%printer_ZZ
            #printer_ZZ.x = p
            #print "p = %s"%printer_ZZ
            #printer_ZZ.x = self.prime_pow.pow_ZZ_top()[0]
            #print "p^n = %s"%printer_ZZ
            #print "p_shift = %s"%p_shift
            #print "a"
            _sig_on
            ZZ_InvMod(leftover, leftover, self.prime_pow.pow_ZZ_tmp(ans_rprec)[0])
            _sig_off
            #print "b"
            ZZ_pX_mul_ZZ_p(to_add, xpow, ZZ_to_ZZ_p(leftover))
            if self.prime_pow.e == 1:
                ZZ_pX_left_pshift(to_add, to_add, self.prime_pow.pow_ZZ_tmp(to_shift - p_shift - ans_ordp)[0], self.prime_pow.get_context(ans_rprec).x)
            else:
                ##printer_ZZ.x = self.prime_pow.pow_ZZ_tmp(p_shift)[0]
                ##print printer_ZZ
                #ZZ_pX_clear(to_mul)
                #ZZ_pX_SetCoeff(to_mul, 0, ZZ_to_ZZ_p(self.prime_pow.pow_ZZ_tmp(p_shift)[0]))
                #self.prime_pow.eis_shift(&to_mul, &to_mul, p_shift * self.prime_pow.e, self.prime_pow.ram_prec_cap)
                ##printer_ZZ_pX.x = to_mul
                ##print printer_ZZ_pX
                ##print "c"
                if PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_small_Eis):
                    ZZ_pX_conv_modulus(to_mul, (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_shifter[0], self.prime_pow.get_context_capdiv(ans_rprec).x)
                elif PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_big_Eis):
                    ZZ_pX_conv_modulus(to_mul, (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_shifter[0], self.prime_pow.get_context_capdiv(ans_rprec).x)
                else:
                    raise RuntimeError, "unrecognized PowComputer type"
                _sig_on
                ZZ_pX_InvMod_newton_ram(to_mul, to_mul, self.prime_pow.get_modulus_capdiv(ans_rprec)[0], self.prime_pow.get_context_capdiv(ans_rprec).x)
                _sig_off
                #print "d"
                ZZ_pX_PowerMod_long_pre(to_mul, to_mul, p_shift, self.prime_pow.get_modulus_capdiv(ans_rprec)[0])
                ZZ_pX_MulMod_pre(to_add, to_add, to_mul, self.prime_pow.get_modulus_capdiv(ans_rprec)[0])
                self.prime_pow.eis_shift_capdiv(&to_add, &to_add, -(to_shift - p_shift * self.prime_pow.e - ans_ordp), ans_rprec)
            ##printer_ZZ.x = j
            ##print "j = %s"%printer_ZZ
            #printer_ZZ_pX.x = to_add
            #print "to_add = %s"%printer_ZZ_pX
            #printer_ZZ_pX.x = xpow
            #print "xpow = %s\n"%printer_ZZ_pX
            ZZ_pX_add(ans.unit, ans.unit, to_add)
            ZZ_pX_MulMod_pre(xpow, xpow, x, self.prime_pow.get_modulus_capdiv(ans_rprec)[0])
            ZZ_add_long(j, j, 1)
        #print "starting L"
        L.sort()
        for m in L:
            p_shift = ZZ_remove(leftover, (<ntl_ZZ>m).x, p)
            ZZ_mul_long(ZZ_tmp2, (<ntl_ZZ>m).x, val)
            ZZ_conv_from_long(ZZ_tmp, self.prime_pow.e)
            ZZ_mul_long(ZZ_tmp, ZZ_tmp, p_shift)
            ZZ_sub(ZZ_tmp, ZZ_tmp2, ZZ_tmp)
            ZZ_conv_to_long(to_shift, ZZ_tmp)
            to_shift = to_shift - ans_ordp
            #print "e"
            _sig_on
            ZZ_InvMod(leftover, leftover, self.prime_pow.pow_ZZ_tmp(ans_rprec)[0])
            _sig_off
            #print "f"
            ZZ_pX_mul_ZZ_p(to_add, xpow, ZZ_to_ZZ_p(leftover))
            if self.prime_pow.e == 1:
                ZZ_pX_left_pshift(to_add, to_add, self.prime_pow.pow_ZZ_tmp(to_shift)[0], self.prime_pow.get_context(ans_rprec).x)
            else:
                #ZZ_pX_clear(to_mul)
                #ZZ_pX_SetCoeff(to_mul, 0, ZZ_to_ZZ_p(self.prime_pow.pow_ZZ_tmp(p_shift)[0]))
                #self.prime_pow.eis_shift(&to_mul, &to_mul, p_shift * self.prime_pow.e, self.prime_pow.ram_prec_cap)
                #print "g"
                if PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_small_Eis):
                    ZZ_pX_conv_modulus(to_mul, (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_shifter[0], self.prime_pow.get_context_capdiv(ans_rprec).x)
                elif PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_big_Eis):
                    ZZ_pX_conv_modulus(to_mul, (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_shifter[0], self.prime_pow.get_context_capdiv(ans_rprec).x)
                else:
                    raise RuntimeError, "unrecognized PowComputer type"
                _sig_on
                ZZ_pX_InvMod_newton_ram(to_mul, to_mul, self.prime_pow.get_modulus_capdiv(ans_rprec)[0], self.prime_pow.get_context_capdiv(ans_rprec).x)
                _sig_off
                #print "h"
                ZZ_pX_PowerMod_long_pre(to_mul, to_mul, p_shift, self.prime_pow.get_modulus_capdiv(ans_rprec)[0])
                ZZ_pX_MulMod_pre(to_add, to_add, to_mul, self.prime_pow.get_modulus_capdiv(ans_rprec)[0])
                self.prime_pow.eis_shift_capdiv(&to_add, &to_add, -to_shift, ans_rprec)
            #print "m = %s"%m
            #printer_ZZ_pX.x = to_add
            #print "to_add = %s"%printer_ZZ_pX
            ZZ_pX_add(ans.unit, ans.unit, to_add)
            ZZ_sub(gap, (<ntl_ZZ>m).x, ZZ_top)
            ZZ_top = (<ntl_ZZ>m).x
            #printer_ZZ.x = gap
            #print "gap = %s"%printer_ZZ
            ZZ_pX_PowerMod_pre(to_mul, x, gap, self.prime_pow.get_modulus_capdiv(ans_rprec)[0])
            ZZ_pX_MulMod_pre(xpow, xpow, to_mul, self.prime_pow.get_modulus_capdiv(ans_rprec)[0])
            #print "loopdone\n"
        if not is_one:
            ZZ_add_long(q, q, -1)
            ZZ_rem(q, q, self.prime_pow.pow_ZZ_tmp(ans_rprec)[0])
            #print "i"
            _sig_on
            ZZ_InvMod(q, q, self.prime_pow.pow_ZZ_tmp(ans_rprec)[0])
            _sig_off
            #print "j"
            ZZ_pX_mul_ZZ_p(ans.unit, ans.unit, ZZ_to_ZZ_p(q))
        if branched:
            return ans + branch_add
        else:
            return ans

    def multiplicative_order(self, prec=None):
        """
        Returns the multiplicative order of self, ie the smallest positive n so that there is an exact p-adic element
        congruent to self modulo self's precision that is an nth root of unity.

        Note: unlike the case for Qp and Zp, it is possible to have non-teichmuller elements with finite orders.
        This can happen only if (p-1) divides the ramification index (see the documentation on __pow__).

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            integer -- the multiplicative order of self
        """
        raise NotImplementedError

    def teichmuller_list(self):
        raise NotImplementedError

    def _teichmuller_set(self):
        """
        Sets self to the teichmuller representative congruent to self modulo pi, with
        the same relative precision as self.

        This function should not be used externally: elements are supposed to be immutable.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: y = W.teichmuller(3); y #indirect doctest
        3 + 3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + 3*w^15 + 2*w^16 + 3*w^17 + w^18 + 3*w^19 + 3*w^20 + 2*w^21 + 2*w^22 + 3*w^23 + 4*w^24 + O(w^25)
        sage: y^5 == y
        True
        sage: g = x^3 + 3*x + 3
        sage: A.<a> = R.ext(g)
        sage: b = A.teichmuller(1 + 2*a - a^2); b
        (4*a^2 + 2*a + 1) + 2*a*5 + (3*a^2 + 1)*5^2 + (a + 4)*5^3 + (a^2 + a + 1)*5^4 + O(5^5)
        sage: b^125 == b
        True
        """
        self._normalize()
        if self.ordp > 0:
            self._set_exact_zero()
        elif self.relprec == 0:
            raise ValueError, "not enough precision known"
        else:
            self.prime_pow.teichmuller_set_c(&self.unit, &self.unit, self.relprec)

    def padded_list(self, n, lift_mode = 'simple'):
        raise NotImplementedError

    def precision_absolute(self):
        """
        Returns the absolute precision of self, ie the power of the uniformizer modulo which
        this element is defined.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(75, 19); a
        3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
        sage: a.valuation()
        10
        sage: a.precision_absolute()
        19
        sage: a.precision_relative()
        9
        sage: a.unit_part()
        3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        cdef Integer ans
        if self.ordp == maxordp:
            return infinity
        else:
            ans = PY_NEW(Integer)
            mpz_set_si(ans.value, self.relprec + self.ordp)
            return ans

    def precision_relative(self):
        """
        Returns the relative precision of self, ie the power of the uniformizer modulo which
        the unit part of self is defined.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(75, 19); a
        3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
        sage: a.valuation()
        10
        sage: a.precision_absolute()
        19
        sage: a.precision_relative()
        9
        sage: a.unit_part()
        3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        self._normalize()
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.relprec)
        return ans

    def residue(self, n):

        raise NotImplementedError

    cdef long valuation_c(self):
        """
        Returns the valuation of self.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(75, 19); a
        3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
        sage: a.valuation() # indirect doctest
        10
        sage: a.precision_absolute()
        19
        sage: a.precision_relative()
        9
        sage: a.unit_part()
        3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        self._normalize()
        return self.ordp

    cpdef pAdicZZpXCRElement unit_part(self):
        """
        Returns the unit part of self, ie self / uniformizer^(self.valuation())

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = W(75, 19); a
        3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
        sage: a.valuation()
        10
        sage: a.precision_absolute()
        19
        sage: a.precision_relative()
        9
        sage: a.unit_part()
        3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        self._normalize()
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        ans.ordp = 0
        ans.unit = self.unit
        return ans

    cdef ext_p_list(self, bint pos):
        """
        Returns a list of integers (in the eisenstein case) or a list of lists of integers (in the unramified case).
        self can be reconstructed as a sum of elements of the list times powers of the uniformiser (in the eisenstein case),
        or as a sum of powers of the p times polynomials in the generator (in the unramified case).
        If pos is True, all integers will be in the range [0,p-1], otherwise they will be in the range [(1-p)/2, p/2].

        Note that zeros are truncated from the returned list, so you must use the valuation() function to completely recover self.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: y = W(775, 19); y
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
        sage: y._ext_p_list(True)
        [1, 0, 4, 0, 2, 1, 2, 4, 1]
        sage: y._ext_p_list(False)
        [1, 0, -1, 0, 2, 1, 2, 0, 1]
        sage: w^10 - w^12 + 2*w^14 + w^15 + 2*w^16 + w^18 + O(w^19)
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
        sage: g = x^3 + 3*x + 3
        sage: A.<a> = R.ext(g)
        sage: y = 75 + 45*a + 1200*a^2; y
        4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
        sage: y._ext_p_list(True)
        [[0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1]]
        sage: y._ext_p_list(False)
        [[0, -1], [-2, 2, -2], [1], [0, 0, 2]]
        sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
        4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
        """
        self._normalize()
        return self.ext_p_list_precs(pos, self.relprec)

def make_ZZpXCRElement(parent, unit, ordp, relprec, version):
    cdef pAdicZZpXCRElement ans
    cdef ZZ_pX_c poly
    if version == 0:
        ans = pAdicZZpXCRElement(parent, [], empty = True)
        ans.prime_pow.restore_context_capdiv(mpz_get_si((<Integer>relprec).value))
        poly = (<ntl_ZZ_pX>unit).x
        ans._set(&poly, mpz_get_si((<Integer>ordp).value), mpz_get_si((<Integer>relprec).value))
        return ans
    else:
        raise ValueError, "unknown unpickling version"
