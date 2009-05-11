include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
include "../../ext/python_list.pxi"
include "../../libs/ntl/decl.pxi"

from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement
from sage.rings.integer_mod import is_IntegerMod
from sage.rings.padics.padic_printing cimport pAdicPrinter_class

zero = Integer(0)
one = Integer(1)
two = Integer(2)
big = two**128 + one
#this should not fit in a long, since it's supposed to be bigger than any valid absolute precision.

cdef class pAdicZZpXElement(pAdicExtElement):
    def __init__(self, parent):
        self.prime_pow = <PowComputer_ZZ_pX>parent.prime_pow
        pAdicExtElement.__init__(self, parent)

    cdef int _set_from_list(self, L) except -1:
        """
        Sets self from a list.

        The list should either be uniform in type, or all of the entries should be coercible to integers.
        If any of the entries in L is a list, L will be cast to a ZZ_pEX

        INPUT:
        L -- a list.
        """
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX((<ntl_ZZX>ntl_ZZX(L)).x)
        else:
            self._set_from_ZZ_pX(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx)
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef int _set_from_list_rel(self, L, long relprec) except -1:
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX_rel((<ntl_ZZX>ntl_ZZX(L)).x, relprec)
        else:
            self._set_from_ZZ_pX_rel(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx, relprec)
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef int _set_from_list_abs(self, L, long absprec) except -1:
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX_abs((<ntl_ZZX>ntl_ZZX(L)).x, absprec)
        else:
            self._set_from_ZZ_pX_abs(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx, absprec)
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef int _set_from_list_both(self, L, long absprec, long relprec) except -1:
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX_both((<ntl_ZZX>ntl_ZZX(L)).x, absprec, relprec)
        else:
            self._set_from_ZZ_pX_both(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx, absprec, relprec)
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef long _check_ZZ_pContext(self, ntl_ZZ_pContext_class ctx) except -1:
        cdef ZZ_c leftover
        cdef long val = ZZ_remove(leftover, ctx.p.x, self.prime_pow.pow_ZZ_tmp(1)[0])
        if ZZ_IsOne(leftover):
            return val
        else:
            raise ValueError, "context must be a power of the appropriate prime"

    cdef ext_p_list_precs(self, bint pos, long prec):
        ans = []
        cdef ntl_ZZ ZZ_coeff = ntl_ZZ()
        cdef Integer coeff = PY_NEW(Integer)
        cdef Integer zero = Integer(0)
        cdef Integer list_elt
        cdef ZZ_c halfp
        cdef Py_ssize_t i, j
        cdef ZZ_p_c const_term_holder
        self.prime_pow.restore_top_context()
        ###ZZ_p_construct(&const_term_holder)
        cdef ntl_ZZ holder = ntl_ZZ()
        cdef ZZ_p_c tmp
        cdef pAdicPrinter_class printer = <pAdicPrinter_class>self.parent()._printer
        cdef ZZ_pX_c shifter = (<ntl_ZZ_pX>self._ntl_rep()).x

        #cdef ntl_ZZ_pContext_class cup = self.prime_pow.get_context(self.prime_pow.prec_cap + (<PowComputer_ZZ_pX_FM_Eis>self.prime_pow).low_length)
        #cdef ntl_ZZ_pX printer = ntl_ZZ_pX([],cup)
        #printer.x = ((<PowComputer_ZZ_pX_FM_Eis>self.prime_pow).low_shifter[0]).val()
        #print printer

        if self.prime_pow.e == 1:
            for j from 0 <= j < self.prime_pow.prec_cap:
                ans.append([])
            for i from 0 <= i < self.prime_pow.deg:
                ZZ_coeff.x = ZZ_p_rep(ZZ_pX_coeff(shifter, i))
                ZZ_to_mpz(&coeff.value, &ZZ_coeff.x)
                L = printer.base_p_list(coeff.value, pos)
                for j from 0 <= j < prec:
                    if j < len(L):
                        ans[j].append(L[j])
                    else:
                        ans[j].append(zero)
            for j from 0 <= j < prec:
                while len(ans[j]) > 0:
                    if ans[j][-1] == 0:
                        ans[j].pop()
                    else:
                        break
            zerotest = []
        else:
            halfp = self.prime_pow.pow_ZZ_tmp(1)[0]
            ZZ_DivRem_long(halfp, halfp, 2)
            i = 0
            while True:
                #print shifter._ntl_rep()
                # It's important that one doesn't normalize in between shifting (for capped relative elements):
                # _const_term doesn't normalize and thus we pick up the zeros
                # since we're throwing away leading zeros, it doesn't matter if we start normalized or not.
                for j from 0 <= j < self.prime_pow.e:
                    list_elt = PY_NEW(Integer)
                    if i + j == prec:
                        break
                    ZZ_rem(holder.x, ZZ_p_rep(ZZ_pX_coeff(shifter, j)), self.prime_pow.pow_ZZ_tmp(1)[0])
                    if not pos and not ZZ_IsZero(holder.x) and ZZ_compare(holder.x, halfp) > 0:
                        ZZ_sub(holder.x, self.prime_pow.pow_ZZ_tmp(1)[0], holder.x)
                        ZZ_p_add(tmp, ZZ_to_ZZ_p(holder.x), ZZ_pX_coeff(shifter, j))
                        ZZ_pX_SetCoeff(shifter, j, tmp)
                        ZZ_negate(holder.x, holder.x)
                    ZZ_to_mpz(&list_elt.value, &holder.x)
                    ans.append(list_elt)
                i += self.prime_pow.e
                if i >= prec:
                    break
                self.prime_pow.eis_shift(&shifter, &shifter, self.prime_pow.e, self.prime_pow.capdiv(prec - i))
            zerotest = 0
        while len(ans) > 0:
            if ans[-1] == zerotest:
                ans.pop()
            else:
                break
        while len(ans) > 0:
            if ans[0] == zerotest:
                ans.pop(0)
            else:
                break
        return ans

    def norm(self, base = None):
        """
        Return the absolute or relative norm of this element.

        NOTE!  This is not the p-adic absolute value.  This is a field theoretic norm down to a ground ring.
        If you want the p-adic absolute value, use the abs() function instead.

        If K is given then K must be a subfield of the parent L of
        self, in which case the norm is the relative norm from L to K.
        In all other cases, the norm is the absolute norm down to Qp or Zp.

        EXAMPLES:
        sage: R = ZpCR(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: ((1+2*w)^5).norm()
        1 + 5^2 + O(5^5)
        sage: ((1+2*w)).norm()^5
        1 + 5^2 + O(5^5)

        TESTS:
        sage: R = ZpCA(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: ((1+2*w)^5).norm()
        1 + 5^2 + O(5^5)
        sage: ((1+2*w)).norm()^5
        1 + 5^2 + O(5^5)
        sage: R = ZpFM(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: ((1+2*w)^5).norm()
        1 + 5^2 + O(5^5)
        sage: ((1+2*w)).norm()^5
        1 + 5^2 + O(5^5)
        """
        if base is not None:
            if base is self.parent():
                return self
            else:
                raise NotImplementedError
        if self._is_exact_zero():
            return self.parent().ground_ring()(0)
        elif self._is_inexact_zero():
            return self.ground_ring(0, self.valuation())
        norm_of_uniformizer = (-1)**self.parent().degree() * self.parent().defining_polynomial()[0]
        if self.valuation() == 0:
            return self.parent().ground_ring()(self.unit_part().matrix_mod_pn().det())
        else:
            return self.parent().ground_ring()(self.unit_part().matrix_mod_pn().det()) * norm_of_uniformizer**self.valuation()

    def trace(self, base = None):
        """
        Return the absolute or relative trace of this element.

        If K is given then K must be a subfield of the parent L of
        self, in which case the norm is the relative norm from L to K.
        In all other cases, the norm is the absolute norm down to Qp or Zp.

        EXAMPLES:
        sage: R = ZpCR(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = (2+3*w)^7
        sage: b = (6+w^3)^5
        sage: a.trace()
        3*5 + 2*5^2 + 3*5^3 + 2*5^4 + O(5^5)
        sage: a.trace() + b.trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        sage: (a+b).trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)

        TESTS:
        sage: R = ZpCA(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = (2+3*w)^7
        sage: b = (6+w^3)^5
        sage: a.trace()
        3*5 + 2*5^2 + 3*5^3 + 2*5^4 + O(5^5)
        sage: a.trace() + b.trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        sage: (a+b).trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        sage: R = ZpFM(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = (2+3*w)^7
        sage: b = (6+w^3)^5
        sage: a.trace()
        3*5 + 2*5^2 + 3*5^3 + 2*5^4 + O(5^5)
        sage: a.trace() + b.trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        sage: (a+b).trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        """
        if base is not None:
            if base is self.parent():
                return self
            else:
                raise NotImplementedError
        if self._is_exact_zero():
            return self.parent().ground_ring()(0)
        elif self._is_inexact_zero():
            return self.ground_ring(0, (self.valuation() - 1) // self.parent().e() + 1)
        if self.valuation() >= 0:
            return self.parent().ground_ring()(self.matrix_mod_pn().trace())
        else:
            shift = -(self.valuation() // self.parent().e())
            return self.parent().ground_ring()((self * self.parent().prime() ** shift).matrix_mod_pn().trace()) / self.parent().prime()**shift

    def _ntl_rep(self):
        raise NotImplementedError

    def _prime_pow(self):
        return self.prime_pow

cdef preprocess_list(pAdicZZpXElement elt, L):
    cdef Py_ssize_t i
    cdef ZZ_c tmp
    cdef ntl_ZZ_pContext_class ctx
    cdef ntl_ZZ pshift_z
    cdef Integer pshift_m
    cdef long aprec
    cdef ntl_ZZ py_tmp
    if not PyList_Check(L):
        raise TypeError, "L must be a list"
    #print "before find_val_aprec"
    min_val, min_aprec, total_type = find_val_aprec(elt, L)
    #return "a","b","c"
    if total_type == two:
        # all integers
        return [ntl_ZZ(a) for a in L], zero, None
    if min_val < 0 and not elt.prime_pow.in_field:
        raise ValueError, "negative valuation"
    if total_type == one:
        # rationals and integers
        py_tmp = PY_NEW(ntl_ZZ)
        py_tmp.x = elt.prime_pow.pow_ZZ_top()[0]
        ctx = ntl_ZZ_pContext(py_tmp)
    else:
        # integers, rationals and things with finite precision
        # note that min_val will be non-positive since things with finite precision return non-positive valuation from get_val_prec
        py_tmp = PY_NEW(ntl_ZZ)
        py_tmp.x = elt.prime_pow.pow_ZZ_tmp(mpz_get_ui((<Integer>(min_aprec - min_val)).value))[0]
        ctx = ntl_ZZ_pContext(py_tmp)
    if min_val < 0:
        pshift_z = PY_NEW(ntl_ZZ)
        pshift_z.x = elt.prime_pow.pow_ZZ_tmp(-mpz_get_si((<Integer>min_val).value))[0]
        pshift_m = elt.prime_pow.pow_Integer(-mpz_get_si((<Integer>min_val).value))
        for i from 0 <= i < len(L):
            if PY_TYPE_CHECK(L[i], ntl_ZZ):
                L[i] = ntl_ZZ_p(L[i]*pshift_z, ctx)
            elif PY_TYPE_CHECK(L[i], Integer) or PY_TYPE_CHECK(L[i], Rational) or isinstance(L[i], (int, long)):
                L[i] = ntl_ZZ_p(L[i]*pshift_m, ctx)
            elif PY_TYPE_CHECK(L[i], pAdicBaseGenericElement):
                L[i] = ntl_ZZ_p((L[i] << min_val).lift(), ctx)
            elif is_IntegerMod(L[i]):
                L[i] = ntl_ZZ_p(L[i].lift()*pshift_m, ctx)
            elif (L[i].modulus_context() is not ctx) or min_val != zero:
                L[i] = ntl_ZZ_p(L[i].lift()*pshift_z, ctx)
    elif elt.parent().is_capped_relative() and min_val > 0:
        pshift_z = PY_NEW(ntl_ZZ)
        pshift_z.x = elt.prime_pow.pow_ZZ_tmp(mpz_get_ui((<Integer>min_val).value))[0]
        pshift_m = elt.prime_pow.pow_Integer(mpz_get_ui((<Integer>min_val).value))
        for i from 0 <= i < len(L):
            if PY_TYPE_CHECK(L[i], ntl_ZZ):
                ZZ_div(tmp, (<ntl_ZZ>L[i]).x, pshift_z.x)
                py_tmp = PY_NEW(ntl_ZZ)
                py_tmp.x = tmp
                L[i] = ntl_ZZ_p(py_tmp, ctx)
            elif PY_TYPE_CHECK(L[i], Integer) or PY_TYPE_CHECK(L[i], Rational) or isinstance(L[i], (int, long)):
                L[i] = ntl_ZZ_p(L[i]//pshift_m, ctx)
            elif PY_TYPE_CHECK(L[i], pAdicBaseGenericElement):
                L[i] = ntl_ZZ_p((L[i] << min_val).lift(), ctx)
            elif is_IntegerMod(L[i]):
                L[i] = ntl_ZZ_p(L[i].lift()//pshift_m, ctx)
            elif (L[i].modulus_context() is not ctx) or min_val != zero:
                ZZ_div(tmp, (<ntl_ZZ>L[i].lift()).x, pshift_z.x)
                py_tmp = PY_NEW(ntl_ZZ)
                py_tmp.x = tmp
                L[i] = ntl_ZZ_p(py_tmp, ctx)
    else:
        for i from 0 <= i < len(L):
            if PY_TYPE_CHECK(L[i], ntl_ZZ) or PY_TYPE_CHECK(L[i], Integer) or PY_TYPE_CHECK(L[i], Rational) or isinstance(L[i], (int, long)):
                L[i] = ntl_ZZ_p(L[i], ctx)
            elif PY_TYPE_CHECK(L[i], pAdicBaseGenericElement) or is_IntegerMod(L[i]) or (L[i].modulus_context() is not ctx):
                L[i] = ntl_ZZ_p(L[i].lift(), ctx)
    return L, min_val, ctx

cdef find_val_aprec(pAdicZZpXElement elt, L):
    cdef Py_ssize_t i
    min_val = big
    min_aprec = big
    total_type = two # we begin by defaulting to the list elements being integers
    for i from 0 <= i < len(L):
        #print "before get_val_prec"
        cur_val, cur_aprec, cur_type = get_val_prec(elt, L[i])
        #return "a","b","c"
        # proc_type == 0 indicates something with finite precision
        # proc_type == 1 indicates a rational, or something that cannot be coerced to an integer
        # proc_type == 2 indicates an integer, or something that can be coerced to an integer
        # but can be coerced to Z/p^n for any n.
        if cur_aprec < min_aprec:
            min_aprec = cur_aprec
        if cur_val < min_val:
            min_val = cur_val
        if cur_type < total_type:
            cur_type = total_type
    return min_val, min_aprec, total_type

cdef get_val_prec(pAdicZZpXElement elt, a):
    cdef ntl_ZZ py_tmp
    #print "pre Integer check"
    if PY_TYPE_CHECK(a, Integer):
        return (a.valuation(elt.prime_pow.prime), big, two)
    #print "pre ntl_ZZ check"
    if PY_TYPE_CHECK(a, ntl_ZZ):
        py_tmp = PY_NEW(ntl_ZZ)
        py_tmp.x = elt.prime_pow.pow_ZZ_tmp(1)[0]
        return (a.valuation(py_tmp), big, two)
    #print "pre int/long check"
    if isinstance(a, (int, long)):
        #print a, type(a)
        #print elt.prime_pow
        #print elt.prime_pow.prime
        #print big, type(big)
        #print two, type(two)
        return (Integer(a).valuation(elt.prime_pow.prime), big, two)
    #print "pre Rational check"
    if PY_TYPE_CHECK(a, Rational):
        val = a.valuation(elt.prime_pow.prime)
        return (val, big, one)
    #print "pre padic-base check"
    if PY_TYPE_CHECK(a, pAdicBaseGenericElement):
        if a.prime() == elt.prime():
            val = a.valuation()
            return (val if val < zero else zero, a.precision_absolute(), zero)
        else:
            raise TypeError, "primes must match"
    cdef mpz_t leftover
    cdef long long_val
    cdef Integer Integer_val
    #print "pre IntegerMod check"
    if is_IntegerMod(a):
        mpz_init(leftover)
        long_val = mpz_remove(leftover, (<Integer>a.modulus()).value, elt.prime_pow.prime.value)
        if long_val > 0 and mpz_cmp_ui(leftover, 1) == 0:
            mpz_clear(leftover)
            Integer_val = PY_NEW(Integer)
            mpz_set_ui(Integer_val.value, long_val)
            # Since we're guaranteed to be in type 0, we don't care about computing the actual valuation
            return (zero, Integer_val, zero)
        else:
            mpz_clear(leftover)
            raise TypeError, "modulus must be a positive power of the appropriate prime"
    cdef ZZ_c leftover_z
    #print "pre ntl_ZZ_p check"
    if PY_TYPE_CHECK(a, ntl_ZZ_p):
        long_val = ZZ_remove(leftover_z, ZZ_p_rep((<ntl_ZZ_p>a).x), (<ntl_ZZ_p>a).c.p.x)
        if long_val > 0 and ZZ_IsOne(leftover_z):
            Integer_val = PY_NEW(Integer)
            mpz_set_ui(Integer_val.value, long_val)
            # Since we're guaranteed to be in type 0, we don't care about computing the actual valuation
            return (zero, Integer_val, zero)
        else:
            raise TypeError, "modulus must be a positive power of the appropriate prime"
    raise TypeError, "unsupported type for list element: %s"%(type(a))


