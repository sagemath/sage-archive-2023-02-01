include "../../libs/ntl/decl.pxi"
include "../../ext/cdefs.pxi"
include "../../ext/stdsage.pxi"

from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.structure.element cimport Element
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX_FM_Eis
from sage.rings.integer cimport Integer

cdef class EisensteinFixedModElement(pAdicZZpXFMElement):
    def __init__(self, parent, x, absprec = None, relprec = None):
        """
        """
        # Need to do preprocessing on x
        pAdicZZpXFMElement.__init__(self, parent, x)

    def __reduce__(self):
        """
        sage: a = ZpFM(5)(-3)
        sage: type(a)
        <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
        sage: loads(dumps(a)) == a
        True
        """
        self.prime_pow.restore_top_context()
        cdef ntl_ZZ_pX holder = PY_NEW(ntl_ZZ_pX)
        holder.c = self.prime_pow.get_top_context()
        holder.x = self.value
        return make_EisensteinFixedModElement, (self.parent(), holder)

    cdef pAdicZZpXFMElement _new_c(self):
        self.prime_pow.restore_top_context()
        cdef EisensteinFixedModElement ans = PY_NEW(EisensteinFixedModElement)
        ans._parent = self._parent
        ZZ_pX_construct(&ans.value)
        ans.prime_pow = self.prime_pow
        return ans

    def __richcmp__(left, right, op):
        return (<Element>left)._richcmp(right, op)

    cdef pAdicZZpXFMElement _lshift_c(self, long n, bint top_zeros):
        if n < 0:
            return self._rshift_c(-n, top_zeros)
        elif n == 0:
            return self
        cdef pAdicZZpXFMElement ans = self._new_c()
        if n < self.prime_pow.ram_prec_cap:
            ZZ_pX_PowerXMod_long_pre(ans.value, n, self.prime_pow.get_top_modulus()[0])
            ZZ_pX_MulMod_pre(ans.value, ans.value, self.value, self.prime_pow.get_top_modulus()[0])
        return ans

    cdef pAdicZZpXFMElement _rshift_c(self, long n, bint top_zeros):
        """
        sage: R = ZpFM(5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: ((w^5  >> 2) - w^3).valuation() >= 98
        True
        """
        if n < 0:
            return self._lshift_c(-n, top_zeros)
        elif n == 0:
            return self
        cdef pAdicZZpXFMElement ans = self._new_c()
        cdef Py_ssize_t i
        cdef long topcut, rem
        cdef ntl_ZZ holder
        if n < self.prime_pow.ram_prec_cap:
            ZZ_pX_eis_shift(ans.value, self.value, n, \
                            (<PowComputer_ZZ_pX_FM_Eis>self.prime_pow).low_shifter, \
                            (<PowComputer_ZZ_pX_FM_Eis>self.prime_pow).high_shifter, \
                            self.prime_pow.get_top_modulus()[0], \
                            self.prime_pow.small_powers[1], \
                            self.prime_pow.get_top_context().x)
            if top_zeros:
                holder = ntl_ZZ()
                rem = self.prime_pow.e - n % self.prime_pow.e
                topcut = self.prime_pow.prec_cap - n / self.prime_pow.e
                #print rem, topcut
                for i from 0 <= i <= ZZ_pX_deg(ans.value):
                    if i >= rem:
                        ZZ_rem(holder.x, ZZ_p_rep(ZZ_pX_coeff(ans.value, i)), self.prime_pow.pow_ZZ_tmp(topcut - 1)[0])
                    else:
                        ZZ_rem(holder.x, ZZ_p_rep(ZZ_pX_coeff(ans.value, i)), self.prime_pow.pow_ZZ_tmp(topcut)[0])
                    ZZ_pX_SetCoeff(ans.value, i, ZZ_to_ZZ_p(holder.x))
        return ans

    def add_bigoh(self, absprec):
        raise NotImplementedError

    cdef long valuation_c(self):
        cdef long valuation, index
        ZZ_pX_min_val_coeff(valuation, index, self.value, self.prime_pow.pow_ZZ_tmp(1)[0])
        if index == -1: # self == 0
            return self.prime_pow.ram_prec_cap
        return index + valuation * self.prime_pow.e

    def _set_uniformizerpow(self, n):
        cdef ntl_ZZ exp = ntl_ZZ(n)
        ZZ_pX_PowerXMod_pre(self.value, exp.x, self.prime_pow.get_top_modulus()[0])

    cdef ext_p_list(self, bint pos):
        ans = []
        cdef Integer list_elt
        cdef Integer halfp
        cdef Py_ssize_t i, j
        cdef ZZ_p_c const_term_holder
        self.prime_pow.restore_top_context()
        ZZ_p_construct(&const_term_holder)
        cdef ntl_ZZ holder = ntl_ZZ()
        halfp = self.prime_pow.prime._rshift(1)
        cdef pAdicZZpXFMElement shifter = self.copy()
        for i from 0 <= i < self.prime_pow.ram_prec_cap:
            list_elt = PY_NEW(Integer)
            ZZ_rem(holder.x, ZZ_p_rep(ZZ_pX_ConstTerm(shifter.value)), self.prime_pow.pow_ZZ_tmp(1)[0])
            ZZ_to_mpz(&list_elt.value, &holder.x)
            ans.append(list_elt)
            ZZ_p_sub(const_term_holder, ZZ_pX_ConstTerm(shifter.value), ZZ_to_ZZ_p(holder.x))
            ZZ_pX_SetCoeff(shifter.value, 0, const_term_holder)
            shifter = shifter._rshift_c(1, 0)
            #print shifter._ntl_rep()
        if not pos:
            # I know this isn't the best way to do this...
            for i from 0 <= i < self.prime_pow.ram_prec_cap:
                if ans[i] > halfp:
                    ans[i] -= self.prime_pow.prime
                    for j from i < j < self.prime_pow.ram_prec_cap:
                        ans[j] += 1
                        if ans[j] == self.prime_pow.prime:
                            ans[j] = 0
                        else:
                            break
        while len(ans) > 0:
            if ans[-1] == 0:
                ans.pop()
            else:
                break
        while len(ans) > 0:
            if ans[0] == 0:
                ans.pop(0)
            else:
                break
        ZZ_p_destruct(&const_term_holder)
        return ans

def make_EisensteinFixedModElement(parent, f):
    return EisensteinFixedModElement(parent, f)
