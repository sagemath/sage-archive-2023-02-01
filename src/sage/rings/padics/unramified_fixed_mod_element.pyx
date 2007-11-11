
include "../../libs/ntl/decl.pxi"
include "../../ext/cdefs.pxi"
include "../../ext/stdsage.pxi"

from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.rings.integer cimport Integer
from sage.structure.element cimport Element
from sage.rings.padics.padic_printing cimport pAdicPrinter_class

cdef class UnramifiedFixedModElement(pAdicZZpXFMElement):
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
        return make_UnramifiedFixedModElement, (self.parent(), holder)

    cdef pAdicZZpXFMElement _new_c(self):
        self.prime_pow.restore_top_context()
        cdef UnramifiedFixedModElement ans = PY_NEW(UnramifiedFixedModElement)
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
        ZZ_pX_left_pshift(ans.value, self.value, self.prime_pow.pow_ZZ_tmp(n)[0], self.prime_pow.get_top_context().x)
        return ans

    cdef pAdicZZpXFMElement _rshift_c(self, long n, bint top_zeros):
        if n < 0:
            return self._lshift_c(-n, top_zeros)
        elif n == 0:
            return self
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_right_pshift(ans.value, self.value, self.prime_pow.pow_ZZ_tmp(n)[0], self.prime_pow.get_top_context().x)
        return ans

    def add_bigoh(self, absprec):
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        cdef pAdicZZpXFMElement ans = self._new_c()
        if mpz_fits_ulong_p((<Integer>absprec).value) == 0:
            if mpz_sgn((<Integer>absprec).value) < 0:
                return ans # assumes _new_c() initializes to 0
            return self # absprec > prec_cap
        cdef unsigned long aprec = mpz_get_ui((<Integer>absprec).value)
        if aprec >= self.prime_pow.prec_cap:
            return self
        cdef ZZ_pX_c tmp
        cdef ntl_ZZ_pContext_class c = self.prime_pow.get_context(aprec)
        c.restore_c()
        ZZ_pX_construct(&tmp)
        ZZ_pX_conv_modulus(tmp, self.value, c.x)
        ZZ_pX_conv_modulus(ans.value, tmp, (<ntl_ZZ_pContext_class>self.prime_pow.get_top_context()).x)
        ZZ_pX_destruct(&tmp)
        return ans

    cdef long valuation_c(self):
        cdef long valuation, index
        ZZ_pX_min_val_coeff(valuation, index, self.value, self.prime_pow.pow_ZZ_tmp(1)[0])
        return valuation

    cdef ext_p_list(self, bint pos):
        ans = []
        cdef Py_ssize_t i, j
        cdef ntl_ZZ ZZ_coeff = ntl_ZZ()
        cdef Integer coeff = PY_NEW(Integer)
        cdef Integer zero = Integer(0)
        cdef pAdicPrinter_class printer = <pAdicPrinter_class>self.parent()._printer
        for j from 0 <= j < self.prime_pow.prec_cap:
            ans.append([])
        for i from 0 <= i < self.prime_pow.deg:
            ZZ_coeff.x = ZZ_p_rep(ZZ_pX_coeff(self.value, i))
            ZZ_to_mpz(&coeff.value, &ZZ_coeff.x)
            L = printer.base_p_list(coeff.value, pos)
            for j from 0 <= j < self.prime_pow.prec_cap:
                if j < len(L):
                    ans[j].append(L[j])
                else:
                    ans[j].append(zero)
        for j from 0 <= j < self.prime_pow.prec_cap:
            while len(ans[j]) > 0:
                if ans[j][-1] == 0:
                    ans[j].pop()
                else:
                    break
        while len(ans) > 0:
            if ans[-1] == []:
                ans.pop()
            else:
                break
        while len(ans) > 0:
            if ans[0] == []:
                ans.pop(0)
            else:
                break
        return ans


def make_UnramifiedFixedModElement(parent, f):
    return UnramifiedFixedModElement(parent, f)



