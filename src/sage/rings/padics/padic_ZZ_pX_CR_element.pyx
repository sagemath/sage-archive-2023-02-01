
include "../../ext/cdefs.pxi"
include "../../ext/stdsage.pxi"
include "../../ext/interrupt.pxi"

from sage.rings.integer cimport Integer
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext

cdef object infinity
from sage.rings.infinity import infinity

cdef class pAdicZZpXCRElement(pAdicZZpXElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, empty = False):
        """
        INPUT:
        parent -- either an EisensteinRingCappedRelative or UnramifiedRingCappedRelative
        x -- an ntl_ZZ_pX with modulus p^prec_cap
        """
        pAdicZZpXElement.__init__(self, parent)
        ZZ_pX_construct(&self.unit)
        if empty:
            self._normalized = 0
            return
        self._normalized = 1
        if relprec is not infinity and not PY_TYPE_CHECK(relprec, Integer):
            relprec = Integer(relprec)
        if (relprec is infinity) or (relprec > parent.precision_cap()):
            relprec = parent.precision_cap()
        if not absprec is infinity and not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        cdef mpz_t tmp
        cdef Py_ssize_t i
        if PY_TYPE_CHECK(x, pAdicGenericElement):
            if sefl.prime_pow.in_field == 0 and x.valuation() < 0:
                raise ValueError, "element has negative valuation"
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes."
        if PY_TYPE_CHECK(x, pAdicBaseGenericElement):
            mpz_init(tmp)
            (<pAdicBaseGenericElement>x)._set_to_mpz(tmp)
            self._set_from_mpz(tmp)
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
            if (<Integer>x.modulus())._is_power_of(<Integer>parent.prime()):
                x = x.lift()
            else:
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"
        cdef ntl_ZZ_pX poly
        if PY_TYPE_CHECK(x, Integer):
            self._set_from_mpz((<Integer>x).value)
        elif PY_TYPE_CHECK(x, Rational):
            self._set_from_mpq((<Rational>x).value)
        elif PY_TYPE_CHECK(x, ntl_ZZ_pX):
            self._set_from_ZZ_pX_c((<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c)
        elif PY_TYPE_CHECK(x, pAdicZZpXFMElement):
            self._set_from_ZZ_pX_c((<pAdicZZpXFMElement>x).value, ((<pAdicZZpXFMElement>x).prime_pow.get_top_context()))
        elif isinstance(x, list):
            poly = ntl_ZZ_pX(x, self.prime_pow.get_top_context())
            self._set_from_ZZ_pX_c(poly.x, poly.c)
        else:
            poly = ntl_ZZ_pX(x.list(), self.prime_pow.get_top_context())
            self._set_from_ZZ_pX_c(poly.x, poly.c)

    cdef int _set_from_mpz_rel(self, mpz_t x, unsigned long relprec) except -1:
        self.prime_pow.restore_top_context()
        cdef ZZ_c tmp
        ZZ_construct(&tmp)
        mpz_to_ZZ(&tmp, &x)
        ZZ_pX_SetCoeff(self.value, 0, ZZ_to_ZZ_p(tmp))
        ZZ_destruct(&tmp)

    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, unsigned long relprec) except -1:


    cdef int _set_from_ZZ_pX_c(self, ZZ_pX_c poly) except -1:
        self.prime_pow.restore_top_context()
        ZZ_pX_conv_modulus(self.value, poly, (<ntl_ZZ_pContext_class>self.prime_pow.get_top_context()).x)

    def __dealloc__(self):
        # Might need to call
        ZZ_pX_destruct(&self.value)

    cdef pAdicZZpXFMElement _new_c(self):
        raise NotImplementedError

    cdef RingElement _invert_c_impl(self):
        if self.valuation_c() > 0:
            raise ValueError, "cannot invert non-unit"
        cdef pAdicZZpXFMElement ans = self._new_c()
        _sig_on
        ZZ_pX_InvMod_newton(ans.value, self.value, self.prime_pow.get_top_modulus()[0], self.prime_pow.get_top_context().x, self.prime_pow.get_context(1).x)
        _sig_off
        return ans

    cdef pAdicZZpXFMElement _lshift_c(self, long n, bint top_zeros):
        raise NotImplementedError

    def __lshift__(pAdicZZpXFMElement self, shift):
        cdef pAdicZZpXFMElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            #Assuming that _new_c() initializes to zero.
            return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value), 1)

    cdef pAdicZZpXFMElement _rshift_c(self, long n, bint top_zeros):
        raise NotImplementedError

    def __rshift__(pAdicZZpXFMElement self, shift):
        cdef pAdicZZpXFMElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            #Assuming that _new_c() initializes to zero.
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value), 1)

    cdef ModuleElement _neg_c_impl(self):
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_negate(ans.value, self.value)
        return ans

    def __pow__(pAdicZZpXFMElement self, right, m): # m ignored
        if not PY_TYPE_CHECK(right, Integer):
            right = Integer(right)
        if not right and not self:
            raise ArithmeticError, "0^0 is undefined"
        cdef pAdicZZpXFMElement ans = self._new_c()
        cdef ntl_ZZ rZZ = PY_NEW(ntl_ZZ)
        mpz_to_ZZ(&rZZ.x, &(<Integer>right).value)
        _sig_on
        ZZ_pX_PowerMod_pre(ans.value, self.value, rZZ.x, self.prime_pow.get_top_modulus()[0])
        _sig_off
        return ans

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_add(ans.value, self.value, (<pAdicZZpXFMElement>right).value)
        return ans

    cdef RingElement _mul_c_impl(self, RingElement right):
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_MulMod_pre(ans.value, self.value, (<pAdicZZpXFMElement>right).value, self.prime_pow.get_top_modulus()[0])
        return ans

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_sub(ans.value, self.value, (<pAdicZZpXFMElement>right).value)
        return ans

    cdef RingElement _div_c_impl(self, RingElement right):
        if self.valuation_c() > 0:
            raise ValueError, "cannot invert non-unit"
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_PowerMod_long_pre(ans.value, (<pAdicZZpXFMElement>right).value, -1, self.prime_pow.get_top_modulus()[0])
        ZZ_pX_MulMod_pre(ans.value, self.value, ans.value, self.prime_pow.get_top_modulus()[0])
        return ans

    def copy(self):
        cdef pAdicZZpXFMElement ans = self._new_c()
        ans.value = self.value # does this actually copy correctly
        return ans

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def is_zero(self, absprec = None):
        cdef bint ans
        if absprec is None:
            ans = ZZ_pX_IsZero(self.value)
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0:
                    return True
                else:
                    ans = ZZ_pX_IsZero(self.value)
            else:
                ans = (self.valuation_c() >= mpz_get_si((<Integer>absprec).value))
        return ans

    def _ntl_rep(self):
        self.prime_pow.restore_top_context()
        cdef ntl_ZZ_pX ans = PY_NEW(ntl_ZZ_pX)
        ans.c = self.prime_pow.get_top_context()
        ans.x = self.value
        return ans

    def is_equal_to(self, right, absprec = None):
        # Should be sped up later
        return (self - right).is_zero(absprec)

    def lift(self):
        raise NotImplementedError

    def lift_to_precision(self, absprec):
        raise NotImplementedError

    def list(self, lift_mode = 'simple'):
        if lift_mode == 'simple':
            return self.parent().printer.ext_p_list(self, 1)
        elif lift_mode == 'smallest':
            return self.parent().printer.ext_p_list(self, 0)
        elif lift_mode == 'teichmuller':
            return self.teichmuller_list()
        else:
            raise ValueError, "lift mode must be one of 'simple', 'smallest' or 'teichmuller'"

    def teichmuller_list(self):
        raise NotImplementedError

    def _teichmuller_set(self):
        raise NotImplementedError

    def multiplicative_order(self):
        raise NotImplementedError

    def padded_list(self, n, lift_mode = 'simple'):
        raise NotImplementedError

    def precision_absolute(self):
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.prime_pow.ram_prec_cap)
        return ans

    def precision_relative(self):
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.prime_pow.ram_prec_cap - self.valuation_c())
        return ans

    def residue(self, n):
        raise NotImplementedError

    def unit_part(self, top_zeros = True):
        return self._unit_part_c(top_zeros)

    cdef pAdicZZpXFMElement _unit_part_c(self, bint top_zeros):
        return self._rshift_c(self.valuation_c(), top_zeros)

