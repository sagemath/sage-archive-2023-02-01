
include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

from sage.rings.padics.padic_fixed_mod_element cimport pAdicFixedModElement
from sage.rings.padics.padic_capped_absolute_element cimport pAdicCappedAbsoluteElement
from sage.rings.padics.padic_capped_relative_element cimport pAdicCappedRelativeElement
from sage.rings.padics.eisenstein_fixed_mod_element cimport EisensteinFixedModElement
from sage.rings.padics.unramified_fixed_mod_element cimport UnramifiedFixedModElement

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX

from sage.rings.padics.padic_ring_fixed_mod import pAdicRingFixedMod
from sage.rings.padics.padic_ring_capped_absolute import pAdicRingCappedAbsolute
from sage.rings.padics.padic_ring_capped_relative import pAdicRingCappedRelative
from sage.rings.padics.padic_field_capped_relative import pAdicFieldCappedRelative
from sage.rings.padics.padic_extension_leaves import EisensteinExtensionRingFixedMod, UnramifiedExtensionRingFixedMod
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.homset import RingHomset

cdef pAdicFixedModElement make_new_ZpFM(parent):
    cdef pAdicFixedModElement ans = PY_NEW(pAdicFixedModElement)
    ans._parent = parent
    ans.prime_pow = <PowComputer_class>parent.prime_pow
    mpz_init(ans.value)
    return ans

cdef pAdicCappedAbsoluteElement make_new_ZpCA(parent):
    cdef pAdicCappedAbsoluteElement ans = PY_NEW(pAdicCappedAbsoluteElement)
    ans._parent = parent
    ans.prime_pow = <PowComputer_class>parent.prime_pow
    mpz_init(ans.value)
    return ans

cdef pAdicCappedRelativeElement make_new_ZQpCR(parent):
    cdef pAdicCappedRelativeElement ans = PY_NEW(pAdicCappedRelativeElement)
    ans._parent = parent
    ans.prime_pow = <PowComputer_class>parent.prime_pow
    mpz_init(ans.unit)
    return ans

cdef EisensteinFixedModElement make_new_EisFM(parent):
    cdef PowComputer_ZZ_pX ppow = <PowComputer_ZZ_pX>parent.prime_pow
    ppow.restore_top_context()
    cdef EisensteinFixedModElement ans = PY_NEW(EisensteinFixedModElement)
    ans._parent = parent
    ans.prime_pow = ppow
    ZZ_pX_construct(&ans.value)
    return ans

cdef UnramifiedFixedModElement make_new_UnrFM(parent):
    cdef PowComputer_ZZ_pX ppow = <PowComputer_ZZ_pX>parent.prime_pow
    ppow.restore_top_context()
    cdef UnramifiedFixedModElement ans = PY_NEW(UnramifiedFixedModElement)
    ans._parent = parent
    ans.prime_pow = ppow
    ZZ_pX_construct(&ans.value)
    return ans

cdef long_min(long a, long b):
    if a <= b:
        return a
    else:
        return b

cdef class Section_ZpFM_ZZ(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be pAdicFixedModElement
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, (<pAdicFixedModElement>x).value)
        return ans

cdef class Morphism_ZZ_ZpFM(RingHomomorphism_coercion):
    def __init__(self, ZpFM):
        if not isinstance(ZpFM, pAdicRingFixedMod):
            raise TypeError, "Must be a fixed mod ring"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZZ, ZpFM), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be an Integer
        """
        cdef pAdicFixedModElement ans = make_new_ZpFM(self._codomain)
        ans._set_from_mpz((<Integer>x).value)
        return ans

cdef class Section_ZpCA_ZZ(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be pAdicCappedAbsoluteElement
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, (<pAdicCappedAbsoluteElement>x).value)
        return ans

cdef class Morphism_ZZ_ZpCA(RingHomomorphism_coercion):
    def __init__(self, ZpCA):
        if not isinstance(ZpCA, pAdicRingCappedAbsolute):
            raise TypeError, "Must be a capped absolute ring"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZZ, ZpCA), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be an Integer
        """
        cdef pAdicCappedAbsoluteElement ans = make_new_ZpCA(self._codomain)
        ans._set_from_mpz_abs((<Integer>x).value, ans.prime_pow.prec_cap)
        return ans

cdef class Section_ZpCR_ZZ(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be pAdicCappedRelativeElement
        """
        cdef Integer ans = PY_NEW(Integer)
        _sig_on
        cdef pAdicCappedRelativeElement _x = <pAdicCappedRelativeElement>x
        _sig_off
        mpz_set(ans.value, _x.unit)
        if mpz_sgn(ans.value) == -1:
            mpz_set_ui(ans.value, 0)
            return ans
        elif mpz_sgn(ans.value) == 0:
            mpz_set_ui(ans.value, 1)
        mpz_mul(ans.value, ans.value, _x.prime_pow.pow_mpz_t_tmp(_x.valuation_c())[0])
        return ans

cdef class Morphism_ZZ_ZpCR(RingHomomorphism_coercion):
    def __init__(self, ZpCR):
        if not isinstance(ZpCR, pAdicRingCappedRelative):
            raise TypeError, "Must be a capped relative ring"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZZ, ZpCR), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be an Integer
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        ans._set_from_mpz_rel((<Integer>x).value, ans.prime_pow.prec_cap)
        return ans

cdef class Section_QpCR_ZZ(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be pAdicCappedRelativeElement
        """
        cdef Integer ans = PY_NEW(Integer)
        (<pAdicCappedRelativeElement>x)._set_to_mpz(ans.value)
        return ans

cdef class Morphism_ZZ_QpCR(RingHomomorphism_coercion):
    def __init__(self, QpCR):
        if not isinstance(QpCR, pAdicFieldCappedRelative):
            raise TypeError, "Must be a capped relative field"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZZ, QpCR), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be an Integer
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        ans._set_from_mpz_rel((<Integer>x).value, ans.prime_pow.prec_cap)
        return ans

cdef class Section_QpCR_QQ(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be pAdicCappedRelativeElement
        """
        cdef Rational ans = PY_NEW(Rational)
        (<pAdicCappedRelativeElement>x)._set_to_mpq(ans.value)
        return ans

cdef class Morphism_QQ_QpCR(RingHomomorphism_coercion):
    def __init__(self, QpCR):
        if not isinstance(QpCR, pAdicFieldCappedRelative):
            raise TypeError, "Must be a capped relative field"
        RingHomomorphism_coercion.__init__(self, RingHomset(QQ, QpCR), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a Rational
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        ans._set_from_mpq_rel((<Rational>x).value, self.prime_pow.prec_cap)
        return ans

cdef class Section_QpCR_ZpCR(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedRelativeElement
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        ans._set_from_CR(x)
        return ans

cdef class Morphism_ZpCR_QpCR(RingHomomorphism_coercion):
    def __init__(self, ZpCR, QpCR):
        if not isinstance(ZpCR, pAdicRingCappedRelative) or not isinstance(QpCR, pAdicFieldCappedRelative):
            raise TypeError, "must be appropriate types"
        if ZpCR.prime() != QpCR.prime():
            raise TypeError, "primes must be the same"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZpCR, QpCR), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedRelativeElement
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        ans._set_from_CR(x)
        return ans

cdef class Section_ZpCA_ZpCR(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedAbsoluteElement
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        cdef pAdicCappedAbsoluteElement _x = <pAdicCappedAbsoluteElement>x
        ans._set_from_mpz_both(_x.value, _x.absprec, ans.prime_pow.prec_cap)
        return ans

cdef class Morphism_ZpCR_ZpCA(RingHomomorphism_coercion):
    def __init__(self, ZpCR, ZpCA):
        if not isinstance(ZpCR, pAdicRingCappedRelative) or not isinstance(ZpCA, pAdicRingCappedAbsolute):
            raise TypeError, "rings of the wrong types"
        if ZpCR.prime() != ZpCA.prime():
            raise TypeError, "primes must be the same"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZpCR, ZpCA), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedRelativeElement
        """
        cdef pAdicCappedAbsoluteElement ans = make_new_ZpCA(self._codomain)
        cdef pAdicCappedRelativeElement _x = <pAdicCappedRelativeElement>x
        cdef mpz_t tmp
        mpz_init(tmp)
        if _x.is_exact_zero():
            mpz_set_ui(tmp, 0)
            ans._set_from_mpz_abs(tmp, ans.prime_pow.prec_cap)
        else:
            _x._set_to_mpz(tmp)
            ans._set_from_mpz_abs(tmp, long_min(ans.prime_pow.prec_cap, _x.valuation_c() + _x.ordp))
        mpz_clear(tmp)
        return ans

cdef class Section_QpCR_ZpCA(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedRelativeElement
        """
        cdef pAdicCappedAbsoluteElement ans = make_new_ZpCA(self._codomain)
        cdef pAdicCappedRelativeElement _x = <pAdicCappedRelativeElement>x
        cdef mpz_t tmp
        mpz_init(tmp)
        if _x.is_exact_zero():
            mpz_set_ui(tmp, 0)
            ans._set_from_mpz_abs(tmp, ans.prime_pow.prec_cap)
        else:
            _x._set_to_mpz(tmp)
            ans._set_from_mpz_abs(tmp, long_min(ans.prime_pow.prec_cap, _x.valuation_c() + _x.ordp))
        mpz_clear(tmp)
        return ans


cdef class Morphism_ZpCA_QpCR(RingHomomorphism_coercion):
    def __init__(self, ZpCA, QpCR):
        if not isinstance(ZpCA, pAdicRingCappedAbsolute) or not isinstance(QpCR, pAdicFieldCappedRelative):
            raise TypeError, "rings of the wrong types"
        if QpCR.prime() != ZpCA.prime():
            raise TypeError, "primes must be the same"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZpCA, QpCR), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedAbsoluteElement
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        cdef pAdicCappedAbsoluteElement _x = <pAdicCappedAbsoluteElement>x
        ans._set_from_mpz_both(_x.value, _x.absprec, ans.prime_pow.prec_cap)
        return ans

cdef class Section_ZpFM_ZpCA(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicFixedModElement
        """
        cdef pAdicCappedAbsoluteElement ans = make_new_ZpCA(self._codomain)
        cdef pAdicFixedModElement _x = <pAdicFixedModElement>x
        ans._set_from_mpz_abs(_x.value, long_min(_x.prime_pow.prec_cap, ans.prime_pow.prec_cap))
        return ans

cdef class Morphism_ZpCA_ZpFM(RingHomomorphism_coercion):
    def __init__(self, ZpCA, ZpFM):
        if not isinstance(ZpCA, pAdicRingCappedAbsolute) or not isinstance(ZpFM, pAdicRingFixedMod):
            raise TypeError, "rings of the wrong types"
        if ZpFM.prime() != ZpCA.prime():
            raise TypeError, "primes must be the same"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZpCA, ZpFM), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedAbsoluteElement
        """
        cdef pAdicFixedModElement ans = make_new_ZpFM(self._codomain)
        ans._set_from_mpz((<pAdicCappedAbsoluteElement>x).value)
        return ans

cdef class Section_ZpFM_ZpCR(Section):
    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicFixedModElement
        """
        cdef pAdicCappedRelativeElement ans = make_new_ZQpCR(self._codomain)
        cdef pAdicFixedModElement _x = <pAdicFixedModElement>x
        ans._set_from_mpz_both(_x.value, _x.prime_pow.prec_cap, ans.prime_pow.prec_cap)
        return ans

cdef class Morphism_ZpCR_ZpFM(RingHomomorphism_coercion):
    def __init__(self, ZpCR, ZpFM):
        if not isinstance(ZpCR, pAdicRingCappedRelative) or not isinstance(ZpFM, pAdicRingFixedMod):
            raise TypeError, "rings of the wrong types"
        if ZpFM.prime() != ZpCR.prime():
            raise TypeError, "primes must be the same"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZpCR, ZpFM), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicCappedRelativeElement
        """
        cdef pAdicFixedModElement ans = make_new_ZpFM(self._codomain)
        cdef pAdicCappedRelativeElement _x = <pAdicCappedRelativeElement>x
        if _x.ordp >= ans.prime_pow.prec_cap:
            mpz_set_ui(ans.value, 0)
        else:
            ans._set_from_mpz(_x.unit)
            _sig_on
            mpz_mul(ans.value, ans.value, _x.prime_pow.pow_mpz_t_tmp(_x.ordp)[0])
            mpz_mod(ans.value, ans.value, ans.prime_pow.pow_mpz_t_top()[0])
            _sig_off



## Extension types
cdef class Morphism_ZZ_EisFM(RingHomomorphism_coercion):
    def __init__(self, EisFM):
        if not isinstance(EisFM, EisensteinExtensionRingFixedMod):
            raise TypeError, "ring of the wrong type"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZZ, EisFM), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be an Integer.
        """
        cdef EisensteinFixedModElement ans = make_new_EisFM(self._codomain)
        ans._set_from_mpz((<Integer>x).value)
        return ans

cdef class Morphism_ZZ_UnrFM(RingHomomorphism_coercion):
    def __init__(self, UnrFM):
        if not isinstance(UnrFM, UnramifiedExtensionRingFixedMod):
            raise TypeError, "ring of the wrong type"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZZ, UnrFM), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be an Integer.
        """
        cdef UnramifiedFixedModElement ans = make_new_UnrFM(self._codomain)
        ans._set_from_mpz((<Integer>x).value)
        return ans

cdef class Morphism_ZpFM_EisFM(RingHomomorphism_coercion):
    def __init__(self, ZpFM, EisFM):
        if not isinstance(EisFM, EisensteinExtensionRingFixedMod) or not isinstance(ZpFM, pAdicRingFixedMod):
            raise TypeError, "rings of the wrong types"
        if ZpFM.prime() != EisFM.prime():
            raise TypeError, "must have the same prime"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZpFM, EisFM), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicFixedModElement
        """
        cdef EisensteinFixedModElement ans = make_new_EisFM(self._codomain)
        ans._set_from_mpz((<pAdicFixedModElement>x).value)
        return ans

cdef class Morphism_ZpFM_UnrFM(RingHomomorphism_coercion):
    def __init__(self, ZpFM, UnrFM):
        if not isinstance(UnrFM, UnramifiedExtensionRingFixedMod) or not isinstance(ZpFM, pAdicRingFixedMod):
            raise TypeError, "rings of the wrong types"
        if ZpFM.prime() != UnrFM.prime():
            raise TypeError, "must have the same prime"
        RingHomomorphism_coercion.__init__(self, RingHomset(ZpFM, UnrFM), False)

    cdef Element _call_c_impl(self, Element x):
        """
        x must be a pAdicFixedModElement
        """
        cdef EisensteinFixedModElement ans = make_new_UnrFM(self._codomain)
        ans._set_from_mpz((<pAdicFixedModElement>x).value)
        return ans

