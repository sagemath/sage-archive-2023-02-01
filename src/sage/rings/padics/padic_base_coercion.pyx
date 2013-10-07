"""
This file contains morphisms implementing coercion between base p-adic rings and fields as well as ZZ and QQ

AUTHOR:
-- David Roe (initial version: 2010-8-25)
"""

include "sage/ext/gmp.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.structure.element cimport Element
from sage.rings.padics.padic_capped_relative_element cimport pAdicCappedRelativeElement
from sage.rings.padics.padic_capped_absolute_element cimport pAdicCappedAbsoluteElement
from sage.rings.padics.padic_fixed_mod_element cimport pAdicFixedModElement
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps

cdef _process_args_and_kwds(args, kwds):
    """
    This function obtains values for absprec and relprec from a combination of positional and keyword arguments.
    """
    if len(args) > 3:
        raise TypeError, "too many positional arguments"
    if len(args) == 3:
        if kwds.has_key("empty"):
            raise TypeError, "_call_with_args() got multiple values for keyword argument 'empty'"
        if args[2]: return True
    elif kwds.has_key("empty") and kwds["empty"]:
        return True
    if len(args) >= 2:
        if kwds.has_key("relprec"):
            raise TypeError, "_call_with_args() got multiple values for keyword argument 'relprec'"
        relprec = args[1]
    elif kwds.has_key("relprec"):
        relprec = kwds["relprec"]
    else:
        relprec = infinity
    if relprec is not infinity and not PY_TYPE_CHECK(relprec, Integer):
        relprec = Integer(relprec)
    if len(args) >= 1:
        if kwds.has_key("absprec"):
            raise TypeError, "_call_with_args() got multiple values for keyword argument 'absprec'"
        absprec = args[0]
    elif kwds.has_key("absprec"):
        absprec = kwds["absprec"]
    else:
        absprec = infinity
    if absprec is not infinity and not PY_TYPE_CHECK(absprec, Integer):
        absprec = Integer(absprec)
    return (absprec, relprec)

cdef class pAdicCoercion_ZZ_CR(RingHomomorphism_coercion):
    """
    The canonical inclusion from ZZ to either ZpCR or QpCR.

    EXAMPLES::

        sage: f = Zp(5).coerce_map_from(ZZ); f
        Ring Coercion morphism:
          From: Integer Ring
          To:   5-adic Ring with capped relative precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicCoercion_ZZ_CR'>
        """
        RingHomomorphism_coercion.__init__(self, ZZ.Hom(R), check=False)
        self.prime_pow = R.prime_pow
        self._section = pAdicConvert_CR_ZZ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Integer Ring
              To:   5-adic Ring with capped relative precision 20
            sage: g == f
            True
            sage: g is f
            False
            sage: g(5)
            5 + O(5^21)
            sage: g(5) == f(5)
            True

        """
        _slots['prime_pow'] = self.prime_pow
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Integer Ring
              To:   5-adic Ring with capped relative precision 20
            sage: g == f
            True
            sage: g is f
            False
            sage: g(5)
            5 + O(5^21)
            sage: g(5) == f(5)
            True

        """
        self.prime_pow = _slots['prime_pow']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring with capped relative precision 20
            sage: f(5)
            5 + O(5^21)
        """
        cdef pAdicCappedRelativeElement ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        ans._set_from_mpz_rel((<Integer>x).value, self.prime_pow.prec_cap)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both), or an empty element is desired.

        See the documentation for pAdicCappedRelativeElement.__init__ for more details.

        EXAMPLES::

            sage: R = Zp(5,4)
            sage: type(R(10,2))
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: R(10,2)
            2*5 + O(5^2)
            sage: R(10,3,1)
            2*5 + O(5^2)
            sage: R(10,absprec=2)
            2*5 + O(5^2)
            sage: R(10,relprec=2)
            2*5 + O(5^3)
            sage: R(10,absprec=1)
            O(5)
            sage: R(10,empty=True)
            O(5^0)
        """
        cdef pAdicCappedRelativeElement ans
        processed = _process_args_and_kwds(args, kwds)
        ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        if processed is True:
            ans._normalized = 0
            return ans
        ans._set_from_Integer(x, processed[0], processed[1])
        return ans

    def section(self):
        """
        Returns a map back to ZZ that approximates an element of Zp or Qp by an integer.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ).section()
            sage: f(Zp(5)(-1)) - 5^20
            -1
        """
        return self._section

cdef class pAdicConvert_CR_ZZ(RingMap):
    """
    The map from Zp or Qp back to ZZ that returns the the smallest non-negative integer approximation to its input which is accurate up to the precision.

    If the input has negative valuation, raises a ValueError.

    EXAMPLES::

        sage: f = Zp(5).coerce_map_from(ZZ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Ring with capped relative precision 20
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(ZZ).section(); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicConvert_CR_ZZ'>
            sage: f.category()
            Category of hom sets in Category of sets with partial maps
            sage: Zp(5).coerce_map_from(ZZ).section().category()
            Category of hom sets in Category of sets
        """
        if R.is_field():
            RingMap.__init__(self, Hom(R, ZZ, SetsWithPartialMaps()))
        else:
            RingMap.__init__(self, Hom(R, ZZ, Sets()))

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(ZZ).section()
            sage: f(Qp(5)(-1)) - 5^20
            -1
            sage: f(Qp(5)(0))
            0
            sage: f(Qp(5)(1/5))
            Traceback (most recent call last):
            ...
            ValueError: negative valuation
        """
        cdef Integer ans = PY_NEW(Integer)
        (<pAdicCappedRelativeElement>x)._set_mpz_into(ans.value)
        return ans

cdef class pAdicCoercion_QQ_CR(RingHomomorphism_coercion):
    """
    The canonical inclusion from QQ to QpCR.

    EXAMPLES::

        sage: f = copy(Qp(5).coerce_map_from(QQ)); f
        Ring Coercion morphism:
          From: Rational Field
          To:   5-adic Field with capped relative precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicCoercion_QQ_CR'>
        """
        RingHomomorphism_coercion.__init__(self, QQ.Hom(R), check=False)
        self.prime_pow = R.prime_pow
        self._section = pAdicConvert_CR_QQ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ)
            sage: g = copy(f)    # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Rational Field
              To:   5-adic Field with capped relative precision 20
            sage: g is f
            False
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True

        """
        _slots['prime_pow'] = self.prime_pow
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Rational Field
              To:   5-adic Field with capped relative precision 20
            sage: g is f
            False
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True

        """
        self.prime_pow = _slots['prime_pow']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ)
            sage: f(0).parent()
            5-adic Field with capped relative precision 20
            sage: f(1/5)
            5^-1 + O(5^19)
            sage: f(1/4)
            4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + 3*5^10 + 3*5^11 + 3*5^12 + 3*5^13 + 3*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 3*5^18 + 3*5^19 + O(5^20)
        """
        cdef pAdicCappedRelativeElement ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        ans._set_from_mpq_rel((<Rational>x).value, self.prime_pow.prec_cap)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both), or an empty element is desired.

        See the documentation for pAdicCappedRelativeElement.__init__ for more details.


        EXAMPLES::

            sage: R = Qp(5,4)
            sage: type(R(10/3,2))
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: R(10/3,2)
            4*5 + O(5^2)
            sage: R(10/3,3,1)
            4*5 + O(5^2)
            sage: R(10/3,absprec=2)
            4*5 + O(5^2)
            sage: R(10/3,relprec=2)
            4*5 + 5^2 + O(5^3)
            sage: R(10/3,absprec=1)
            O(5)
            sage: R(10/3,empty=True)
            O(5^0)
            sage: R(3/100,absprec=-1)
            2*5^-2 + O(5^-1)
        """
        cdef pAdicCappedRelativeElement ans
        processed = _process_args_and_kwds(args, kwds)
        ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        if processed is True:
            ans._normalized = 0
            return ans
        ans._set_from_Rational(x, processed[0], processed[1])
        return ans

    def section(self):
        """
        Returns a map back to QQ that approximates an element of Qp by a rational number.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ).section()
            sage: f(Qp(5)(1/4))
            1/4
            sage: f(Qp(5)(1/5))
            1/5
        """
        return self._section

cdef class pAdicConvert_CR_QQ(RingMap):
    """
    The map from Qp back to QQ that returns the the smallest non-negative integer approximation to its input which is accurate up to the precision, divided by a power of p if the valuation of the input is negative.

    EXAMPLES::

        sage: f = Qp(5).coerce_map_from(QQ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Field with capped relative precision 20
          To:   Rational Field
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ).section(); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicConvert_CR_QQ'>
            sage: f.category()
            Category of hom sets in Category of sets
        """
        RingMap.__init__(self, Hom(R, QQ, Sets()))

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ).section()
            sage: f(Qp(5)(-1))
            -1
            sage: f(Qp(5)(0))
            0
            sage: f(Qp(5)(1/5))
            1/5
        """
        cdef Rational ans = PY_NEW(Rational)
        cdef pAdicCappedRelativeElement x = <pAdicCappedRelativeElement> _x
        if mpz_sgn(x.unit) <= 0:
            mpq_set_ui(ans.value, 0, 1)
        else:
            mpq_rational_reconstruction(ans.value, x.unit, x.prime_pow.pow_mpz_t_tmp(x.relprec)[0])
            if x.ordp > 0:
                mpz_mul(mpq_numref(ans.value), mpq_numref(ans.value), x.prime_pow.pow_mpz_t_tmp(x.ordp)[0])
            elif x.ordp < 0:
                mpz_mul(mpq_denref(ans.value), mpq_denref(ans.value), x.prime_pow.pow_mpz_t_tmp(-x.ordp)[0])
        return ans

cdef class pAdicConvert_QQ_CR(Morphism):
    """
    The inclusion map from QQ to Zp that is defined on all elements with non-negative p-adic valuation.

    EXAMPLES::

        sage: f = Zp(5).convert_map_from(QQ); f
        Generic morphism:
          From: Rational Field
          To:   5-adic Ring with capped relative precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Zp(5).convert_map_from(QQ); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicConvert_QQ_CR'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self.prime_pow = R.prime_pow
        self._section = pAdicConvert_CR_QQ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).convert_map_from(QQ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f       # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True

        """
        _slots['prime_pow'] = self.prime_pow
        _slots['_section'] = self._section
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).convert_map_from(QQ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f       # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True

        """
        self.prime_pow = _slots['prime_pow']
        self._section = _slots['_section']
        Morphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Zp(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: f(0)
            0
        """
        cdef pAdicCappedRelativeElement ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        ans._set_from_mpq_rel((<Rational>x).value, self.prime_pow.prec_cap)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both), or an empty element is desired.

        See the documentation for pAdicCappedRelativeElement.__init__ for more details.

        EXAMPLES::

            sage: R = Zp(5,4)
            sage: type(R(10/3,2))
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: R(10/3,2)
            4*5 + O(5^2)
            sage: R(10/3,3,1)
            4*5 + O(5^2)
            sage: R(10/3,absprec=2)
            4*5 + O(5^2)
            sage: R(10/3,relprec=2)
            4*5 + 5^2 + O(5^3)
            sage: R(10/3,absprec=1)
            O(5)
            sage: R(10/3,empty=True)
            O(5^0)
            sage: R(3/100,relprec=3)
            Traceback (most recent call last):
            ...
            ValueError: p divides the denominator
        """
        cdef pAdicCappedRelativeElement ans
        processed = _process_args_and_kwds(args, kwds)
        ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        if processed is True:
            ans._normalized = 0
            return ans
        ans._set_from_Rational(x, processed[0], processed[1])
        return ans

    def section(self):
        """
        Returns the map from Zp back to QQ that returns the smallest non-negative integer approximation to its input which is accurate up to the precision.

        EXAMPLES::

            sage: f = Zp(5,4).convert_map_from(QQ).section()
            sage: f(Zp(5,4)(-1))
            -1
        """
        return self._section

cdef class pAdicCoercion_ZZ_CA(RingHomomorphism_coercion):
    """
    The canonical inclusion from ZZ to ZpCA.

    EXAMPLES::

        sage: f = ZpCA(5).coerce_map_from(ZZ); f
        Ring Coercion morphism:
          From: Integer Ring
          To:   5-adic Ring with capped absolute precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicCoercion_ZZ_CA'>
        """
        RingHomomorphism_coercion.__init__(self, ZZ.Hom(R), check=False)
        self.prime_pow = R.prime_pow
        self._section = pAdicConvert_CA_ZZ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: f(6) == g(6)
            True

        """
        _slots['prime_pow'] = self.prime_pow
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: f(6) == g(6)
            True

        """
        self.prime_pow = _slots['prime_pow']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring with capped absolute precision 20
            sage: f(5)
            5 + O(5^20)
        """
        cdef pAdicCappedAbsoluteElement ans = PY_NEW(pAdicCappedAbsoluteElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        ans._set_from_mpz_abs((<Integer>x).value, self.prime_pow.prec_cap)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both), or an empty element is desired.

        See the documentation for pAdicCappedAbsoluteElement.__init__ for more details.

        EXAMPLES::

            sage: R = ZpCA(5,4)
            sage: type(R(10,2))
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(10,2)
            2*5 + O(5^2)
            sage: R(10,3,1)
            2*5 + O(5^2)
            sage: R(10,absprec=2)
            2*5 + O(5^2)
            sage: R(10,relprec=2)
            2*5 + O(5^3)
            sage: R(10,absprec=1)
            O(5)
            sage: R(10,empty=True)
            O(5^0)
        """
        cdef pAdicCappedAbsoluteElement ans
        processed = _process_args_and_kwds(args, kwds)
        ans = PY_NEW(pAdicCappedAbsoluteElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        if processed is True:
            return ans
        ans._set_from_Integer(x, processed[0], processed[1])
        return ans

    def section(self):
        """
        Returns a map back to ZZ that approximates an element of Zp by an integer.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ).section()
            sage: f(ZpCA(5)(-1)) - 5^20
            -1
        """
        return self._section

cdef class pAdicConvert_CA_ZZ(RingMap):
    """
    The map from ZpCA back to ZZ that returns the the smallest non-negative integer approximation to its input which is accurate up to the precision.

    EXAMPLES::

        sage: f = ZpCA(5).coerce_map_from(ZZ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Ring with capped absolute precision 20
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ).section(); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicConvert_CA_ZZ'>
            sage: f.category()
            Category of hom sets in Category of sets
        """
        RingMap.__init__(self, Hom(R, ZZ, Sets()))

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ).section()
            sage: f(ZpCA(5)(-1)) - 5^20
            -1
            sage: f(ZpCA(5)(0))
            0
        """
        cdef Integer ans = PY_NEW(Integer)
        (<pAdicCappedAbsoluteElement>x)._set_mpz_into(ans.value)
        return ans

cdef class pAdicConvert_QQ_CA(Morphism):
    """
    The inclusion map from QQ to ZpCA that is defined on all elements with non-negative p-adic valuation.

    EXAMPLES::

        sage: f = ZpCA(5).convert_map_from(QQ); f
        Generic morphism:
          From: Rational Field
          To:   5-adic Ring with capped absolute precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpCA(5).convert_map_from(QQ); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicConvert_QQ_CA'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self.prime_pow = R.prime_pow

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).convert_map_from(QQ)
            sage: g = copy(f) # indirect doctest
            sage: g == f      # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True

        """
        _slots['prime_pow'] = self.prime_pow
        try:
            _slots['_section'] = self._section
        except AttributeError:
            pass
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).convert_map_from(QQ)
            sage: g = copy(f) # indirect doctest
            sage: g == f      # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True

        """
        self.prime_pow = _slots['prime_pow']
        try:
            self._section = _slots['_section']
        except KeyError:
            pass
        Morphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpCA(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: f(0)
            O(5^4)
        """
        cdef pAdicCappedAbsoluteElement ans = PY_NEW(pAdicCappedAbsoluteElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        ans._set_from_mpq_abs((<Rational>x).value, self.prime_pow.prec_cap)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both), or an empty element is desired.

        See the documentation for pAdicCappedAbsoluteElement.__init__ for more details.

        EXAMPLES::

            sage: R = ZpCA(5,4)
            sage: type(R(10/3,2))
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(10/3,2)
            4*5 + O(5^2)
            sage: R(10/3,3,1)
            4*5 + O(5^2)
            sage: R(10/3,absprec=2)
            4*5 + O(5^2)
            sage: R(10/3,relprec=2)
            4*5 + 5^2 + O(5^3)
            sage: R(10/3,absprec=1)
            O(5)
            sage: R(10/3,empty=True)
            O(5^0)
            sage: R(3/100,relprec=3)
            Traceback (most recent call last):
            ...
            ValueError: p divides denominator
        """
        cdef pAdicCappedAbsoluteElement ans
        processed = _process_args_and_kwds(args, kwds)
        ans = PY_NEW(pAdicCappedAbsoluteElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        if processed is True:
            return ans
        ans._set_from_Rational(x, processed[0], processed[1])
        return ans

cdef class pAdicCoercion_ZZ_FM(RingHomomorphism_coercion):
    """
    The canonical inclusion from ZZ to ZpFM.

    EXAMPLES::

        sage: f = ZpFM(5).coerce_map_from(ZZ); f
        Ring Coercion morphism:
          From: Integer Ring
          To:   5-adic Ring of fixed modulus 5^20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicCoercion_ZZ_FM'>
        """
        RingHomomorphism_coercion.__init__(self, ZZ.Hom(R), check=False)
        self.prime_pow = R.prime_pow
        self._section = pAdicConvert_FM_ZZ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True

        """
        _slots['prime_pow'] = self.prime_pow
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True

        """
        self.prime_pow = _slots['prime_pow']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring of fixed modulus 5^20
            sage: f(5)
            5 + O(5^20)
        """
        cdef pAdicFixedModElement ans = PY_NEW(pAdicFixedModElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        ans._set_from_mpz((<Integer>x).value)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both), or an empty element is desired.

        For speed reasons (since unlike other p-adic types, this function can still be relatively fast compared to just _call_()), the empty keyword has no effect.

        EXAMPLES::

            sage: R = ZpFM(5,4)
            sage: type(R(10,2))
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R(30,2)
            5 + 5^2 + O(5^4)
            sage: R(30,3,1)
            5 + 5^2 + O(5^4)
            sage: R(30,absprec=2)
            5 + 5^2 + O(5^4)
            sage: R(30,relprec=2)
            5 + 5^2 + O(5^4)
            sage: R(30,absprec=1)
            5 + 5^2 + O(5^4)
            sage: R(30,empty=True)
            5 + 5^2 + O(5^4)
        """
        cdef pAdicFixedModElement ans = PY_NEW(pAdicFixedModElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        ans._set_from_mpz((<Integer>x).value)
        return ans

    def section(self):
        """
        Returns a map back to ZZ that approximates an element of Zp by an integer.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ).section()
            sage: f(ZpFM(5)(-1)) - 5^20
            -1
        """
        return self._section

cdef class pAdicConvert_FM_ZZ(RingMap):
    """
    The map from ZpFM back to ZZ that returns the the smallest non-negative integer approximation to its input which is accurate up to the precision.

    EXAMPLES::

        sage: f = ZpFM(5).coerce_map_from(ZZ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Ring of fixed modulus 5^20
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ).section(); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicConvert_FM_ZZ'>
            sage: f.category()
            Category of hom sets in Category of sets
        """
        RingMap.__init__(self, Hom(R, ZZ, Sets()))

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ).section()
            sage: f(ZpFM(5)(-1)) - 5^20
            -1
            sage: f(ZpFM(5)(0))
            0
        """
        cdef Integer ans = PY_NEW(Integer)
        (<pAdicCappedAbsoluteElement>x)._set_mpz_into(ans.value)
        return ans

cdef class pAdicConvert_QQ_FM(Morphism):
    """
    The inclusion map from QQ to ZpFM that is defined on all elements with non-negative p-adic valuation.

    EXAMPLES::

        sage: f = ZpFM(5).convert_map_from(QQ); f
        Generic morphism:
          From: Rational Field
          To:   5-adic Ring of fixed modulus 5^20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFM(5).convert_map_from(QQ); type(f)
            <type 'sage.rings.padics.padic_base_coercion.pAdicConvert_QQ_FM'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self.prime_pow = R.prime_pow

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).convert_map_from(QQ)
            sage: g = copy(f) # indirect doctest
            sage: g == f      # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True

        """
        _slots['prime_pow'] = self.prime_pow
        try:
            _slots['_section'] = self._section
        except AttributeError:
            pass
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).convert_map_from(QQ)
            sage: g = copy(f) # indirect doctest
            sage: g == f      # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True

        """
        self.prime_pow = _slots['prime_pow']
        try:
            self._section = _slots['_section']
        except KeyError:
            pass
        Morphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpFM(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: f(0)
            O(5^4)
        """
        cdef pAdicFixedModElement ans = PY_NEW(pAdicFixedModElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        ans._set_from_mpq((<Rational>x).value)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both), or an empty element is desired.

        For speed reasons (since unlike other p-adic types, this function can still be relatively fast compared to just _call_()), the empty keyword has no effect.

        EXAMPLES::

            sage: R = ZpFM(5,4)
            sage: type(R(1/7,2))
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R(1/7,2)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,3,1)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,absprec=2)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,relprec=2)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,absprec=1)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,empty=True)
            3 + 3*5 + 2*5^3 + O(5^4)
        """
        cdef pAdicFixedModElement ans = PY_NEW(pAdicFixedModElement)
        ans._parent = self.codomain()
        ans.prime_pow = self.prime_pow
        mpz_init(ans.value)
        ans._set_from_mpq((<Rational>x).value)
        return ans

