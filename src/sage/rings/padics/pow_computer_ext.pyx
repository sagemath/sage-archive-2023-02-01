"""
AUTHOR:
    -- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include "../../ext/python_list.pxi"
include "../../ext/python_dict.pxi"

import weakref
from sage.misc.misc import cputime
from sage.rings.infinity import infinity
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_factory
from sage.libs.ntl.ntl_ZZ_pContext import ZZ_pContext_factory
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX, ntl_ZZ_pX_Modulus
from sage.rings.integer cimport Integer

cdef class PowComputer_ext(PowComputer_class):
    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        """
        Constructs the storage for powers of prime as ZZ_c's.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "small") #indirect doctest
        """
        self._initialized = 0
        _sig_on
        self.small_powers = <ZZ_c *>sage_malloc(sizeof(ZZ_c) * (cache_limit + 1))
        _sig_off
        if self.small_powers == NULL:
            raise MemoryError, "out of memory allocating power storing"
        ZZ_construct(&self.top_power)

        cdef Py_ssize_t i
        cdef Integer x

        ZZ_construct(self.small_powers)
        ZZ_conv_from_int(self.small_powers[0], 1)

        if cache_limit > 0:
            ZZ_construct(&(self.small_powers[1]))
            mpz_to_ZZ(&(self.small_powers[1]), &prime.value)

        _sig_on
        for i from 2 <= i <= cache_limit:
            ZZ_construct(&(self.small_powers[i]))
            ZZ_mul(self.small_powers[i], self.small_powers[i-1], self.small_powers[1])
        mpz_to_ZZ(&self.top_power, &prime.value)
        ZZ_power(self.top_power, self.top_power, prec_cap)
        _sig_off
        mpz_init(self.temp_m)
        ZZ_construct(&self.temp_z)

    def __init__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        """
        Initializes prime, cache_limit, prec_cap, in_field.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "small") #indirect doctest
        """
        PowComputer_class.__init__(self, prime, cache_limit, prec_cap, in_field)

    def __dealloc__(self):
        """
        Frees allocated memory.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: del PC # indirect doctest
        """
        if (<PowComputer_class>self)._initialized:
            self.cleanup_ext()

    def __repr__(self):
        """
        Returns a string representation of self.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "small")
        sage: PC # indirect doctest
        PowComputer_ext for 5, with polynomial [9765620 0 1]
        """
        return "PowComputer_ext for %s, with polynomial %s"%(self.prime, self.polynomial())

    cdef void cleanup_ext(self):
        """
        Frees memory allocated in PowComputer_ext.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: del PC # indirect doctest
        """
        cdef Py_ssize_t i
        for i from 0 <= i <= self.cache_limit:
            ZZ_destruct(&(self.small_powers[i]))
        sage_free(self.small_powers)
        ZZ_destruct(&self.top_power)
        mpz_clear(self.temp_m)
        ZZ_destruct(&self.temp_z)

    cdef mpz_t* pow_mpz_t_tmp(self, unsigned long n):
        """
        Provides fast access to an mpz_t* pointing to self.prime^n.

        The location pointed to depends on the underlying representation.
        In no circumstances should you mpz_clear the result.
        The value pointed to may be an internal temporary variable for the class.
        In particular, you should not try to refer to the results of two
        pow_mpz_t_tmp calls at the same time, because the second
        call may overwrite the memory pointed to by the first.

        In the case of PowComputer_exts, the mpz_t pointed to will always
        be a temporary variable.

        See pow_mpz_t_tmp_demo for an example of this phenomenon.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "small")
        sage: PC.pow_mpz_t_tmp_test(4) #indirect doctest
        625
        """
        # READ THE DOCSTRING
        if n <= self.cache_limit:
            ZZ_to_mpz(&self.temp_m, &(self.small_powers[n]))
        elif n == self.prec_cap:
            ZZ_to_mpz(&self.temp_m, &self.top_power)
        else:
            mpz_pow_ui(self.temp_m, self.prime.value, n)
        return &self.temp_m

    cdef ZZ_c* pow_ZZ_tmp(self, unsigned long n):
        """
        Provides fast access to a ZZ_c* pointing to self.prime^n.

        The location pointed to depends on the underlying representation.
        In no circumstances should you ZZ_destruct the result.
        The value pointed to may be an internal temporary variable for the class.
        In particular, you should not try to refer to the results of two
        pow_ZZ_tmp calls at the same time, because the second
        call may overwrite the memory pointed to by the first.

        See pow_ZZ_tmp_demo for an example of this phenomenon.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "small")
        sage: PC.pow_mpz_t_tmp_test(4) #indirect doctest
        625
        """
        if n <= self.cache_limit:
            return &(self.small_powers[n])
        if n == self.prec_cap:
            return &self.top_power
        ZZ_power(self.temp_z, self.small_powers[1], n)
        return &self.temp_z

    def pow_ZZ_tmp_test(self, n):
        """
        Tests the pow_ZZ_tmp function

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 6, 6, False, ntl.ZZ_pX([-5,0,1],5^6),"small")
        sage: PC.pow_ZZ_tmp_test(4)
        625
        sage: PC.pow_ZZ_tmp_test(7)
        78125
        """
        cdef Integer _n = Integer(n)
        if _n < 0: raise ValueError
        cdef ntl_ZZ ans = PY_NEW(ntl_ZZ)
        ans.x = self.pow_ZZ_tmp(mpz_get_ui(_n.value))[0]
        return ans

    def pow_ZZ_tmp_demo(self, m, n):
        """
        This function demonstrates a danger in using pow_ZZ_tmp.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")

        When you cal pow_ZZ_tmp with an input that is not stored
        (ie n > self.cache_limit and n != self.prec_cap),
        it stores the result in self.temp_z and returns a pointer
        to that ZZ_c.  So if you try to use the results of two
        calls at once, things will break.
        sage: PC.pow_ZZ_tmp_demo(6, 8)
        244140625
        sage: 5^6*5^8
        6103515625
        sage: 5^6*5^6
        244140625

        Note that this does not occur if you try a stored value,
        because the result of one of the calls points to that
        stored value.
        sage: PC.pow_ZZ_tmp_demo(6, 10)
        152587890625
        sage: 5^6*5^10
        152587890625
        """
        m = Integer(m)
        n = Integer(n)
        if m < 0 or n < 0:
            raise ValueError, "m, n must be non-negative"
        cdef ntl_ZZ ans = PY_NEW(ntl_ZZ)
        ZZ_mul(ans.x, self.pow_ZZ_tmp(mpz_get_ui((<Integer>m).value))[0], self.pow_ZZ_tmp(mpz_get_ui((<Integer>n).value))[0])
        return ans


    cdef mpz_t* pow_mpz_t_top(self):
        """
        Returns self.prime^self.prec_cap as an mpz_t*.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 6, 6, False, ntl.ZZ_pX([-5,0,1],5^6),"small")
        sage: PC.pow_mpz_t_top_test() #indirect doctest
        15625
        """
        ZZ_to_mpz(&self.temp_m, &self.top_power)
        return &self.temp_m

    cdef ZZ_c* pow_ZZ_top(self):
        """
        Returns self.prime^self.prec_cap as a ZZ_c.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 6, 6, False, ntl.ZZ_pX([-5,0,1],5^6),"small")
        sage: PC.pow_ZZ_top_test() #indirect doctest
        15625
        """
        return &self.top_power

    def pow_ZZ_top_test(self):
        """
        Tests the pow_ZZ_top function.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 6, 6, False, ntl.ZZ_pX([-5,0,1],5^6),"small")
        sage: PC.pow_ZZ_top_test()
        15625
        """
        cdef ntl_ZZ ans = PY_NEW(ntl_ZZ)
        ans.x = self.pow_ZZ_top()[0]
        return ans

cdef class PowComputer_ZZ_pX(PowComputer_ext):
    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        if not PY_TYPE_CHECK(poly, ntl_ZZ_pX):
            raise TypeError
        self.deg = ZZ_pX_deg((<ntl_ZZ_pX>poly).x)

    def polynomial(self):
        """
        Returns the polynomial (with coefficient precision prec_cap) associated to this PowComputer.

        The polynomial is output as an ntl_ZZ_pX.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.polynomial()
        [9765620 0 1]
        """
        self.restore_top_context()
        cdef ntl_ZZ_pX r = PY_NEW(ntl_ZZ_pX)
        r.c = self.get_top_context()
        r.x = (self.get_top_modulus()[0]).val()
        return r

    cdef ntl_ZZ_pContext_class get_context(self, unsigned long n):
        """
        Returns a ZZ_pContext for self.prime^n.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "FM")
        sage: PC.get_context_test(15) #indirect doctest
        NTL modulus 30517578125
        """
        cdef ntl_ZZ pn = PY_NEW(ntl_ZZ)
        pn.x = self.pow_ZZ_tmp(n)[0]
        cdef ntl_ZZ_pContext_class context = (<ntl_ZZ_pContext_factory>ZZ_pContext_factory).make_c(pn)
        return context

    def get_context_test(self, n):
        """
        Returns a ZZ_pContext for self.prime^n.

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "FM")
        sage: PC.get_context_test(15)
        NTL modulus 30517578125
        """
        cdef Integer _n = Integer(n)
        if _n < 1: raise ValueError
        return self.get_context(mpz_get_ui(_n.value))

    def speed_test(self, n, runs):
        """
        Runs a speed test.

        INPUT:
        n -- input to a function to be tested (the function needs to be set in the source code).
        runs -- The number of runs of that function
        OUTPUT:
        The time in seconds that it takes to call the function on n, runs times.

        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "small")
        sage: PC.speed_test(10, 10^6) # random
        0.0090679999999991878
        """
        cdef Py_ssize_t i, end, _n
        end = mpz_get_ui((<Integer>Integer(runs)).value)
        _n = mpz_get_ui((<Integer>Integer(n)).value)
        t = cputime()
        for i from 0 <= i < end:
            # Put the function you want speed tested here.
            self.get_modulus(_n)
        return cputime(t)

    cdef ntl_ZZ_pContext_class get_top_context(self):
        """
        Returns a ZZ_pContext for self.prime^self.prec_cap

        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.get_top_context_test() #indirect doctest
        NTL modulus 9765625
        """
        return self.get_context(self.prec_cap)

    def get_top_context_test(self):
        """
        Returns a ZZ_pContext for self.prime^self.prec_cap

        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.get_top_context_test()
        NTL modulus 9765625
        """
        return self.get_top_context()

    cdef restore_context(self, unsigned long n):
        """
        Restores the contest corresponding to self.prime^n

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.restore_context_test(4) #indirect doctest
        """
        self.get_context(n).restore_c()

    def restore_context_test(self, n):
        """
        Restores the contest corresponding to self.prime^n

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.restore_context_test(4)
        """
        cdef Integer _n = Integer(n)
        if _n < 0: raise ValueError
        self.restore_context(mpz_get_ui(_n.value))

    cdef void restore_top_context(self):
        """
        Restores the context corresponding to self.prime^self.prec_cap

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.restore_top_context_test()
        """
        (<ntl_ZZ_pContext_class>self.get_top_context()).restore_c()

    def restore_top_context_test(self):
        """
        Restores the context corresponding to self.prime^self.prec_cap

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.restore_top_context_test()
        """
        self.restore_top_context()

    cdef ZZ_pX_Modulus_c* get_modulus(self, unsigned long n):
        """
        Returns the modulus corresponding to self.polynomial() (mod self.prime^n)

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 1000, False, ntl.ZZ_pX([-5,0,1],5^1000), "big")
        sage: a = ntl.ZZ_pX([4,2],5^2)
        sage: b = ntl.ZZ_pX([6,3],5^2)
        sage: A.get_modulus_test(a, b, 2) # indirect doctest
        [4 24]
        """
        raise NotImplementedError

    def get_modulus_test(self, ntl_ZZ_pX a, ntl_ZZ_pX b, Integer n):
        """
        Multiplies a and b modulo the modulus corresponding to self.polynomial() (mod self.prime^n).

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 1000, False, ntl.ZZ_pX([-5,0,1],5^1000), "big")
        sage: a = ntl.ZZ_pX([4,2],5^2)
        sage: b = ntl.ZZ_pX([6,3],5^2)
        sage: A.get_modulus_test(a, b, 2)
        [4 24]
        sage: a * b
        [24 24 6]
        sage: mod(6 * 5 + 24, 25)
        4
        """
        cdef ntl_ZZ_pX r = (<ntl_ZZ_pX>a)._new()
        ZZ_pX_MulMod_pre(r.x, a.x, b.x, self.get_modulus(mpz_get_ui(n.value))[0])
        return r

    cdef ZZ_pX_Modulus_c* get_top_modulus(self):
        """
        Returns the modulus corresponding to self.polynomial() (mod self.prime^self.prec_cap)

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: a = ntl.ZZ_pX([129223,1231],5^10)
        sage: b = ntl.ZZ_pX([289741,323],5^10)
        sage: A.get_top_modulus_test(a, b) #indirect doctest
        [1783058 7785200]
        """
        raise NotImplementedError

    def get_top_modulus_test(self, ntl_ZZ_pX a, ntl_ZZ_pX b):
        """
        Multiplies a and b modulo the modulus corresponding to self.polynomial() (mod self.prime^self.prec_cap)

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: a = ntl.ZZ_pX([129223,1231],5^10)
        sage: b = ntl.ZZ_pX([289741,323],5^10)
        sage: A.get_top_modulus_test(a, b)
        [1783058 7785200]
        sage: a*b
        [9560618 7785200 397613]
        sage: mod(397613 * 5 + 9560618, 5^10)
        1783058
        """
        cdef ntl_ZZ_pX ans = a._new()
        ZZ_pX_MulMod_pre(ans.x, a.x, b.x, self.get_top_modulus()[0])
        return ans

cdef class PowComputer_ZZ_pX_FM(PowComputer_ZZ_pX):
    """
    This class only caches a context and modulus for p^prec_cap.
    Designed for use with fixed modulus p-adic rings, in Eisenstein and unramified extensions of $\mathbb{Z}_p$.
    """

    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        """
        Caches a context and modulus for prime^prec_cap

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM") #indirect doctest
        sage: A
        PowComputer_ext for 5, with polynomial [9765620 0 1]
        """

        # The __new__ method for PowComputer_ext has already run, so we have access to small_powers, top_power.

        # We use ntl_ZZ_pContexts so that contexts are cached centrally.

        self.c = self.get_context(prec_cap)
        self.c.restore_c()
        # For now, we don't do anything complicated with poly
        if PY_TYPE_CHECK(poly, ntl_ZZ_pX) and (<ntl_ZZ_pX>poly).c is self.c:
            ZZ_pX_Modulus_construct(&self.mod)
            ZZ_pX_Modulus_build(self.mod, (<ntl_ZZ_pX>poly).x)
            # These will be reset if we're actually eisenstein
            self.e = 1
            self.f = ZZ_pX_deg((<ntl_ZZ_pX>poly).x)
            self.ram_prec_cap = self.prec_cap
        else:
            print "NOT IMPLEMENTED IN PowComputer_ZZ_pX_FM"
            raise NotImplementedError

    #def __init__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
    #    """
    #    Caches a context and modulus for prime^prec_cap
    #
    #    EXAMPLES:
    #    sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM") #indirect doctest
    #    sage: A
    #    PowComputer_ext for 5, with polynomial [9765620 0 1]
    #    """
    #    PowComputer_ZZ_pX.__init__(self, prime, cache_limit, prec_cap, in_field, poly)

    def __dealloc__(self):
        """
        Cleans up the memory for self.mod

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM") #indirect doctest
        sage: del A # indirect doctest
        """
        if self._initialized:
            self.cleanup_ZZ_pX_FM()

    cdef void cleanup_ZZ_pX_FM(self):
        """
        Cleans up the memory for self.mod

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM") #indirect doctest
        sage: del A # indirect doctest
        """
        ZZ_pX_Modulus_destruct(&self.mod)

    cdef ntl_ZZ_pContext_class get_top_context(self):
        """
        Returns a ZZ_pContext for self.prime^self.prec_cap

        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.get_top_context_test() # indirect doctest
        NTL modulus 9765625
        """
        return self.c

    cdef void restore_top_context(self):
        """
        Restores the context corresponding to self.prime^self.prec_cap

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: PC.restore_top_context_test() #indirect doctest
        """
        self.c.restore_c()

    cdef ZZ_pX_Modulus_c* get_top_modulus(self):
        """
        Returns the modulus corresponding to self.polynomial() (mod self.prime^self.prec_cap)

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "FM")
        sage: a = ntl.ZZ_pX([129223,1231],5^10)
        sage: b = ntl.ZZ_pX([289741,323],5^10)
        sage: A.get_top_modulus_test(a, b) #indirect doctest
        [1783058 7785200]
        """
        return &self.mod

cdef class PowComputer_ZZ_pX_FM_Eis(PowComputer_ZZ_pX_FM):
    """
    This class computes and stores low_shifter and high_shifter, which aid in right shifting elements.
    """

    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        # The __new__ method for PowComputer_ZZ_pX_FM has already run, so we have access to self.mod

        # self.low_shifter stores multipliers for p/x^(2^i)
        # If self.deg is one more than a power of 2, we need to store p/x, p/x^2, up to p/x^(2^n) where 2^n is that power of 2.
        #print "here"
        #cdef ntl_ZZ_pX printer = ntl_ZZ_pX([], self.c)
        ZZ_pX_Eis_init(self, self.low_shifter, self.high_shifter)

    def _low_shifter(self, i):
        cdef long _i = i
        cdef ntl_ZZ_pX ans
        if _i >= 0 and _i < self.low_length:
            ans = ntl_ZZ_pX([], self.get_top_context())
            ans.x = self.low_shifter[i].val()
            return ans
        else:
            raise IndexError

    def _high_shifter(self, i):
        cdef long _i = i
        cdef ntl_ZZ_pX ans
        if _i >= 0 and _i < self.high_length:
            ans = ntl_ZZ_pX([], self.get_top_context())
            ans.x = self.high_shifter[i].val()
            return ans
        else:
            raise IndexError

    def __dealloc__(self):
        if self._initialized:
            self.cleanup_ZZ_pX_FM_Eis()

    cdef void cleanup_ZZ_pX_FM_Eis(self):
        pass
        # I may or may not need to deallocate these:
        #
        #cdef int i # yes, an int is good enough
        #for i from 0 <= i < self.low_length:
        #    ZZ_pX_Multiplier_destruct(self.low_shifter[i])
        #sage_free(self.low_shifter)
        #for i from 0 <= i < self.high_length:
        #    ZZ_pX_Multiplier_destruct(self.high_shifter[i])
        #sage_free(self.high_shifter)


cdef class PowComputer_ZZ_pX_small(PowComputer_ZZ_pX):
    """
    This class caches contexts and moduli densely between 1 and cache_limit.  It requires cache_limit == prec_cap.

    It is intended for use with capped relative and capped absolute rings and fields, in Eisenstein and unramified
    extensions of the base p-adic fields.
    """

    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        """
        Caches contexts and moduli densely between 1 and cache_limit.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small") # indirect doctest
        sage: A
        PowComputer_ext for 5, with polynomial [9765620 0 1]
        """
        # The __new__ method for PowComputer_ext has already run, so we have access to small_powers, top_power.

        # We use ntl_ZZ_pContexts so that contexts are cached centrally.

        if not PY_TYPE_CHECK(poly, ntl_ZZ_pX):
            self.cleanup_ext()
            self.cleanup_ZZ_pX()
            raise TypeError

        if cache_limit != prec_cap:
            self.cleanup_ext()
            self.cleanup_ZZ_pX()
            raise ValueError, "prec_cap and cache_limit must be equal in the small case"

        self.c = []
        #if self.c == NULL:
        #    self.cleanup_ext()
        #    self.cleanup_ZZ_pX()
        #    raise MemoryError, "out of memory allocating contexts"
        _sig_on
        self.mod = <ZZ_pX_Modulus_c *>sage_malloc(sizeof(ZZ_pX_Modulus_c) * (cache_limit + 1))
        _sig_off
        if self.mod == NULL:
            self.cleanup_ext()
            self.cleanup_ZZ_pX()
            raise MemoryError, "out of memory allocating moduli"

        cdef Py_ssize_t i
        cdef ZZ_pX_c tmp, pol
        ZZ_pX_construct(&tmp)
        ZZ_pX_construct(&pol)
        pol = (<ntl_ZZ_pX>poly).x
        self.c.append(None)
        for i from 1 <= i <= cache_limit:
            self.c.append(PowComputer_ZZ_pX.get_context(self,i))
            ZZ_pX_Modulus_construct(&(self.mod[i]))
            ZZ_pX_conv_modulus(tmp, pol, (<ntl_ZZ_pContext_class>self.c[i]).x)
            ZZ_pX_Modulus_build(self.mod[i], tmp)
        ## I get malloc errors if I leave the following in: I don't know why.
        #ZZ_pX_destruct(&tmp)
        #ZZ_pX_destruct(&pol)

    #def __init__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
    #    """
    #    Initializes prime, cache_limit, prec_cap, in_field.
    #
    #    EXAMPLES:
    #    sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small") # indirect doctest
    #    sage: A
    #    PowComputer_ext for 5, with polynomial [9765620 0 1]
    #    """
    #    PowComputer_ZZ_pX.__init__(self, prime, cache_limit, prec_cap, in_field, poly)

    def __dealloc__(self):
        """
        Deallocates cache of contexts, moduli.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: del A # indirect doctest
        """
        if self._initialized:
            self.cleanup_ZZ_pX_small()

    cdef void cleanup_ZZ_pX_small(self):
        """
        Deallocates cache of contexts, moduli.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: del A # indirect doctest
        """
        pass
        ## These cause segfaults: I don't know why.
        #cdef Py_ssize_t i
        #for i from 0 <= i <= self.cache_limit:
        #    ZZ_pX_Modulus_destruct(&(self.mod[i]))
        #sage_free(self.mod)

    cdef ntl_ZZ_pContext_class get_context(self, unsigned long n):
        """
        Returns the context for p^n.

        Note that this function will raise an Index error if n > self.cache_limit.
        Also, it will return None on input 0

        INPUT:
        n -- A long between 1 and self.cache_limit, inclusive
        OUTPUT:
        A context for p^n

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: A.get_context_test(4) #indirect doctest
        NTL modulus 625
        """
        return self.c[n]

    cdef restore_context(self, unsigned long n):
        """
        Restores the context for p^n.

        INPUT:
        n -- A long between 1 and self.cache_limit, inclusive

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: A.restore_context_test(4) #indirect doctest
        """
        _sig_on
        (<ntl_ZZ_pContext_class>self.c[n]).restore_c()
        _sig_off

    cdef ntl_ZZ_pContext_class get_top_context(self):
        """
        Returns a ZZ_pContext for self.prime^self.prec_cap

        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: PC.get_top_context_test() # indirect doctest
        NTL modulus 9765625
        """
        return self.c[self.prec_cap]

    cdef void restore_top_context(self):
        """
        Restores the context corresponding to self.prime^self.prec_cap

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: PC.restore_top_context_test() #indirect doctest
        """
        (<ntl_ZZ_pContext_class>self.c[self.prec_cap]).restore_c()

    cdef ZZ_pX_Modulus_c* get_modulus(self, unsigned long n):
        """
        Returns the modulus corresponding to self.polynomial() (mod self.prime^n).

        INPUT:
        n -- A long between 1 and self.cache_limit, inclusive

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: a = ntl.ZZ_pX([4,2],5^2)
        sage: b = ntl.ZZ_pX([6,3],5^2)
        sage: A.get_modulus_test(a, b, 2)
        [4 24]
        """
        return &(self.mod[n])

    cdef ZZ_pX_Modulus_c* get_top_modulus(self):
        """
        Returns the modulus corresponding to self.polynomial() (mod self.prime^self.prec_cap)

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "small")
        sage: a = ntl.ZZ_pX([129223,1231],5^10)
        sage: b = ntl.ZZ_pX([289741,323],5^10)
        sage: A.get_top_modulus_test(a, b) #indirect doctest
        [1783058 7785200]
        """
        return &(self.mod[self.prec_cap])

cdef class PowComputer_ZZ_pX_small_Eis(PowComputer_ZZ_pX_small):
    """
    This class computes and stores low_shifter and high_shifter, which aid in right shifting elements.
    These are only stored at maximal precision: in order to get lower precision versions just reduce mod p^n.
    """
    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        ZZ_pX_Eis_init(self, self.low_shifter, self.high_shifter)

    def _low_shifter(self, i):
        cdef long _i = i
        cdef ntl_ZZ_pX ans
        if _i >= 0 and _i < self.low_length:
            ans = ntl_ZZ_pX([], self.get_top_context())
            ans.x = self.low_shifter[i].val()
            return ans
        else:
            raise IndexError

    def _high_shifter(self, i):
        cdef long _i = i
        cdef ntl_ZZ_pX ans
        if _i >= 0 and _i < self.high_length:
            ans = ntl_ZZ_pX([], self.get_top_context())
            ans.x = self.high_shifter[i].val()
            return ans
        else:
            raise IndexError


    def __dealloc__(self):
        if self._initialized:
            self.cleanup_ZZ_pX_small_Eis()

    cdef void cleanup_ZZ_pX_small_Eis(self):
        pass
        # I may or may not need to deallocate these:
        #
        #cdef int i # yes, an int is good enough
        #for i from 0 <= i < self.low_length:
        #    ZZ_pX_Multiplier_destruct(self.low_shifter[i])
        #sage_free(self.low_shifter)
        #for i from 0 <= i < self.high_length:
        #    ZZ_pX_Multiplier_destruct(self.high_shifter[i])
        #sage_free(self.high_shifter)

cdef class PowComputer_ZZ_pX_big(PowComputer_ZZ_pX):
    """
    This class caches all contexts and moduli between 1 and cache_limit, and also caches for prec_cap.  In addition, it stores
    a dictionary of contexts and moduli of
    """

    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        """
        Caches contexts and moduli densely between 1 and cache_limit.  Caches a context and modulus for prec_cap.
        Also creates the dictionaries.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big") # indirect doctest
        sage: A
        PowComputer_ext for 5, with polynomial [9765620 0 1]
        """
        # The __new__ method for PowComputer_ext has already run, so we have access to small_powers, top_power.

        # We use ntl_ZZ_pContexts so that contexts are cached centrally.

        if not PY_TYPE_CHECK(poly, ntl_ZZ_pX):
            self.cleanup_ext()
            self.cleanup_ZZ_pX()
            raise TypeError

        self.context_list = []
        #if self.c == NULL:
        #    self.cleanup_ext()
        #    self.cleanup_ZZ_pX()
        #    raise MemoryError, "out of memory allocating contexts"
        _sig_on
        self.modulus_list = <ZZ_pX_Modulus_c *>sage_malloc(sizeof(ZZ_pX_Modulus_c) * (cache_limit + 1))
        _sig_off
        if self.modulus_list == NULL:
            self.cleanup_ext()
            self.cleanup_ZZ_pX()
            raise MemoryError, "out of memory allocating moduli"

        cdef Py_ssize_t i
        cdef ZZ_pX_c tmp, pol
        ZZ_pX_construct(&tmp)
        ZZ_pX_construct(&pol)
        pol = (<ntl_ZZ_pX>poly).x
        self.context_list.append(None)
        for i from 1 <= i <= cache_limit:
            self.context_list.append(PowComputer_ZZ_pX.get_context(self,i))
            ZZ_pX_Modulus_construct(&(self.modulus_list[i]))
            ZZ_pX_conv_modulus(tmp, pol, (<ntl_ZZ_pContext_class>self.context_list[i]).x)
            ZZ_pX_Modulus_build(self.modulus_list[i], tmp)
        self.top_context = PowComputer_ZZ_pX.get_context(self, prec_cap)
        ZZ_pX_Modulus_construct(&(self.top_mod))
        ZZ_pX_conv_modulus(tmp, pol, self.top_context.x)
        ZZ_pX_Modulus_build(self.top_mod, tmp)
        ## I get malloc errors if I leave the following in: I don't know why.
        #ZZ_pX_destruct(&tmp)
        #ZZ_pX_destruct(&pol)
        self.context_dict = {}
        self.modulus_dict = {}

    #def __init__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
    #    """
    #    Initializes prime, cache_limit, prec_cap, in_field.
    #
    #    EXAMPLES:
    #    sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big") # indirect doctest
    #    sage: A
    #    PowComputer_ext for 5, with polynomial [9765620 0 1]
    #    """
    #    PowComputer_ZZ_pX.__init__(self, prime, cache_limit, prec_cap, in_field, poly)

    def __dealloc__(self):
        """
        Deallocates the stored moduli and contexts.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: del A # indirect doctest
        """
        if self._initialized:
            self.cleanup_ZZ_pX_big()

    cdef void cleanup_ZZ_pX_big(self):
        """
        Deallocates the stored moduli and contexts.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: del A # indirect doctest
        """
        pass
        ## These cause a segfault.  I don't know why.
        #cdef Py_ssize_t i
        #for i from 0 <= i <= self.cache_limit:
        #    ZZ_pX_Modulus_destruct(&(self.modulus_list[i]))
        #sage_free(self.modulus_list)
        #ZZ_pX_Modulus_destruct(&self.top_mod)

    def reset_dictionaries(self):
        """
        Resets the dictionaries.  Note that if there are elements lying around that need access to these dictionaries, calling this function and then doing arithmetic with those elements could cause trouble (if the context object gets garbage collected for example.  The bugs introduced could be very subtle, because NTL will generate a new context object and use it, but there's the potential for the object to be incompatible with the different context object).

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: P = A.get_context_test(8)
        sage: A._context_dict()
        {8L: NTL modulus 390625}
        sage: A.reset_dictionaries()
        """
        self.context_dict = {}
        self.modulus_dict = {}

    def _context_dict(self):
        """
        Returns the context dictionary.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: P = A.get_context_test(8)
        sage: A._context_dict()
        {8L: NTL modulus 390625}
        """
        return self.context_dict

    def _modulus_dict(self):
        """
        Returns the context dictionary.

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: P = A.get_context_test(8)
        sage: A._modulus_dict()
        {}
        sage: a = ntl.ZZ_pX([4,2],5^8)
        sage: b = ntl.ZZ_pX([6,3],5^8)
        sage: A.get_modulus_test(a, b, 8)
        [54 24]
        sage: A._modulus_dict()
        {8L: NTL ZZ_pXModulus [390620 0 1] (mod 390625)}
        """
        return self.modulus_dict

    cdef ntl_ZZ_pContext_class get_context(self, unsigned long n):
        """
        Returns the context for p^n.

        Note that this function will raise an Index error if n > self.cache_limit.
        Also, it will return None on input 0

        INPUT:
        n -- A long between 1 and self.cache_limit, inclusive
        OUTPUT:
        A context for p^n

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 6, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: A.get_context_test(4) #indirect doctest
        NTL modulus 625
        sage: A.get_context_test(8) #indirect doctest
        NTL modulus 390625
        """
        if n == 0:
            raise ValueError, "n must be positive"
        if n <= self.cache_limit:
            return self.context_list[n]
        elif n == self.prec_cap:
            return self.top_context
        else:
            try:
                return self.context_dict[n]
            except KeyError:
                self.context_dict[n] = PowComputer_ZZ_pX.get_context(self, n)
                return self.context_dict[n]

    cdef ntl_ZZ_pContext_class get_top_context(self):
        """
        Returns a ZZ_pContext for self.prime^self.prec_cap

        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: PC.get_top_context_test() # indirect doctest
        NTL modulus 9765625
        """
        return self.top_context

    cdef void restore_top_context(self):
        """
        Restores the context corresponding to self.prime^self.prec_cap

        EXAMPLES:
        sage: PC = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: PC.restore_top_context_test() #indirect doctest
        """
        self.top_context.restore_c()

    cdef ZZ_pX_Modulus_c* get_modulus(self, unsigned long n):
        """
        Returns the modulus corresponding to self.polynomial() (mod self.prime^n).

        INPUT:
        n -- A long between 1 and self.cache_limit, inclusive

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 3, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: a = ntl.ZZ_pX([4,2],5^2)
        sage: b = ntl.ZZ_pX([6,3],5^2)
        sage: A.get_modulus_test(a, b, 2) # indirect doctest
        [4 24]
        sage: a = ntl.ZZ_pX([4,2],5^6)
        sage: b = ntl.ZZ_pX([6,3],5^6)
        sage: A.get_modulus_test(a, b, 6) # indirect doctest
        [54 24]
        """
        cdef ntl_ZZ_pX tmp
        cdef ntl_ZZ_pX_Modulus holder
        cdef ntl_ZZ_pContext_class c
        if n == 0:
            raise ValueError, "n must be positive"
        elif n <= self.cache_limit:
            return &(self.modulus_list[n])
        elif n == self.prec_cap:
            return &self.top_mod
        else:
            try:
                holder = self.modulus_dict[n]
                return &(holder.x)
            except KeyError:
                c = self.get_context(n)
                c.restore_c()
                tmp = PY_NEW(ntl_ZZ_pX)
                tmp.c = c
                ZZ_pX_conv_modulus(tmp.x, self.top_mod.val(), c.x)
                holder = ntl_ZZ_pX_Modulus(tmp)
                self.modulus_dict[n] = holder
                return &(holder.x)

    cdef ZZ_pX_Modulus_c* get_top_modulus(self):
        """
        Returns the modulus corresponding to self.polynomial() (mod self.prime^self.prec_cap)

        EXAMPLES:
        sage: A = PowComputer_ext_maker(5, 5, 10, False, ntl.ZZ_pX([-5,0,1],5^10), "big")
        sage: a = ntl.ZZ_pX([129223,1231],5^10)
        sage: b = ntl.ZZ_pX([289741,323],5^10)
        sage: A.get_top_modulus_test(a, b) #indirect doctest
        [1783058 7785200]
        """
        return &self.top_mod

cdef class PowComputer_ZZ_pX_big_Eis(PowComputer_ZZ_pX_big):
    """
    This class computes and stores low_shifter and high_shifter, which aid in right shifting elements.
    These are only stored at maximal precision: in order to get lower precision versions just reduce mod p^n.
    """
    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        ZZ_pX_Eis_init(self, self.low_shifter, self.high_shifter)

    def _low_shifter(self, i):
        cdef long _i = i
        cdef ntl_ZZ_pX ans
        if _i >= 0 and _i < self.low_length:
            ans = ntl_ZZ_pX([], self.get_top_context())
            ans.x = self.low_shifter[i].val()
            return ans
        else:
            raise IndexError

    def _high_shifter(self, i):
        cdef long _i = i
        cdef ntl_ZZ_pX ans
        if _i >= 0 and _i < self.high_length:
            ans = ntl_ZZ_pX([], self.get_top_context())
            ans.x = self.high_shifter[i].val()
            return ans
        else:
            raise IndexError


    def __dealloc__(self):
        if self._initialized:
            self.cleanup_ZZ_pX_big_Eis()

    cdef void cleanup_ZZ_pX_big_Eis(self):
        pass
        # I may or may not need to deallocate these:
        #
        #cdef int i # yes, an int is good enough
        #for i from 0 <= i < self.low_length:
        #    ZZ_pX_Multiplier_destruct(self.low_shifter[i])
        #sage_free(self.low_shifter)
        #for i from 0 <= i < self.high_length:
        #    ZZ_pX_Multiplier_destruct(self.high_shifter[i])
        #sage_free(self.high_shifter)


def PowComputer_ext_maker(prime, cache_limit, prec_cap, in_field, poly, prec_type = "small", ext_type = "u"):
    """
    Returns a PowComputer that caches the values $1, prime, prime^2, \ldots, prime^cache_limit$.

    Once you create a PowComputer, merely call it to get values out.
    You can input any integer, even if it's outside of the precomputed range.

    INPUT:
    prime -- An integer, the base that you want to exponentiate.
    cache_limit -- A positive integer that you want to cache powers up to.

    EXAMPLES:
    sage: PC = PowComputer_ext_maker(5, 10, 10, False, ntl.ZZ_pX([-5, 0, 1], 5^10), "small")
    sage: PC
    PowComputer_ext for 5, with polynomial [9765620 0 1]
    """
    cdef Integer _prime = Integer(prime)
    cdef Integer _cache_limit = Integer(cache_limit)
    cdef Integer _prec_cap = Integer(prec_cap)
    cdef bint inf = in_field
    if prec_type == "small" and ext_type == "u":
        return PowComputer_ZZ_pX_small(_prime, mpz_get_ui(_cache_limit.value), mpz_get_ui(_prec_cap.value), inf, poly)
    elif prec_type == "small" and ext_type == "e":
        return PowComputer_ZZ_pX_small_Eis(_prime, mpz_get_ui(_cache_limit.value), mpz_get_ui(_prec_cap.value), inf, poly)
    elif prec_type == "big" and ext_type == "u":
        return PowComputer_ZZ_pX_big(_prime, mpz_get_ui(_cache_limit.value), mpz_get_ui(_prec_cap.value), inf, poly)
    elif prec_type == "big" and ext_type == "e":
        return PowComputer_ZZ_pX_big_Eis(_prime, mpz_get_ui(_cache_limit.value), mpz_get_ui(_prec_cap.value), inf, poly)
    elif prec_type == "FM" and ext_type == "u":
        return PowComputer_ZZ_pX_FM(_prime, mpz_get_ui(_cache_limit.value), mpz_get_ui(_prec_cap.value), inf, poly)
    elif prec_type == "FM" and ext_type == "e":
        return PowComputer_ZZ_pX_FM_Eis(_prime, mpz_get_ui(_cache_limit.value), mpz_get_ui(_prec_cap.value), inf, poly)
    else:
        ValueError, "prec_type must be one of 'small', 'big' or 'FM' and ext_type must be one of 'u' or 'e'"


cdef int ZZ_pX_Eis_init(PowComputer_ZZ_pX prime_pow, ZZ_pX_Multiplier_c* low_shifter, ZZ_pX_Multiplier_c* high_shifter) except -1:
    if prime_pow.deg <= 1:
        raise ValueError, "Eisenstein extension must have degree at least 2"
    prime_pow.e = ZZ_pX_deg(prime_pow.get_top_modulus()[0].val())
    prime_pow.f = 1
    prime_pow.ram_prec_cap = prime_pow.prec_cap * prime_pow.e
    cdef unsigned long D = prime_pow.deg - 1
    cdef int low_length = 0
    cdef int high_length = 0
    if sizeof(long) > 4 and D > 4294967295: # 2^32 - 1
        low_length += 32
        D = D >> 32
    if D >= 65536: # 2^16
        low_length += 16
        D = D >> 16
    if D >= 256: # 2^8
        low_length += 8
        D = D >> 8
    if D >= 16: # 2^4
        low_length += 4
        D = D >> 4
    if D >= 4: # 2^2
        low_length += 2
        D = D >> 2
    if D >= 2: # 2^1
        low_length += 1
        D = D >> 1
    low_length += 1
    # low_length is the number of elements in the list we need to store.
    # if deg = 2, low_length = 1 (store p/x)
    # if deg = 3,4, low_length = 2 (store p/x, p/x^2)
    # if deg = 5,6,7,8, low_length = 3 (store p/x, p/x^2, p/x^4)
    # if deg = 9,...,16, low_length = 4 (store p/x, p/x^2, p/x^4, p/x^8)

    # Now we do the same process for powers of p, ie storing p^(2^k)/x^(e*2^k)
    D = prime_pow.prec_cap - 1
    high_length = 0
    if sizeof(long) > 4 and D > 4294967295: # 2^32 - 1
        high_length += 32
        D = D >> 32
    if D >= 65536: # 2^16
        high_length += 16
        D = D >> 16
    if D >= 256: # 2^8
        high_length += 8
        D = D >> 8
    if D >= 16: # 2^4
        high_length += 4
        D = D >> 4
    if D >= 4: # 2^2
        high_length += 2
        D = D >> 2
    if D >= 2: # 2^1
        high_length += 1
        D = D >> 1
    high_length += 1
    # high_length is the number of elements in the list we need to store.
    # if prec_cap = 2, high_length = 1 (store p/x^e)
    # if prec_cap = 3,4, high_length = 2 (store p/x^e, p^2/x^(2e))
    # if prec_cap = 5,6,7,8, high_length = 3 (store p/x^e, p^2/x^(2e), p^4/x^(4e))
    # if prec_cap = 9,...,16, high_length = 4 (store p/x, p^2/x^(2e), p^4/x^(4e), p^8/x^(8e))

    prime_pow.low_length = low_length
    prime_pow.high_length = high_length

    _sig_on
    low_shifter = <ZZ_pX_Multiplier_c *>sage_malloc(sizeof(ZZ_pX_Multiplier_c) * low_length)
    high_shifter = <ZZ_pX_Multiplier_c *>sage_malloc(sizeof(ZZ_pX_Multiplier_c) * high_length)
    _sig_off
    cdef long i
    cdef ZZ_pX_c tmp, modup, into_multiplier
    cdef ZZ_c a
    ZZ_construct(&a)
    # We obtain successive p/x^(2^i) by squaring and then dividing by p.  So we need one extra digit of precision.
    prime_pow.restore_top_context()
    ZZ_pX_construct(&into_multiplier)
    cdef ntl_ZZ_pContext_class cup = prime_pow.get_context(prime_pow.prec_cap + prime_pow.low_length)
    cup.restore_c()
    ZZ_pX_construct(&tmp)
    ZZ_pX_construct(&modup)
    ZZ_pX_conv_modulus(modup, prime_pow.get_top_modulus()[0].val(), cup.x)
    ZZ_div(a, ZZ_p_rep(ZZ_pX_ConstTerm(modup)), prime_pow.small_powers[1])
    ZZ_InvMod(a, a, prime_pow.pow_ZZ_tmp(prime_pow.prec_cap + prime_pow.low_length)[0])
    ZZ_negate(a, a)
    #cdef ntl_ZZ_pX printer = ntl_ZZ_pX([],cup)
    #printer.x = modup
    #print printer
    # Note that we're losing one digit of precision here.
    # This is correct because right shifting does not preserve precision.
    # a is now the negative of the inverse of the unit part of the constant of the defining polynomial (there's a mouthful)
    ZZ_pX_RightShift(tmp, modup, 1)
    ## printer.x = modup
    ## print printer
    ZZ_pX_mul_ZZ_p(tmp, tmp, ZZ_to_ZZ_p(a))
    # tmp is now p/x
    ZZ_pX_conv_modulus(into_multiplier, tmp, prime_pow.get_top_context().x)
    ZZ_pX_Multiplier_construct(low_shifter)
    ZZ_pX_Multiplier_build(low_shifter[0], into_multiplier, prime_pow.get_top_modulus()[0])
    for i from 1 <= i < prime_pow.low_length:
        # Currently tmp = p / x^(2^(i-1)).  Squaring yields p^2 / x^(2^i)
        ZZ_pX_SqrMod(tmp, tmp, modup)
        # Now we divide by p.  We don't really have the extra digit of precision that cup would imply, but we're about to reduce, and then square.
        # This should give us one digit of precision loss, as expected.
        ZZ_pX_right_pshift(tmp, tmp, prime_pow.small_powers[1], cup.x)
        ZZ_pX_conv_modulus(into_multiplier, tmp, prime_pow.get_top_context().x)
        ZZ_pX_Multiplier_construct(&(low_shifter[i]))
        ZZ_pX_Multiplier_build(low_shifter[i], into_multiplier, prime_pow.get_top_modulus()[0])

    # Now we handle high_shifter.
    # We can obtain p/x^e by computing the inverse of x^e/p.
    # Note that modup is still defined from before
    cup.restore_c()

    ZZ_pX_conv_modulus(modup, prime_pow.get_top_modulus()[0].val(), cup.x)
    ZZ_pX_SetCoeff_long(modup, prime_pow.deg, 0)
    ZZ_pX_negate(modup, modup)
    ZZ_pX_right_pshift(into_multiplier, modup, prime_pow.small_powers[1], prime_pow.get_top_context().x)

    # into_multiplier now holds x^e/p
    # prime_pow.c.x should have been restored, but we make sure
    prime_pow.restore_top_context()
    ZZ_pX_InvMod_newton(into_multiplier, into_multiplier, prime_pow.get_top_modulus()[0], prime_pow.get_top_context().x, (<ntl_ZZ_pContext_class>prime_pow.get_context(1)).x)
    ZZ_pX_Multiplier_construct(high_shifter)
    ZZ_pX_Multiplier_build(high_shifter[0], into_multiplier, prime_pow.get_top_modulus()[0])
    # Now we cache powers of p/x^e.  This is a unit, so we don't have to worry about precision issues (yay!)
    for i from 1 <= i < high_length:
        ZZ_pX_SqrMod_pre(into_multiplier, into_multiplier, prime_pow.get_top_modulus()[0])
        ZZ_pX_Multiplier_construct(&(high_shifter[i]))
        ZZ_pX_Multiplier_build(high_shifter[i], into_multiplier, prime_pow.get_top_modulus()[0])

    # I'm not sure whether I need to destruct the temporary variables.
    # ZZ_pX_destruct(&tmp)
    # ZZ_pX_destruct(&modup)
    # ZZ_pX_destruct(&into_multiplier)
    # ZZ_destruct(&a)
