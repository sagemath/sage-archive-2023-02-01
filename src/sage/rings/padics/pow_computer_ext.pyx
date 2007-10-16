from sage.rings.integer cimport Integer

"""
A class for computing and caching powers of the same integer.

This class is designed to be used as a field of p-adic rings and fields.  Since elements of p-adic rings and fields need to use powers of p over and over, this class precomputes and stores powers of p.  There is no reason that the base has to be prime however.

EXAMPLES:
sage: X = PowComputer(3, 4)
sage: X(3)
27
sage: X(20) == 3^20
True

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



import weakref
from sage.rings.infinity import infinity
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext, ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ import ntl_ZZ

#from sage.rings.integer import Integer

include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"

cdef class PowComputer_ext(PowComputer_class):
    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        (<PowComputer_class>self)._initialized = 0
        _sig_on
        self.small_powers = <ZZ_c *>sage_malloc(sizeof(ZZ_c) * (cache_limit + 1))
        _sig_off
        if self.small_powers == NULL:
            raise MemoryError, "out of memory allocating power storing"
        ZZ_construct(&self.top_power)

        cdef Py_ssize_t i
        cdef Integer x

        self.cache_limit = cache_limit
        self.prec_cap = prec_cap

        ZZ_construct(self.small_powers)
        ZZ_conv_int(self.small_powers[0], 1)

        if cache_limit > 0:
            ZZ_construct(&(self.small_powers[1]))
            mpz_to_ZZ(&(self.small_powers[1]), &prime.value)

        _sig_on
        for i from 2 <= i <= cache_limit:
            ZZ_construct(&(self.small_powers[i]))
            mul_ZZ(self.small_powers[i], self.small_powers[i-1], self.small_powers[1])
        mpz_to_ZZ(&self.top_power, &prime.value)
        power_ZZ(self.top_power, self.top_power, prec_cap)
        _sig_off


    def __init__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        PowComputer_class.__init__(self, prime, in_field)

    def __dealloc__(self):
        cdef Py_ssize_t i
        if (<PowComputer_class>self)._initialized:
            for i from 0 <= i <= self.cache_limit:
                ZZ_destruct(&(self.small_powers[i]))
            sage_free(self.small_powers)
            ZZ_destruct(&self.top_power)

    def _cache_limit(self):
        """
        Returns the limit to which powers of prime are computed.
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.cache_limit)
        return ans

    def _top_power(self):
        """
        Returns self._prime()^self._prec_cap()
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.top_power)
        return ans

    def _prec_cap(self):
        """
        Returns prec_cap, a single value that for which self._prime()^prec_cap is stored
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.prec_cap)
        return ans

    cdef Integer pow_Integer(self, unsigned long n):
        cdef Integer ans = PY_NEW(Integer)
        if n <= self.cache_limit:
            ZZ_to_mpz(&ans.value, &(self.small_powers[n]))
        elif n == self.prec_cap:
            ZZ_to_mpz(&ans.value, &self.top_power)
        else:
            mpz_pow_ui(ans.value, self.prime.value, n)
        return ans

    cdef mpz_t pow_mpz_t(self, unsigned long n):
        ## You always need to call mpz_clear on the return value of this function
        cdef mpz_t ans
        mpz_init(ans)
        if n <= self.cache_limit:
            ZZ_to_mpz(&ans, &(self.small_powers[n]))
        elif n == self.prec_cap:
            ZZ_to_mpz(&ans, &self.top_power)
        else:
            mpz_pow_ui(ans, self.prime.value, n)
        return ans

    cdef ZZ_c pow_ZZ(self, unsigned long n):
        #################### WARNING ######################
        ## If you use this function, you MAY need to     ##
        ## call ZZ_destruct on the returned value after  ##
        ## you are finished with it.  You need to do     ##
        ## so if and only if (n > self.cache_limit and   ##
        ## n != self.prec_cap).  So your code should     ##
        ## look something like the following:            ##
        ##                                               ##
        ## cdef PowComputer_ext ppow                     ##
        ## ppow = PowComputer(5, 10, 1000, poly)         ##
        ## cdef unsigned long n = get_long_somehow()     ##
        ## cdef ZZ_c foo = ppow.pow_ZZ(n)                ##
        ## # Do stuff with foo                           ##
        ## # Now done with foo                           ##
        ## if n > ppow.cache_limit and n != ppow.prec_cap##
        ##     ZZ_destruct(foo)                          ##
        ###################################################
        if n <= self.cache_limit:
            return self.small_powers[n]
        if n == self.prec_cap:
            return self.top_power
        cdef ZZ_c ans
        ZZ_construct(&ans)
        power_ZZ(ans, self.small_powers[1], n)
        return ans

cdef class PowComputer_ZZ_pX(PowComputer_ext):
    cdef ntl_ZZ_pContext get_context(self, unsigned long n):
        cdef ZZ_c pn = self.pow_ZZ(n)
        cdef ntl_ZZ_pContext_class context = <ntl_ZZ_pContext_class> ntl_ZZ_pContext(pn)
        if n > self.cache_limit and n != self.prec_cap:
            ZZ_destruct(pn)
        return context

    cdef ntl_ZZ_pContext get_top_context(self):
        return ntl_ZZ_pContext(self.top_power)

    cdef void restore_context(self, unsigned long n):
        self.get_context(n).restore_c()

    cdef void restore_top_context(self):
        self.get_top_context().restore_c()

    cdef ZZ_pX_Modulus_c get_modulus(self, unsigned long n):
        raise NotImplementedError

    cdef ZZ_pX_Modulus_c get_top_modulus(self):
        raise NotImplementedError

cdef class PowComputer_ZZ_pX_FM(PowComputer_ZZ_pX):
    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        # The __new__ method for PowComputer_ext has already run, so we have access to small_powers, top_power.

        # We use ntl_ZZ_pContexts so that contexts are cached centrally.
        self.c = ntl_ZZ_pContext(self.top_power)
        self.c.restore_c()

        # For now, we don't do anything complicated with poly
        if PY_TYPE_CHECK(poly, ntl_ZZ_pX) and (<ntl_ZZ_pX>poly).c is self.c:
            self.poly = (<ntl_ZZ_pX>poly).x
        else:
            raise NotImplementedError
        ZZ_pX_Modulus_construct(&self.mod)
        ZZ_pX_Modulus_build(self.mod, self.poly)

    def __init__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field, poly):
        PowComputer_ZZ_pX.__init__(self, prime, cache_limit, prec_cap, in_field, poly)

    cdef ntl_ZZ_pContext get_top_context(self):
        return self.c

    cdef void restore_top_context(self):
        self.c.restore_c()

    cdef ZZ_pX_Modulus_c get_modulus(self, unsigned long n):
        ## You always need to call ZZ_pX_Modulus_destruct on the return value of this function once you're done.
        cdef ZZX_c temp
        cdef ZZ_pX_c tempZZ_pX
        self.c.restore_c()
        ZZX_construct(&temp)
        ZZ_pX_to_ZZX(temp, self.poly)
        self.restore_context(n)
        ZZ_pX_construct(&tempZZ_pX)
        ZZX_to_ZZ_pX(tempZZ_pX, temp)
        cdef ZZ_pX_Modulus ans
        ZZ_pX_Modulus_construct(&ans)
        ZZ_pX_Modulus_build(ans, tempZZ_pX)
        ZZX_destruct(&temp)
        ZZ_pX_destruct(&tempZZ_pX)
        return ans

    cdef ZZ_pX_Modulus_c get_top_modulus(self):
        ## Never call ZZ_pX_Modulus_destruct on the return value of this function
        return self.mod


