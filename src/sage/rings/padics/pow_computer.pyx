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

#from sage.rings.integer import Integer

include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"

cdef class PowComputer_class(SageObject):
    def __init__(self, Integer prime, bint in_field):
        self.prime = prime
        self.in_field = in_field
        self._initialized = 1

    cdef Integer pow_Integer(self, unsigned long n):
        cdef Integer ans = PY_NEW(Integer)
        _sig_on
        mpz_pow_ui(ans.value, self.prime.value, n)
        _sig_off
        return ans

    cdef mpz_t pow_mpz_t(self, unsigned long n):
        raise NotImplementedError

    cdef mpz_t pow_mpz_t_tmp(self, unsigned long n):
        raise NotImplementedError

    cdef ZZ_c pow_ZZ(self, unsigned long n):
        raise NotImplementedError

    def _prime(self):
        """
        Returns the base that the PowComputer is exponentiating.

        EXAMPLES:
        sage: P = PowComputer(6, 10)
        sage: P._prime()
        6
        """
        return self.prime

    def _in_field(self):
        """
        For use by p-adic rings and fields.  Feel free to ignore if you're using a PowComputer otherwise.
        """
        return self.in_field

    def __call__(self, n):
        cdef Integer z, _n
        cdef mpz_t tmp
        if n is infinity:
            return Integer(0)
        if not PY_TYPE_CHECK(n, Integer):
            _n = Integer(n)
        else:
            _n = <Integer>n
        if mpz_cmp_ui(_n.value, 0) >= 0:
            return self.pow_Integer(mpz_get_ui(_n.value))
        else:
            return self.prime.__pow__(_n)

cdef class PowComputer_base(PowComputer_class):
    def __new__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field):
        (<PowComputer_class>self)._initialized = 0
        _sig_on
        self.small_powers = <mpz_t *>sage_malloc(sizeof(mpz_t) * (cache_limit + 1))
        _sig_off
        if self.small_powers == NULL:
            raise MemoryError, "out of memory allocating power storing"
        mpz_init(self.top_power)
        mpz_init(self.temp)

    def __init__(self, Integer prime, unsigned long cache_limit, unsigned long prec_cap, bint in_field):
        cdef Py_ssize_t i
        cdef Integer x

        self.cache_limit = cache_limit
        self.prec_cap = prec_cap

        mpz_init_set_ui(self.small_powers[0], 1)

        if cache_limit > 0:
            mpz_init_set(self.small_powers[1], prime.value)

        for i from 2 <= i <= cache_limit:
            mpz_init(self.small_powers[i])
            mpz_mul(self.small_powers[i], self.small_powers[i - 1], prime.value)
        mpz_pow_ui(self.top_power, prime.value, prec_cap)
        PowComputer_class.__init__(self, prime, in_field)
        (<PowComputer_class>self)._initialized = 1

    def __dealloc__(self):
        cdef Py_ssize_t i
        if (<PowComputer_class>self)._initialized:
            for i from 0 <= i <= self.cache_limit:
                mpz_clear(self.small_powers[i])
            sage_free(self.small_powers)
            mpz_clear(self.top_power)
            mpz_clear(self.temp)

    def __cmp__(self, other):
        if isinstance(other, PowComputer_base):
            if self.prime < (<PowComputer_class>other).prime:
                return -1
            elif self.prime > (<PowComputer_class>other).prime:
                return 1
            elif self.prec_cap < (<PowComputer_base>other).prec_cap:
                return -1
            elif self.prec_cap > (<PowComputer_base>other).prec_cap:
                return 1
            elif self.cache_limit < (<PowComputer_base>other).cache_limit:
                return -1
            elif self.cache_limit > (<PowComputer_base>other).cache_limit:
                return 1
            elif self.in_field < (<PowComputer_class>other).in_field:
                return -1
            elif self.in_field > (<PowComputer_class>other).in_field:
                return 1
            else:
                return 0
        else:
            return cmp(type(self), type(other))

    def __reduce__(self):
        """
        sage: P = PowComputer(5, 7)
        sage: Q = PowComputer(5, 10)
        sage: R = loads(dumps(P))
        sage: P == R
        True
        """
        return PowComputer, (self.prime, self.cache_limit, self.prec_cap, self.in_field)

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
            mpz_set(ans.value, self.small_powers[n])
        elif n == self.prec_cap:
            mpz_set(ans.value, self.top_power)
        else:
            mpz_pow_ui(ans.value, self.prime.value, n)
        return ans

    cdef mpz_t pow_mpz_t(self, unsigned long n):
        #################### WARNING ######################
        ## If you use this function, you MAY need to     ##
        ## call mpz_clear on the returned value after    ##
        ## you are finished with it.  You need to do     ##
        ## so if and only if (n > self.cache_limit and   ##
        ## n != self.prec_cap).  So your code should     ##
        ## look something like the following:            ##
        ##                                               ##
        ## cdef PowComputer_base ppow                    ##
        ## ppow = PowComputer(5, 10, 1000)               ##
        ## cdef unsigned long n = get_long_somehow()     ##
        ## cdef mpz_t foo = ppow.pow_mpz_t(n)            ##
        ## # Do stuff with foo                           ##
        ## # Now done with foo                           ##
        ## if n > ppow.cache_limit and n != ppow.prec_cap##
        ##     mpz_clear(foo)                            ##
        ###################################################
        if n <= self.cache_limit:
            return self.small_powers[n]
        if n == self.prec_cap:
            return self.top_power
        cdef mpz_t ans
        mpz_init(ans)
        mpz_pow_ui(ans, self.prime.value, n)
        return ans

    cdef mpz_t pow_mpz_t_tmp(self, unsigned long n):
        ## Solves the problem noted in the above warning by storing the returned mpz_t in a temporary variable
        ## Each call to this function overwrites that temporary variable.
        if n <= self.cache_limit:
            return self.small_powers[n]
        if n == self.prec_cap:
            return self.top_power
        mpz_pow_ui(self.temp, self.prime.value, n)
        return self.temp


    cdef ZZ_c pow_ZZ(self, unsigned long n):
        ## You always need to call ZZ_destruct on the return value of this function
        cdef ZZ_c ans
        cdef mpz_t temp
        ZZ_construct(&ans)
        if n <= self.cache_limit:
            mpz_to_ZZ(&ans, &(self.small_powers[n]))
        elif n == self.prec_cap:
            mpz_to_ZZ(&ans, &self.top_power)
        else:
            mpz_init(temp)
            mpz_pow_ui(temp, self.prime.value, n)
            mpz_to_ZZ(&ans, &temp)
            mpz_clear(temp)
        return ans


pow_comp_cache = {}
cdef PowComputer_base PowComputer_c(Integer m, Integer cache_limit, Integer prec_cap, in_field):
    if cache_limit < 0:
        raise ValueError, "cache_limit must be non-negative."
    if prec_cap < 0:
        raise ValueError, "prec_cap must be non-negative."
    key = (m, cache_limit, prec_cap, in_field)
    if pow_comp_cache.has_key(key):
        PC = pow_comp_cache[key]()
        if PC is not None:
            return PC
    PC = PowComputer_class(m, mpz_get_ui(cache_limit.value), mpz_get_ui(prec_cap.value), in_field)
    pow_comp_cache[key] = weakref.ref(PC)
    return PC

# To speed up the creation of PowComputers with the same m, we might eventually want to copy over data from an existing PowComputer.

def PowComputer(m, cache_limit, prec_cap, in_field = False):
    """
    Returns a PowComputer that caches the values $1, m, m^2, \ldots, m^cache_limit$.

    Once you create a PowComputer, merely call it to get values out.
    You can input any integer, even if it's outside of the precomputed range.

    INPUT:
    m -- An integer, the base that you want to exponentiate.
    cache_limit -- A positive integer that you want to cache powers up to.

    EXAMPLES:
    sage: PC = PowComputer(3, 5)
    sage: PC(4)
    81
    sage: PC(6)
    729
    sage: PC(-1)
    1/3
    """
    if not PY_TYPE_CHECK(m, Integer):
        m = Integer(m)
    if not PY_TYPE_CHECK(cache_limit, Integer):
        cache_limit = Integer(cache_limit)
    if not PY_TYPE_CHECK(prec_cap, Integer):
        prec_cap = Integer(prec_cap)
    return PowComputer_c(m, cache_limit, prec_cap, in_field)

