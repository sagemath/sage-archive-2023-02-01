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
    def __new__(self, Integer prime, unsigned long cache_limit):
        self._initialized = 0
        _sig_on
        self.dense_list = <mpz_t *>sage_malloc(sizeof(mpz_t) * (cache_limit + 1))
        _sig_off
        if self.dense_list == NULL:
            raise MemoryError, "out of memory allocating power storing"
        mpz_init(self.modulus)

    def __init__(self, Integer prime, unsigned long cache_limit):
        cdef Py_ssize_t i
        cdef Integer x

        self._cache_limit = cache_limit
        self.prime = prime
        self.in_field = 0

        #self.dense_list_Integer = PyList_New(self.cache_limit + 1)
        self.dense_list_Integer = []

        mpz_init_set_ui(self.dense_list[0], 1)
        x = PY_NEW(Integer)
        mpz_set_ui(x.value, 1)
        #PyList_SET_ITEM(self.dense_list_Integer, 0, x)
        self.dense_list_Integer.append(x)

        mpz_init_set(self.dense_list[1], prime.value)
        x = PY_NEW(Integer)
        mpz_set(x.value, prime.value)
        #PyList_SET_ITEM(self.dense_list_Integer, 1, x)
        self.dense_list_Integer.append(x)

        for i from 2 <= i <= cache_limit:
            mpz_init(self.dense_list[i])
            mpz_mul(self.dense_list[i], self.dense_list[i - 1], prime.value)
            x = PY_NEW(Integer)
            mpz_set(x.value, self.dense_list[i])
            #PyList_SET_ITEM(self.dense_list_Integer, i, x)
            self.dense_list_Integer.append(x)
        mpz_set(self.modulus, self.dense_list[cache_limit])
        self._initialized = 1

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i <= self._cache_limit:
                mpz_clear(self.dense_list[i])
            sage_free(self.dense_list)
            mpz_clear(self.modulus)

    def __cmp__(self, other):
        if isinstance(other, PowComputer_class):
            if self.prime < (<PowComputer_class>other).prime:
                return -1
            elif self.prime > (<PowComputer_class>other).prime:
                return 1
            elif self._cache_limit < (<PowComputer_class>other)._cache_limit:
                return -1
            elif self._cache_limit > (<PowComputer_class>other)._cache_limit:
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
        return PowComputer, (self.prime, self._cache_limit)

    def _prime(self):
        """
        Returns the base that the PowComputer is exponentiating.

        EXAMPLES:
        sage: P = PowComputer(6, 10)
        sage: P._prime()
        6
        """
        return self.prime

    def cache_limit(self):
        """
        Returns the limit to which powers of prime are computed.
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self._cache_limit)
        return ans

    def _modulus(self):
        """
        Returns self.prime()^self.cache_limit()
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.modulus)
        return ans

    def _in_field(self):
        """
        For use by p-adic rings and fields.  Feel free to ignore if you're using a PowComputer otherwise.
        """
        if self.in_field == 1:
            return True
        else:
            return False

    def _set_in_field(self, value):
        """
        For use by p-adic rings and fields.  Feel free to ignore if you're using a PowComputer otherwise.
        """
        if value:
            self.in_field = 1
        else:
            self.in_field = 0

    #cdef void pow_mpz_ui(self, mpz_t ans, unsigned int n):
    #    cdef int sparse_count, cur_mask, sparse_lim
    #    cdef mpz_t tmp
    #    if n <= self.dense_mask:
    #        mpz_set(ans, self.dense_list[n])
    #    elif self._cache.has_key(n):
    #        mpz_set(ans, (<Integer>self._cache[n]).value)
    #    else:
    #        sparse_lim = self.log_of_dense_limit
    #        mpz_init(tmp)
    #        mpz_set(ans, self.dense_list[n & self.dense_mask])
    #        n = n >> self.log_of_dense_limit
    #        while n != 0:
    #            cur_mask = n & self.dense_mask
    #            if cur_mask != 0:
    #                mpz_set(tmp, self.dense_list[cur_mask])
    #                for sparse_count from 0 <= sparse_count < sparse_lim:
    #                    _sig_on
    #                    mpz_mul(tmp, tmp, tmp)
    #                    _sig_off
    #                _sig_on
    #                mpz_mul(ans, ans, tmp)
    #                _sig_off
    #            n = n >> self.log_of_dense_limit
    #            sparse_lim += self.log_of_dense_limit
    #        mpz_clear(tmp)

    #cdef void pow_mpz_mpz(self, mpz_t ans, mpz_t n):
    #    cdef int sparse_count, sparse_lim
    #    cdef mpz_t tmp, cur_mask
    #    #cdef Integer printing
    #    #printing = PY_NEW(Integer)
    #    #if mpz_cmp_ui(n, 0) < 0:
    #    #    raise ValueError "n cannot be negative for this function"
    #    if mpz_cmp_si(n, self.dense_mask) <= 0:
    #        mpz_set(ans, self.dense_list[mpz_get_si(n)])
    #    elif mpz_fits_uint_p(n) != 0 and self._cache.has_key(mpz_get_ui(n)):
    #        mpz_set(ans, (<Integer>self._cache[mpz_get_ui(n)]).value)
    #    else:
    #        sparse_lim = self.log_of_dense_limit
    #        mpz_init(tmp)
    #        mpz_init(cur_mask)
    #        mpz_and(cur_mask, n, self.dense_mask_mpz)
    #        mpz_set(ans, self.dense_list[mpz_get_si(cur_mask)])
    #        mpz_fdiv_q_2exp(n, n, self.log_of_dense_limit)
    #        #mpz_set(printing.value, ans)
    #        #print printing
    #        #mpz_set(printing.value, n)
    #        #print printing
    #        while mpz_cmp_ui(n, 0) != 0:
    #            mpz_and(cur_mask, n, self.dense_mask_mpz)
    #            if mpz_cmp_ui(cur_mask, 0) != 0:
    #                mpz_set(tmp, self.dense_list[mpz_get_si(cur_mask)])
    #                for sparse_count from 0 <= sparse_count < sparse_lim:
    #                    _sig_on
    #                    mpz_mul(tmp, tmp, tmp)
    #                    _sig_off
    #                _sig_on
    #                mpz_mul(ans, ans, tmp)
    #                _sig_off
    #            mpz_fdiv_q_2exp(n, n, self.log_of_dense_limit)
    #            sparse_lim += self.log_of_dense_limit
    #        mpz_clear(tmp)
    #        mpz_clear(cur_mask)

    #def cache(self, n):
    #    cdef int i
    #    n = Integer(n)
    #    if mpz_fits_uint_p((<Integer>n).value) != 0:
    #        if mpz_cmp((<Integer>n).value, self.dense_mask_mpz) > 0:
    #            i = n
    #            self._cache[i] = self(n)
    #    else:
    #        raise ValueError, "Too big to cache"

    #def get_cache(self):
    #    return self._cache

    #def dense_limit(self):
    #    return Integer(1 << self.log_of_dense_limit)

    def __call__(self, n):
        cdef Integer z, _n
        cdef mpz_t tmp
        if n is infinity:
            return Integer(0)
        if not PY_TYPE_CHECK(n, Integer):
            _n = Integer(n)
        else:
            _n = <Integer>n
        if mpz_cmp_ui(_n.value, 0) >= 0 and mpz_cmp_ui(_n.value, self._cache_limit) <= 0:
            return self.dense_list_Integer[mpz_get_ui(_n.value)]
        else:
            return self.prime.__pow__(_n)

pow_comp_cache = {}
cdef PowComputer_class PowComputer_c(Integer m, Integer L, in_field):
    if L < 0:
        raise ValueError, "L must be non-negative."
    key = (m, L, in_field)
    if pow_comp_cache.has_key(key):
        PC = pow_comp_cache[key]()
        if PC is not None:
            return PC
    PC = PowComputer_class(m, mpz_get_ui(L.value))
    PC._set_in_field(in_field)
    pow_comp_cache[key] = weakref.ref(PC)
    return PC

# To speed up the creation of PowComputers with the same m, we might eventually want to copy over data from an existing PowComputer.

def PowComputer(m, L, in_field = False):
    """
    Returns a PowComputer that caches the values $1, m, m^2, \ldots, m^L$.

    Once you create a PowComputer, merely call it to get values out.
    You can input any integer, even if it's outside of the precomputed range.

    INPUT:
    m -- An integer, the base that you want to exponentiate.
    L -- A positive integer that you want to cache powers up to.

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
    if not PY_TYPE_CHECK(L, Integer):
        L = Integer(L)
    if L < 0:
        raise ValueError, "L must be positive"
    return PowComputer_c(m, L, in_field)

