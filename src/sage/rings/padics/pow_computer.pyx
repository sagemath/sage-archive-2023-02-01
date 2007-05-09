import weakref
from sage.rings.infinity import infinity

cimport sage.rings.rational
from sage.rings.rational cimport Rational

include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"

#TODO: doctests and explination of what this file is

cdef class PowComputer_class(SageObject):
    def __new__(self, Integer base, unsigned int log_of_dense_limit):
        # TODO: what here could/should be moved to __init__?
        cdef Py_ssize_t i
        cdef Integer x
        self._initialized = 0
        #if log_of_dense_limit < 1:
        #    raise ValueError, "You can only create a PowComputer with dense_limit >= 2"
        self.log_of_dense_limit = log_of_dense_limit
        self.dense_mask = (1 << self.log_of_dense_limit) #this is a temporary value: it is decremented by 1 below.
        _sig_on
        self.dense_list = <mpz_t *>sage_malloc(sizeof(mpz_t) * self.dense_mask)
        _sig_off
        if self.dense_list == NULL:
            raise MemoryError, "out of memory allocating power storing"
        #self.dense_list_Integer = PyList_New(self.dense_mask)
        self.dense_list_Integer = []
        self.dense_mask -= 1 # fixing self.dense_mask to its correct value
        mpz_init_set_si(self.dense_mask_mpz, self.dense_mask)
        self.basex = base

        mpz_init_set_ui(self.dense_list[0], 1)
        x = PY_NEW(Integer)
        mpz_set_ui(x.value, 1)
        #PyList_SET_ITEM(self.dense_list_Integer, 0, x)
        self.dense_list_Integer.append(x)

        mpz_init_set(self.dense_list[1], self.basex.value)
        x = PY_NEW(Integer)
        mpz_set(x.value, self.basex.value)
        #PyList_SET_ITEM(self.dense_list_Integer, 1, x)
        self.dense_list_Integer.append(x)

        for i from 2 <= i <= self.dense_mask:
            mpz_init(self.dense_list[i])
            mpz_mul(self.dense_list[i], self.dense_list[i - 1], self.basex.value)
            x = PY_NEW(Integer)
            mpz_set(x.value, self.dense_list[i])
            #PyList_SET_ITEM(self.dense_list_Integer, i, x)
            self.dense_list_Integer.append(x)

        self.zero = PY_NEW(Integer)
        mpz_set_ui(self.zero.value, 0)

    def __init__(self, Integer base, unsigned int log_of_dense_limit):
        self._cache = {}
        self._initialized = 1

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._initialized:
            mpz_clear(self.dense_mask_mpz)
            for i from 0 <= i <= self.dense_mask:
                mpz_clear(self.dense_list[i])
            sage_free(self.dense_list)

    def __reduce__(self):
        """
        P = PowComputer(5, 7)
        loads(dumps(P))
        """
        return PowComputer, (self.basex, self.log_of_dense_limit)

    def base(self):
        return self.basex

    cdef void pow_mpz_ui(self, mpz_t ans, unsigned int n):
        cdef int sparse_count, cur_mask, sparse_lim
        cdef mpz_t tmp
        if n <= self.dense_mask:
            mpz_set(ans, self.dense_list[n])
        elif self._cache.has_key(n):
            mpz_set(ans, (<Integer>self._cache[n]).value)
        else:
            sparse_lim = self.log_of_dense_limit
            mpz_init(tmp)
            mpz_set(ans, self.dense_list[n & self.dense_mask])
            n = n >> self.log_of_dense_limit
            while n != 0:
                cur_mask = n & self.dense_mask
                if cur_mask != 0:
                    mpz_set(tmp, self.dense_list[cur_mask])
                    for sparse_count from 0 <= sparse_count < sparse_lim:
                        _sig_on
                        mpz_mul(tmp, tmp, tmp)
                        _sig_off
                    _sig_on
                    mpz_mul(ans, ans, tmp)
                    _sig_off
                n = n >> self.log_of_dense_limit
                sparse_lim += self.log_of_dense_limit
            mpz_clear(tmp)

    cdef void pow_mpz_mpz(self, mpz_t ans, mpz_t n):
        cdef int sparse_count, sparse_lim
        cdef mpz_t tmp, cur_mask
        #cdef Integer printing
        #printing = PY_NEW(Integer)
        #if mpz_cmp_ui(n, 0) < 0:
        #    raise ValueError "n cannot be negative for this function"
        if mpz_cmp_si(n, self.dense_mask) <= 0:
            mpz_set(ans, self.dense_list[mpz_get_si(n)])
        elif mpz_fits_uint_p(n) != 0 and self._cache.has_key(mpz_get_ui(n)):
            mpz_set(ans, (<Integer>self._cache[mpz_get_ui(n)]).value)
        else:
            sparse_lim = self.log_of_dense_limit
            mpz_init(tmp)
            mpz_init(cur_mask)
            mpz_and(cur_mask, n, self.dense_mask_mpz)
            mpz_set(ans, self.dense_list[mpz_get_si(cur_mask)])
            mpz_fdiv_q_2exp(n, n, self.log_of_dense_limit)
            #mpz_set(printing.value, ans)
            #print printing
            #mpz_set(printing.value, n)
            #print printing
            while mpz_cmp_ui(n, 0) != 0:
                mpz_and(cur_mask, n, self.dense_mask_mpz)
                if mpz_cmp_ui(cur_mask, 0) != 0:
                    mpz_set(tmp, self.dense_list[mpz_get_si(cur_mask)])
                    for sparse_count from 0 <= sparse_count < sparse_lim:
                        _sig_on
                        mpz_mul(tmp, tmp, tmp)
                        _sig_off
                    _sig_on
                    mpz_mul(ans, ans, tmp)
                    _sig_off
                mpz_fdiv_q_2exp(n, n, self.log_of_dense_limit)
                sparse_lim += self.log_of_dense_limit
            mpz_clear(tmp)
            mpz_clear(cur_mask)

    def cache(self, n):
        cdef int i
        n = Integer(n)
        if mpz_fits_uint_p((<Integer>n).value) != 0:
            if mpz_cmp((<Integer>n).value, self.dense_mask_mpz) > 0:
                i = n
                self._cache[i] = self(n)
        else:
            raise ValueError, "Too big to cache"

    def get_cache(self):
        return self._cache

    def dense_limit(self):
        return Integer(1 << self.log_of_dense_limit)

    def __call__(self, n):
        cdef Rational r
        cdef Integer z, _n
        cdef mpz_t tmp
        if n is infinity:
            return self.zero
        if not PY_TYPE_CHECK(n, Integer):
            _n = Integer(n)
        else:
            _n = <Integer>n
        if mpz_cmp_ui(_n.value, 0) < 0:
            r = <Rational> PY_NEW(Rational)
            mpz_set_ui(mpq_numref(r.value), 1)
            mpz_init(tmp)
            mpz_neg(tmp, _n.value)
            self.pow_mpz_mpz(mpq_denref(r.value), tmp)
            mpz_clear(tmp)
            return r
        elif mpz_cmp(_n.value, self.dense_mask_mpz) <= 0:
            return self.dense_list_Integer[mpz_get_si(_n.value)]
        else:
            z = PY_NEW(Integer)
            self.pow_mpz_mpz(z.value, _n.value)
            return z

pow_comp_cache = {}
cdef PowComputer_class PowComputer_c(Integer m, int L):
    key = (m, L)
    if pow_comp_cache.has_key(key):
        PC = pow_comp_cache[key]()
        if PC is not None:
            return PC
    PC = PowComputer_class(m, L)
    pow_comp_cache[key] = weakref.ref(PC)
    return PC

def PowComputer(m, L):
    cdef int _L
    m = Integer(m)
    _L = L
    key = (m, _L)
    if pow_comp_cache.has_key(key):
        PC = pow_comp_cache[key]()
        if PC is not None:
            return PC
    PC = PowComputer_class(m, _L)
    pow_comp_cache[key] = weakref.ref(PC)
    return PC

