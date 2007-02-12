include "../ext/cdefs.pxi"

from sage.rings.integer import Integer
from sage.rings.integer cimport Integer


cdef extern from "multi_modular.h":
    ctypedef unsigned long mod_int
    mod_int MOD_INT_MAX
    mod_int START_PRIME_MAX

cdef class MultiModularBasis_base:
    cdef int      n
    cdef mod_int* moduli
    cdef mpz_t*   partial_products
    cdef mod_int* C # precomputed values for CRT

    cdef mod_int last_prime(self)
    cdef int _extend_moduli_to_height_c(self, mpz_t height) except -1
    cdef void _refresh_products(self, int start)
    cdef void _refresh_precomputations(self, int start)
    cdef int min_moduli_count(self, mpz_t height) except -1
    cdef int moduli_list_c(self, mod_int** moduli)

    cdef int mpz_reduce_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1
    cdef int mpz_reduce_vec_tail(self, mpz_t* z, mod_int** b, int vn, int offset, int len) except -1
    cdef int mpz_crt_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1
    cdef int mpz_crt_vec_tail(self, mpz_t* z, mod_int** b, int vn, int offset, int len) except -1

cdef class MultiModularBasis(MultiModularBasis_base):
    cdef int mpz_reduce(self, mpz_t z, mod_int* b) except -1
    cdef int mpz_reduce_vec(self, mpz_t* z, mod_int** b, int vn) except -1
    cdef int mpz_crt(self, mpz_t z, mod_int* b) except -1
    cdef int mpz_crt_vec(self, mpz_t* z, mod_int** b, int vn) except -1

cdef class MutableMultiModularBasis(MultiModularBasis):
    cdef mod_int __last_prime
    cdef mod_int next_prime_c(self) except -1
    cdef mod_int replace_prime_c(self, int i) except -1
