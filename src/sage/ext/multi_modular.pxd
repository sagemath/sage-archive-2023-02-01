include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"


cdef extern from "multi_modular.h":
    ctypedef unsigned long mod_int
    mod_int MOD_INT_MAX
    mod_int START_PRIME_MAX

cdef class MultiModularBasis:
    cdef int      moduli_count
    cdef mod_int* moduli
    cdef mpz_t*   moduli_partial_product
    cdef mod_int* C # precomputed values for CRT

    cdef int _extend_moduli_list(self, mpz_t height) except -1
    cdef int moduli_list_c(self, mod_int** moduli, mpz_t height) except -1
    cdef int mpz_crt(self, mpz_t* z, mod_int* b, int n) except -1
    cdef int mpz_crt_vec(self, mpz_t** z, mod_int** b, int n, int vn) except -1
