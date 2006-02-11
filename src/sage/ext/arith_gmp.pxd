include "cdefs.pxi"

cdef class functions:
    cdef int mpz_crt(self, mpz_t z, mpz_t x, mpz_t y, mpz_t m, mpz_t n) except -1
    cdef int mpz_vec(self, mpz_t **z, int *v, int n) except -1
    cdef int mpzvec_to_intmod(self, unsigned long **z, mpz_t *v, int n, unsigned long p) except -1
    cdef int intmodvec_to_mpz(self, mpz_t **z, unsigned long *v, int n) except -1
    cdef int allocate_mpz_zero_array(self, mpz_t **z, int n) except -1
    cdef int mpzvec_clear(self, mpz_t *z, int n) except -1
    cdef int mpz_height_vec(self, mpz_t H, mpz_t *v, int n) except -1

