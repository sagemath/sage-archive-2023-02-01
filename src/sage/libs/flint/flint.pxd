cdef extern from "flint/flint.h":
    cdef unsigned long FLINT_BIT_COUNT(unsigned long)
    void flint_free(void * ptr)

cdef extern from "flint/fmpz.h":
    void _fmpz_cleanup()
    void _fmpz_cleanup_mpz_content()
