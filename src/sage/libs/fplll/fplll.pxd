include "sage/ext/cdefs.pxi"
include "fplll.pxi"

cdef class FP_LLL:
    cdef object fp_map
    cdef ZZ_mat[mpz_t] *_lattice

cdef extern from "flint/fmpz.h":
    ctypedef void* fmpz_t
    void fmpz_get_mpz(mpz_t x, const fmpz_t f)
    void fmpz_set_mpz(fmpz_t f, const mpz_t val)

cdef extern from "flint/fmpz_mat.h":
    ctypedef void* fmpz_mat_t
    fmpz_t fmpz_mat_entry(fmpz_mat_t mat ,long i ,long j)
