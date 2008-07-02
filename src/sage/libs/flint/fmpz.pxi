include "../ntl/decl.pxi"

cdef extern from "FLINT/fmpz.h":

    ctypedef void* fmpz_t

    fmpz_t fmpz_init(unsigned long limbs)
    void fmpz_clear(fmpz_t f)
    void fmpz_print(fmpz_t f)
    int fmpz_is_one(fmpz_t f)
