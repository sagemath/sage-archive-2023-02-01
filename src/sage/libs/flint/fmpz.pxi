include "../ntl/decl.pxi"

cdef extern from "FLINT/fmpz.h":

    ctypedef void * fmpz_t
    ctypedef void * mpz_t

    fmpz_t fmpz_init(unsigned long limbs)
    void fmpz_clear(fmpz_t f)
    void fmpz_print(fmpz_t f)
    int fmpz_is_one(fmpz_t f)

    void fmpz_add_ui_inplace(fmpz_t output, unsigned long x)
    void fmpz_sub_ui_inplace(fmpz_t output, unsigned long x)

    void fmpz_to_mpz(mpz_t rop, fmpz_t op)
