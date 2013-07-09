include "fmpz.pxi"

cdef extern from "flint/padic.h":

    ctypedef void * padic_t
    ctypedef void * padic_ctx_t
    ctypedef void * padic_inv_t

    cdef enum padic_print_mode:
        PADIC_TERSE
        PADIC_SERIES
        PADIC_VAL_UNIT

    ctypedef struct padic_ctx_struct:
        fmpz_t p
        long N
        double pinv
        fmpz* pow
        long min
        long max

    fmpz_t padic_unit(padic_t)
    long padic_val(padic_t)

    void padic_ctx_init(padic_ctx_t ctx, fmpz_t p, long N, padic_print_mode mode)
    void padic_ctx_clear(padic_ctx_t ctx)
    int _padic_ctx_pow_ui(fmpz_t rop, unsigned long e, padic_ctx_t ctx)
