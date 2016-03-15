# distutils: libraries = flint
# This file was (manually) generated from FLINT's padic.h.
#*****************************************************************************
#       Copyright (C) 2011 Sebastian Pancratz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.flint.types cimport *

cdef extern from "flint/padic.h":
    fmpz_t padic_unit(padic_t)
    long padic_val(padic_t)

    # Context ******************************************************************
    void padic_ctx_init(padic_ctx_t ctx, fmpz_t p, long N, padic_print_mode mode)
    void padic_ctx_clear(padic_ctx_t ctx)
    int _padic_ctx_pow_ui(fmpz_t rop, unsigned long e, padic_ctx_t ctx)

    # Memory management ********************************************************
    void _padic_init(padic_t rop)
    void padic_init(padic_t rop, padic_ctx_t ctx)
    void _padic_clear(padic_t rop)
    void padic_clear(padic_t rop, padic_ctx_t ctx)
    void _padic_canonicalise(padic_t rop, padic_ctx_t ctx)
    void _padic_reduce(padic_t rop, padic_ctx_t ctx)
    void padic_reduce(padic_t rop, padic_ctx_t ctx)

    # Assignments and conversions **********************************************
    void _padic_set(padic_t rop, padic_t op)
    void padic_set(padic_t rop, padic_t op, padic_ctx_t ctx)
    void _padic_set_si(padic_t rop, long op, padic_ctx_t ctx)
    void padic_set_si(padic_t rop, long op, padic_ctx_t ctx)
    void _padic_set_ui(padic_t rop, unsigned long op, padic_ctx_t ctx)
    void padic_set_ui(padic_t rop, unsigned long op, padic_ctx_t ctx)
    void _padic_set_fmpz(padic_t rop, fmpz_t op, padic_ctx_t ctx)
    void padic_set_fmpz(padic_t rop, fmpz_t op, padic_ctx_t ctx)
    void padic_set_fmpq(padic_t rop, fmpq_t op, padic_ctx_t ctx)
    void _padic_set_mpz(padic_t rop, mpz_t op, padic_ctx_t ctx)
    void padic_set_mpz(padic_t rop, mpz_t op, padic_ctx_t ctx)
    void padic_set_mpq(padic_t rop, mpq_t op, padic_ctx_t ctx)
    void _padic_get_fmpz(fmpz_t rop, padic_t op, padic_ctx_t ctx)
    void padic_get_fmpz(fmpz_t rop, padic_t op, padic_ctx_t ctx)
    void _padic_get_fmpq(fmpq_t rop, padic_t op, padic_ctx_t ctx)
    void padic_get_fmpq(fmpq_t rop, padic_t op, padic_ctx_t ctx)
    void _padic_get_mpz(mpz_t rop, padic_t op, padic_ctx_t ctx)
    void padic_get_mpz(mpz_t rop, padic_t op, padic_ctx_t ctx)
    void _padic_get_mpq(mpq_t rop, padic_t op, padic_ctx_t ctx)
    void padic_get_mpq(mpq_t rop, padic_t op, padic_ctx_t ctx)
    void padic_swap(padic_t op1, padic_t op2)
    void padic_zero(padic_t rop)
    void _padic_one(padic_t rop)
    void padic_one(padic_t rop, padic_ctx_t ctx)

    # Arithmetic operations ****************************************************
    void _padic_add(padic_t rop, padic_t op1, padic_t op2, padic_ctx_t ctx)
    void padic_add(padic_t rop, padic_t op1, padic_t op2, padic_ctx_t ctx)
    void _padic_sub(padic_t rop, padic_t op1, padic_t op2, padic_ctx_t ctx)
    void padic_sub(padic_t rop, padic_t op1, padic_t op2, padic_ctx_t ctx)
    void _padic_neg(padic_t rop, padic_t op)
    void padic_neg(padic_t rop, padic_t op, padic_ctx_t ctx)
    void _padic_mul(padic_t rop, padic_t op1, padic_t op2)
    void padic_mul(padic_t rop, padic_t op1, padic_t op2, padic_ctx_t ctx)
    void padic_shift(padic_t rop, padic_t op, long v, padic_ctx_t ctx)
    void padic_div(padic_t rop, padic_t op1, padic_t op2, padic_ctx_t ctx)
    void _padic_inv_precompute(padic_inv_t S, fmpz_t p, long N)
    void _padic_inv_clear(padic_inv_t S)
    void _padic_inv_precomp(fmpz_t rop, fmpz_t op, padic_inv_t S)
    void _padic_inv(fmpz_t rop, fmpz_t op, fmpz_t p, long N)
    void padic_inv(padic_t rop, padic_t op, padic_ctx_t ctx)
    int padic_sqrt(padic_t rop, padic_t op, padic_ctx_t ctx)
    void padic_pow_si(padic_t rop, padic_t op, long e, padic_ctx_t ctx)

    # Comparison ***************************************************************
    int _padic_is_zero(padic_t op)
    int padic_is_zero(padic_t op, padic_ctx_t ctx)
    int _padic_is_one(padic_t op)
    int padic_is_one(padic_t op, padic_ctx_t ctx)
    int _padic_equal(padic_t op1, padic_t op2)
    int padic_equal(padic_t op1, padic_t op2, padic_ctx_t ctx)

    # Special functions ********************************************************
    void _padic_teichmuller(fmpz_t rop, fmpz_t op, fmpz_t p, long N)
    void padic_teichmuller(padic_t rop, padic_t op, padic_ctx_t ctx)
    void _padic_exp(fmpz_t rop, fmpz_t u, long v, fmpz_t p, long N)
    void _padic_exp_naive(fmpz_t rop, fmpz_t u, long v, fmpz_t p, long N)
    void _padic_exp_rectangular(fmpz_t rop, fmpz_t u, long v, fmpz_t p, long N)
    void _padic_exp_balanced(fmpz_t rop, fmpz_t u, long v, fmpz_t p, long N)
    int padic_exp(padic_t rop, padic_t op, padic_ctx_t ctx)
    int padic_exp_rectangular(padic_t rop, padic_t op, padic_ctx_t ctx)
    int padic_exp_balanced(padic_t rop, padic_t op, padic_ctx_t ctx)
    long _padic_log_bound(long v, long N, fmpz_t p)
    void _padic_log(fmpz_t z, fmpz_t y, long v, fmpz_t p, long N)
    void _padic_log_rectangular(fmpz_t z, fmpz_t y, long v, fmpz_t p, long N)
    void _padic_log_satoh(fmpz_t z, fmpz_t y, long v, fmpz_t p, long N)
    void _padic_log_balanced(fmpz_t z, fmpz_t y, long v, fmpz_t p, long N)
    int padic_log(padic_t rop, padic_t op, padic_ctx_t ctx)
    int padic_log_rectangular(padic_t rop, padic_t op, padic_ctx_t ctx)
    int padic_log_satoh(padic_t rop, padic_t op, padic_ctx_t ctx)
    int padic_log_balanced(padic_t rop, padic_t op, padic_ctx_t ctx)
    unsigned long padic_val_fac_ui2(unsigned long N)
    unsigned long padic_val_fac_ui(unsigned long N, fmpz_t p)
    void padic_val_fac(fmpz_t rop, fmpz_t op, fmpz_t p)

    # Input and output *********************************************************
    char * _padic_get_str(char * str, padic_t op, padic_ctx_t ctx)
    char * padic_get_str(char * str, padic_t op, padic_ctx_t ctx)
    int _padic_fprint(FILE * file, fmpz_t u, long v, padic_ctx_t ctx)
    int padic_fprint(FILE * file, padic_t op, padic_ctx_t ctx)
    int padic_fprint(FILE * file, padic_t op, padic_ctx_t ctx)
    int _padic_print(fmpz_t u, long v, padic_ctx_t ctx)
    int padic_print(padic_t op, padic_ctx_t ctx)
    void _padic_debug(padic_t op)
