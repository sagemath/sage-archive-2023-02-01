###############################################################################
#          Copyright (C) 2010 Sebastian Pancratz <sfp@pancratz.org>           #
#                                                                             #
#     Distributed under the terms of the GNU General Public License (GPL)     #
#                                                                             #
#                        http://www.gnu.org/licenses/                         #
###############################################################################

#include "fmpz_poly.pxi"
#include "fmpz.pxi"

cdef extern from "gmp.h":
    ctypedef void * mpz_t
    ctypedef void * mpq_t

cdef extern from "fmpq_poly.h":
    ctypedef void * fmpz_t
    ctypedef void * fmpz_poly_p
    struct fmpq_poly:
        fmpz_poly_p num
        fmpz_t den

    ctypedef fmpq_poly fmpq_poly_struct
    ctypedef fmpq_poly_struct fmpq_poly_t[1]

    void * fmpq_poly_canonicalize(fmpq_poly_t, fmpz_t)

    void * fmpq_poly_numref(fmpq_poly_t)
    void * fmpq_poly_denref(fmpq_poly_t)

    void fmpq_poly_init(fmpq_poly_t)
    void fmpq_poly_clear(fmpq_poly_t)

    void fmpq_poly_set(fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_set_fmpz_poly(fmpq_poly_t, fmpz_poly_t)
    void fmpq_poly_set_si(fmpq_poly_t, long)
    void fmpq_poly_set_mpz(fmpq_poly_t, mpz_t)
    void fmpq_poly_set_mpq(fmpq_poly_t, mpq_t)
    void fmpq_poly_swap(fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_zero(fmpq_poly_t)
    void fmpq_poly_neg(fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_get_coeff_mpq(mpq_t, fmpq_poly_t, unsigned long)
    void fmpq_poly_set_coeff_si(fmpq_poly_t, unsigned long, long)
    void fmpq_poly_set_coeff_mpq(fmpq_poly_t, unsigned long, mpq_t)
    void fmpq_poly_set_coeff_mpz(fmpq_poly_t, unsigned long, mpz_t)

    int fmpq_poly_equal(fmpq_poly_t, fmpq_poly_t)
    int fmpq_poly_cmp(fmpq_poly_t, fmpq_poly_t)
    int fmpq_poly_is_zero(fmpq_poly_t)

    long fmpq_poly_degree(fmpq_poly_t)
    unsigned long fmpq_poly_length(fmpq_poly_t)

    void fmpq_poly_add(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_sub(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_scalar_mul_mpq(fmpq_poly_t, fmpq_poly_t, mpq_t)
    void fmpq_poly_scalar_div_mpq(fmpq_poly_t, fmpq_poly_t, mpq_t)

    void fmpq_poly_mul(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_divrem(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_floordiv(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_mod(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_power(fmpq_poly_t, fmpq_poly_t, unsigned long)

    void fmpq_poly_left_shift(fmpq_poly_t, fmpq_poly_t, unsigned long)
    void fmpq_poly_right_shift(fmpq_poly_t, fmpq_poly_t, unsigned long)

    void fmpq_poly_gcd(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_xgcd(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_lcm(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_derivative(fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_evaluate_mpz(mpq_t, fmpq_poly_t, mpz_t)
    void fmpq_poly_evaluate_mpq(mpq_t, fmpq_poly_t, mpq_t)

    void fmpq_poly_content(mpq_t, fmpq_poly_t)
    void fmpq_poly_primitive_part(fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_resultant(mpq_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_compose(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_getslice(fmpq_poly_t, fmpq_poly_t, unsigned long, unsigned long)
    void fmpq_poly_truncate(fmpq_poly_t, fmpq_poly_t, unsigned long)
    void fmpq_poly_reverse(fmpq_poly_t, fmpq_poly_t, unsigned long)

    void _fmpq_poly_from_list(fmpq_poly_t, mpq_t *, unsigned long)
    void fmpq_poly_from_string(fmpq_poly_t, char *)
    char * fmpq_poly_to_string(fmpq_poly_t, char *)
    char * fmpq_poly_to_string_pretty(fmpq_poly_t, char *)

