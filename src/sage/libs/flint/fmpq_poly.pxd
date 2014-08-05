###############################################################################
#          Copyright (C) 2010 Sebastian Pancratz <sfp@pancratz.org>           #
#                                                                             #
#     Distributed under the terms of the GNU General Public License (GPL)     #
#                                                                             #
#                        http://www.gnu.org/licenses/                         #
###############################################################################

cdef extern from "gmp.h":
    ctypedef void * mpz_t
    ctypedef void * mpq_t

cdef extern from "flint/fmpq.h":
    ctypedef void * fmpq_t
    void fmpq_init(fmpq_t)
    void fmpq_clear(fmpq_t)
    void fmpq_get_mpq(mpq_t, fmpq_t)
    void fmpq_set_mpq(fmpq_t, mpq_t)

cdef extern from "flint/fmpz_vec.h":
    long _fmpz_vec_max_limbs(void * c, long n)

cdef extern from "flint/fmpq_poly.h":
    ctypedef void * fmpz_t
    ctypedef void * fmpq_poly_t

    void fmpq_poly_canonicalise(fmpq_poly_t)

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

    void fmpq_poly_get_numerator(fmpz_poly_t, fmpq_poly_t)

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
    void fmpq_poly_div(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_rem(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_pow(fmpq_poly_t, fmpq_poly_t, unsigned long)

    void fmpq_poly_shift_left(fmpq_poly_t, fmpq_poly_t, unsigned long)
    void fmpq_poly_shift_right(fmpq_poly_t, fmpq_poly_t, unsigned long)

    void fmpq_poly_gcd(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)
    void fmpq_poly_xgcd(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t, fmpq_poly_t, \
            fmpq_poly_t)
    void fmpq_poly_lcm(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_derivative(fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_evaluate_mpz(mpq_t, fmpq_poly_t, mpz_t)
    void fmpq_poly_evaluate_mpq(mpq_t, fmpq_poly_t, mpq_t)

    void fmpq_poly_content(fmpq_t, fmpq_poly_t)
    void fmpq_poly_primitive_part(fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_resultant(fmpq_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_compose(fmpq_poly_t, fmpq_poly_t, fmpq_poly_t)

    void fmpq_poly_get_slice(fmpq_poly_t, fmpq_poly_t, long, long)
    void fmpq_poly_truncate(fmpq_poly_t, unsigned long)
    void fmpq_poly_reverse(fmpq_poly_t, fmpq_poly_t, unsigned long)
    void fmpq_poly_revert_series(fmpq_poly_t, fmpq_poly_t, unsigned long)

    void fmpq_poly_set_array_mpq(fmpq_poly_t, mpq_t *, unsigned long)
    void fmpq_poly_from_string(fmpq_poly_t, char *)
    char * fmpq_poly_to_string(fmpq_poly_t, char *)
    char * fmpq_poly_to_string_pretty(fmpq_poly_t, char *)
