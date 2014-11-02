include "fmpz.pxi"
include "sage/libs/ntl/decl.pxi"

from sage.libs.flint.nmod_poly cimport nmod_poly_t, mp_srcptr, nmod_poly_factor_t

cdef extern from "flint/fmpz_poly.h":
    ctypedef struct fmpz_poly_struct:
        fmpz * coeffs
        long alloc
        long length

    ctypedef fmpz_poly_struct[1] fmpz_poly_t

    ctypedef struct fmpz_poly_powers_precomp_struct:
        fmpz ** powers
        long len

    ctypedef fmpz_poly_powers_precomp_struct fmpz_poly_powers_precomp_t[1]

    ctypedef struct fmpz_poly_factor_struct:
        fmpz c
        fmpz_poly_struct *p
        long *exp
        long num
        long alloc

    ctypedef fmpz_poly_factor_struct fmpz_poly_factor_t[1]

    void fmpz_poly_init(fmpz_poly_t poly)

    void fmpz_poly_init2(fmpz_poly_t poly, long alloc)

    void fmpz_poly_realloc(fmpz_poly_t poly, long alloc)

    void fmpz_poly_fit_length(fmpz_poly_t poly, long len)

    void fmpz_poly_clear(fmpz_poly_t poly)

    void _fmpz_poly_normalise(fmpz_poly_t poly)

    void _fmpz_poly_set_length(fmpz_poly_t poly, long newlen)

    long fmpz_poly_length(fmpz_poly_t poly)

    long fmpz_poly_degree(fmpz_poly_t poly)

    #  Assignment and basic manipulation

    void fmpz_poly_set(fmpz_poly_t poly1, fmpz_poly_t poly2)

    void fmpz_poly_set_ui(fmpz_poly_t poly, unsigned long c)

    void fmpz_poly_set_si(fmpz_poly_t poly, long c)

    void fmpz_poly_set_fmpz(fmpz_poly_t poly, fmpz_t c)

    void fmpz_poly_set_mpz(fmpz_poly_t poly, mpz_t c)

    int _fmpz_poly_set_str(fmpz * poly, char * str)

    int fmpz_poly_set_str(fmpz_poly_t poly, char * str)

    char * _fmpz_poly_get_str(fmpz * poly, long len)

    char * fmpz_poly_get_str(fmpz_poly_t poly)

    char * _fmpz_poly_get_str_pretty(fmpz * poly, long len, char * x)

    char * fmpz_poly_get_str_pretty(fmpz_poly_t poly, char * x)

    void fmpz_poly_zero(fmpz_poly_t poly)

    void fmpz_poly_one(fmpz_poly_t poly)

    void fmpz_poly_zero_coeffs(fmpz_poly_t poly, long i, long j)

    void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_reverse(fmpz * res, fmpz * poly, long len, long n)

    void fmpz_poly_reverse(fmpz_poly_t res, fmpz_poly_t poly, long n)

    void fmpz_poly_truncate(fmpz_poly_t poly, long newlen)

    #  Getting and setting coefficients

    long fmpz_poly_get_coeff_si(fmpz_poly_t poly, long n)

    void fmpz_poly_set_coeff_si(fmpz_poly_t poly, long n, long x)

    unsigned long fmpz_poly_get_coeff_ui(fmpz_poly_t poly, long n)

    void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, long n, unsigned long x)

    void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, long n, fmpz_t x)

    void fmpz_poly_get_coeff_fmpz(fmpz_t x, fmpz_poly_t poly, long n)

    fmpz *fmpz_poly_get_coeff_ptr(fmpz_poly_t poly, long n)

    fmpz *fmpz_poly_lead(fmpz_poly_t poly)

    #  Comparison

    int fmpz_poly_equal(fmpz_poly_t poly1, fmpz_poly_t poly2)

    int fmpz_poly_is_zero(fmpz_poly_t poly)

    int _fmpz_poly_is_one(fmpz *poly, long len)

    int fmpz_poly_is_one(fmpz_poly_t op)

    int fmpz_poly_is_unit(fmpz_poly_t op)

    int fmpz_poly_equal_fmpz(fmpz_poly_t poly, fmpz_t c)

    # Addition and subtraction

    void _fmpz_poly_add(fmpz * res, fmpz * poly1, long len1,
                        fmpz * poly2, long len2)

    void fmpz_poly_add(fmpz_poly_t res, fmpz_poly_t poly1,
                       fmpz_poly_t poly2)

    void _fmpz_poly_sub(fmpz * res, fmpz * poly1, long len1,
                        fmpz * poly2, long len2)

    void fmpz_poly_sub(fmpz_poly_t res, fmpz_poly_t poly1,
                       fmpz_poly_t poly2)

    void fmpz_poly_neg(fmpz_poly_t res, fmpz_poly_t poly)

    #  Scalar multiplication and division

    void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1,
                                 fmpz_poly_t poly2, unsigned long x)

    void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1,
                                 fmpz_poly_t poly2, long x)

    void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1,
                                   fmpz_poly_t poly2, fmpz_t x)

    void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1,
                                      fmpz_poly_t poly2, fmpz_t x)

    void fmpz_poly_scalar_submul_fmpz(fmpz_poly_t poly1,
                                      fmpz_poly_t poly2, fmpz_t x)

    void fmpz_poly_scalar_fdiv_ui(fmpz_poly_t poly1,
                                  fmpz_poly_t poly2, unsigned long x)

    void fmpz_poly_scalar_fdiv_si(fmpz_poly_t poly1,
                                  fmpz_poly_t poly2, long x)

    void fmpz_poly_scalar_fdiv_fmpz(fmpz_poly_t poly1,
                                    fmpz_poly_t poly2, fmpz_t x)

    void fmpz_poly_scalar_tdiv_ui(fmpz_poly_t poly1,
                                  fmpz_poly_t poly2, unsigned long x)

    void fmpz_poly_scalar_tdiv_si(fmpz_poly_t poly1,
                                  fmpz_poly_t poly2, long x)

    void fmpz_poly_scalar_tdiv_fmpz(fmpz_poly_t poly1,
                                    fmpz_poly_t poly2, fmpz_t x)

    void fmpz_poly_scalar_divexact_ui(fmpz_poly_t poly1,
                                      fmpz_poly_t poly2, unsigned long x)

    void fmpz_poly_scalar_divexact_si(fmpz_poly_t poly1,
                                      fmpz_poly_t poly2, long x)

    void fmpz_poly_scalar_divexact_fmpz(fmpz_poly_t poly1,
                                        fmpz_poly_t poly2, fmpz_t x)

    void fmpz_poly_scalar_fdiv_2exp(fmpz_poly_t poly1, fmpz_poly_t poly2,
                                    unsigned long exp)

    void fmpz_poly_scalar_tdiv_2exp(fmpz_poly_t poly1, fmpz_poly_t poly2,
                                    unsigned long exp)

    void fmpz_poly_scalar_mul_2exp(fmpz_poly_t poly1, fmpz_poly_t poly2,
                                   unsigned long exp)

    void fmpz_poly_scalar_mod_fmpz(fmpz_poly_t poly1,
                                   fmpz_poly_t poly2, fmpz_t x)

    void fmpz_poly_scalar_smod_fmpz(fmpz_poly_t poly1,
                                    fmpz_poly_t poly2, fmpz_t x)

    #  Multiplication

    void _fmpz_poly_mul_classical(fmpz * res, fmpz * poly1, long len1,
                                  fmpz * poly2, long len2)

    void fmpz_poly_mul_classical(fmpz_poly_t res,
                                 fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_mullow_classical(fmpz * res, fmpz * poly1, long len1,
                                     fmpz * poly2, long len2, long n)

    void fmpz_poly_mullow_classical(fmpz_poly_t res, fmpz_poly_t poly1,
                                    fmpz_poly_t poly2, long n)

    void _fmpz_poly_mulhigh_classical(fmpz * res, fmpz * poly1,
                                      long len1, fmpz * poly2, long len2, long start)

    void fmpz_poly_mulhigh_classical(fmpz_poly_t res,
                                     fmpz_poly_t poly1, fmpz_poly_t poly2, long start)

    void _fmpz_poly_mulmid_classical(fmpz * res, fmpz * poly1,
                                     long len1, fmpz * poly2, long len2)

    void fmpz_poly_mulmid_classical(fmpz_poly_t res,
                                    fmpz_poly_t poly1, fmpz_poly_t poly2)

    void fmpz_poly_mul_karatsuba(fmpz_poly_t res,
                                 fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_mul_karatsuba(fmpz * res, fmpz * poly1,
                                  long len1, fmpz * poly2, long len2)

    void _fmpz_poly_mullow_karatsuba_n(fmpz * res, fmpz * poly1,
                                       fmpz * poly2, long n)

    void fmpz_poly_mullow_karatsuba_n(fmpz_poly_t res,
                                      fmpz_poly_t poly1, fmpz_poly_t poly2, long n)

    void _fmpz_poly_mulhigh_karatsuba_n(fmpz * res, fmpz * poly1,
                                        fmpz * poly2, long len)

    void fmpz_poly_mulhigh_karatsuba_n(fmpz_poly_t res,
                                       fmpz_poly_t poly1, fmpz_poly_t poly2, long length)

    void _fmpz_poly_mul_KS(fmpz * res, fmpz * poly1, long len1,
                           fmpz * poly2, long len2)

    void fmpz_poly_mul_KS(fmpz_poly_t res,
                          fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_mullow_KS(fmpz * res, fmpz * poly1, long len1,
                              fmpz * poly2, long len2, long n)

    void fmpz_poly_mullow_KS(fmpz_poly_t res, fmpz_poly_t poly1,
                             fmpz_poly_t poly2, long n)

    void _fmpz_poly_mul_SS(fmpz * output, fmpz * input1, long length1,
                           fmpz * input2, long length2)

    void fmpz_poly_mul_SS(fmpz_poly_t res,
                          fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_mullow_SS(fmpz * output, fmpz * input1, long length1,
                              fmpz * input2, long length2, long n)

    void fmpz_poly_mullow_SS(fmpz_poly_t res,
                             fmpz_poly_t poly1, fmpz_poly_t poly2, long n)

    void _fmpz_poly_mul(fmpz * res, fmpz * poly1,
                        long len1, fmpz * poly2, long len2)

    void fmpz_poly_mul(fmpz_poly_t res,
                       fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_mullow(fmpz * res, fmpz * poly1, long len1,
                           fmpz * poly2, long len2, long n)

    void fmpz_poly_mullow(fmpz_poly_t res,
                          fmpz_poly_t poly1, fmpz_poly_t poly2, long n)

    void fmpz_poly_mulhigh_n(fmpz_poly_t res,
                             fmpz_poly_t poly1, fmpz_poly_t poly2, long n)

    # Squaring

    void _fmpz_poly_sqr_KS(fmpz * rop, fmpz * op, long len)

    void fmpz_poly_sqr_KS(fmpz_poly_t rop, fmpz_poly_t op)

    void fmpz_poly_sqr_karatsuba(fmpz_poly_t rop, fmpz_poly_t op)

    void _fmpz_poly_sqr_karatsuba(fmpz * rop, fmpz * op, long len)

    void _fmpz_poly_sqr_classical(fmpz * rop, fmpz * op, long len)

    void fmpz_poly_sqr_classical(fmpz_poly_t rop, fmpz_poly_t op)

    void _fmpz_poly_sqr(fmpz * rop, fmpz * op, long len)

    void fmpz_poly_sqr(fmpz_poly_t rop, fmpz_poly_t op)

    void _fmpz_poly_sqrlow_KS(fmpz * res, fmpz * poly, long len, long n)

    void fmpz_poly_sqrlow_KS(fmpz_poly_t res, fmpz_poly_t poly, long n)

    void _fmpz_poly_sqrlow_karatsuba_n(fmpz * res, fmpz * poly, long n)

    void fmpz_poly_sqrlow_karatsuba_n(fmpz_poly_t res, fmpz_poly_t poly, long n)

    void _fmpz_poly_sqrlow_classical(fmpz * res, fmpz * poly, long len, long n)

    void fmpz_poly_sqrlow_classical(fmpz_poly_t res, fmpz_poly_t poly, long n)

    void _fmpz_poly_sqrlow(fmpz * res, fmpz * poly, long len, long n)

    void fmpz_poly_sqrlow(fmpz_poly_t res, fmpz_poly_t poly, long n)

    #  Powering

    void _fmpz_poly_pow_multinomial(fmpz * res, fmpz * poly, long len, unsigned long e)

    void fmpz_poly_pow_multinomial(fmpz_poly_t res, fmpz_poly_t poly, unsigned long e)

    void _fmpz_poly_pow_binomial(fmpz * res, fmpz * poly, unsigned long e)

    void fmpz_poly_pow_binomial(fmpz_poly_t res, fmpz_poly_t poly, unsigned long e)

    void _fmpz_poly_pow_binexp(fmpz * res, fmpz * poly, long len, unsigned long e)

    void fmpz_poly_pow_binexp(fmpz_poly_t res, fmpz_poly_t poly, unsigned long e)

    void _fmpz_poly_pow_addchains(fmpz * res, fmpz * poly, long len, int * a, int n)

    void fmpz_poly_pow_addchains(fmpz_poly_t res, fmpz_poly_t poly, unsigned long e)

    void _fmpz_poly_pow_small(fmpz * res, fmpz * poly, long len, unsigned long e)

    void _fmpz_poly_pow(fmpz * res, fmpz * poly, long len, unsigned long e)

    void fmpz_poly_pow(fmpz_poly_t res, fmpz_poly_t poly, unsigned long e)

    void _fmpz_poly_pow_trunc(fmpz * res, fmpz * poly, unsigned long e, long n)

    void fmpz_poly_pow_trunc(fmpz_poly_t res, fmpz_poly_t poly, unsigned long e, long n)

    #  Shifting

    void _fmpz_poly_shift_left(fmpz * res, fmpz * poly, long len, long n)

    void _fmpz_poly_shift_right(fmpz * res, fmpz * poly, long len, long n)

    void fmpz_poly_shift_left(fmpz_poly_t res, fmpz_poly_t poly, long n)

    void fmpz_poly_shift_right(fmpz_poly_t res, fmpz_poly_t poly, long n)

    #  Norms

    void _fmpz_poly_2norm(fmpz_t res, fmpz * poly, long len)

    void fmpz_poly_2norm(fmpz_t res, fmpz_poly_t poly)

    mp_bitcnt_t _fmpz_poly_2norm_normalised_bits(fmpz * poly, long len)

    unsigned long fmpz_poly_max_limbs(fmpz_poly_t poly)

    long fmpz_poly_max_bits(fmpz_poly_t poly)

    void fmpz_poly_height(fmpz_t res, fmpz_poly_t poly)

    #  Greatest common divisor

    void _fmpz_poly_gcd_subresultant(fmpz * res, fmpz * poly1, long len1,
                                     fmpz * poly2, long len2)

    void fmpz_poly_gcd_subresultant(fmpz_poly_t res, fmpz_poly_t poly1,
                                    fmpz_poly_t poly2)

    int _fmpz_poly_gcd_heuristic(fmpz * res, fmpz * poly1, long len1,
                                 fmpz * poly2, long len2)

    int fmpz_poly_gcd_heuristic(fmpz_poly_t res,
                                fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_gcd_modular(fmpz * res, fmpz * poly1, long len1,
                                fmpz * poly2, long len2)

    void fmpz_poly_gcd_modular(fmpz_poly_t res,
                               fmpz_poly_t poly1, fmpz_poly_t poly2)

    void _fmpz_poly_gcd(fmpz * res, fmpz * poly1, long len1,
                        fmpz * poly2, long len2)

    void fmpz_poly_gcd(fmpz_poly_t res, fmpz_poly_t poly1,
                       fmpz_poly_t poly2)

    void _fmpz_poly_lcm(fmpz * res, fmpz * poly1, long len1,
                        fmpz * poly2, long len2)

    void fmpz_poly_lcm(fmpz_poly_t res, fmpz_poly_t poly1,
                       fmpz_poly_t poly2)

    void _fmpz_poly_resultant(fmpz_t res, fmpz * poly1, long len1,
                              fmpz * poly2, long len2)

    void fmpz_poly_resultant(fmpz_t res, fmpz_poly_t poly1,
                             fmpz_poly_t poly2)

    void _fmpz_poly_xgcd_modular(fmpz_t r, fmpz * s, fmpz * t,
                                 fmpz * poly1, long len1, fmpz * poly2, long len2)

    void fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                                fmpz_poly_t poly1, fmpz_poly_t poly2)


    void _fmpz_poly_xgcd(fmpz_t r, fmpz * s, fmpz * t,
                         fmpz * poly1, long len1, fmpz * poly2, long len2)

    void fmpz_poly_xgcd(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                        fmpz_poly_t poly1, fmpz_poly_t poly2)

    #  Gaussian content

    void _fmpz_poly_content(fmpz_t res, fmpz * poly, long len)

    void fmpz_poly_content(fmpz_t res, fmpz_poly_t poly)

    void _fmpz_poly_primitive_part(fmpz * res, fmpz * poly, long len)

    void fmpz_poly_primitive_part(fmpz_poly_t res, fmpz_poly_t poly)

    #  Square-free

    int _fmpz_poly_is_squarefree(fmpz * poly, long len)

    int fmpz_poly_is_squarefree(fmpz_poly_t poly)

    #  Euclidean division

    void _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, fmpz * A,
                                       long lenA, fmpz * B, long lenB)

    void fmpz_poly_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R,
                                   fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W,
                                                fmpz * A, fmpz * B, long lenB)

    void _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R,
                                      fmpz * A, long lenA, fmpz * B, long lenB)

    void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R,
                                     fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_divrem(fmpz * Q, fmpz * R, fmpz * A, long lenA,
                           fmpz * B, long lenB)

    void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A,
                          fmpz_poly_t B)

    void _fmpz_poly_div_basecase(fmpz * Q, fmpz * R, fmpz * A, long lenA,
                                 fmpz * B, long lenB)

    void fmpz_poly_div_basecase(fmpz_poly_t Q,
                                fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_divremlow_divconquer_recursive(fmpz * Q, fmpz * QB,
                                                   fmpz * A, fmpz * B, long lenB)

    void _fmpz_poly_div_divconquer_recursive(fmpz * Q, fmpz * temp,
                                             fmpz * A, fmpz * B, long lenB)

    void _fmpz_poly_div_divconquer(fmpz * Q, fmpz * A, long lenA,
                                   fmpz * B, long lenB)

    void fmpz_poly_div_divconquer(fmpz_poly_t Q,
                                  fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_div(fmpz * Q, fmpz * A, long lenA,
                        fmpz * B, long lenB)

    void fmpz_poly_div(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_preinvert(fmpz * B_inv, fmpz * B, long n)

    void fmpz_poly_preinvert(fmpz_poly_t B_inv, fmpz_poly_t B)

    void _fmpz_poly_div_preinv(fmpz * Q, fmpz * A, long len1,
                               fmpz * B, fmpz * B_inv, long len2)

    void fmpz_poly_div_preinv(fmpz_poly_t Q, fmpz_poly_t A,
                              fmpz_poly_t B, fmpz_poly_t B_inv)

    void _fmpz_poly_divrem_preinv(fmpz * Q, fmpz * A, long len1,
                                  fmpz * B, fmpz * B_inv, long len2)

    void fmpz_poly_divrem_preinv(fmpz_poly_t Q, fmpz_poly_t R,
                                 fmpz_poly_t A, fmpz_poly_t B, fmpz_poly_t B_inv)

    fmpz ** _fmpz_poly_powers_precompute(fmpz * B, long len)

    void fmpz_poly_powers_precompute(fmpz_poly_powers_precomp_t pinv,
                                     fmpz_poly_t poly)

    void _fmpz_poly_powers_clear(fmpz ** powers, long len)

    void fmpz_poly_powers_clear(fmpz_poly_powers_precomp_t pinv)

    void _fmpz_poly_rem_powers_precomp(fmpz * A, long m,
                                       fmpz * B, long n, fmpz ** powers)

    void fmpz_poly_rem_powers_precomp(fmpz_poly_t R,
                                      fmpz_poly_t A, fmpz_poly_t B,
                                      fmpz_poly_powers_precomp_t B_inv)

    void _fmpz_poly_rem_basecase(fmpz * Q, fmpz * A, long lenA,
                                 fmpz * B, long lenB)

    void fmpz_poly_rem_basecase(fmpz_poly_t R,
                                fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_rem(fmpz * R, fmpz * A, long lenA,
                        fmpz * B, long lenB)

    void fmpz_poly_rem(fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)

    void fmpz_poly_div_root(fmpz_poly_t Q, fmpz_poly_t A, fmpz_t c)

    void _fmpz_poly_div_root(fmpz * Q, fmpz * A, long len, fmpz_t c)

    #  Power series division

    void _fmpz_poly_inv_series_newton(fmpz * Qinv, fmpz * Q, long n)

    void fmpz_poly_inv_series_newton(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)

    void  _fmpz_poly_inv_series(fmpz * Qinv, fmpz * Q, long n)

    void  fmpz_poly_inv_series(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)

    void _fmpz_poly_div_series(fmpz * Q, fmpz * A, fmpz * B, long n)

    void fmpz_poly_div_series(fmpz_poly_t Q, fmpz_poly_t A,
                              fmpz_poly_t B, long n)

    #  Divisibility testing

    int _fmpz_poly_divides(fmpz * q, fmpz * a,
                           long len1, fmpz * b, long len2)

    int fmpz_poly_divides(fmpz_poly_t q, fmpz_poly_t a, fmpz_poly_t b)


    #  Pseudo division

    void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R,
                                           unsigned long * d, fmpz * A, long A_len,
                                           fmpz * B, long B_len, fmpz_preinvn_t inv)

    void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R,
                                          unsigned long * d, fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_pseudo_divrem_divconquer(fmpz * Q, fmpz * R,
                                             unsigned long * d, fmpz * A, long lenA,
                                             fmpz * B, long lenB, fmpz_preinvn_t inv)

    void fmpz_poly_pseudo_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R,
                                            unsigned long * d, fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_pseudo_divrem_cohen(fmpz * Q, fmpz * R, fmpz * A,
                                        long lenA, fmpz * B, long lenB)

    void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R,
                                       fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_pseudo_rem_cohen(fmpz * R, fmpz * A, long lenA,
                                     fmpz * B, long lenB)

    void fmpz_poly_pseudo_rem_cohen(fmpz_poly_t R, fmpz_poly_t A,
                                    fmpz_poly_t B)


    void _fmpz_poly_pseudo_divrem(fmpz * Q, fmpz * R,
                                  unsigned long * d, fmpz * A, long A_len,
                                  fmpz * B, long B_len, fmpz_preinvn_t inv)

    void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R,
                                 unsigned long * d, fmpz_poly_t A, fmpz_poly_t B)

    void _fmpz_poly_pseudo_div(fmpz * Q, unsigned long * d, fmpz * A, long lenA,
                               fmpz * B, long lenB, fmpz_preinvn_t inv)

    void fmpz_poly_pseudo_div(fmpz_poly_t Q, unsigned long * d, fmpz_poly_t A,
                              fmpz_poly_t B)

    void _fmpz_poly_pseudo_rem(fmpz * R, unsigned long * d, fmpz * A, long lenA,
                               fmpz * B, long lenB, fmpz_preinvn_t inv)

    void fmpz_poly_pseudo_rem(fmpz_poly_t R, unsigned long * d, fmpz_poly_t A,
                              fmpz_poly_t B)

    #  Derivative

    void _fmpz_poly_derivative(fmpz * rpoly, fmpz * poly, long len)

    void fmpz_poly_derivative(fmpz_poly_t res, fmpz_poly_t poly)

    #  Evaluation

    void  _fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, fmpz * poly, long len,
                                              fmpz_t a)

    void fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, fmpz_poly_t poly,
                                            fmpz_t a)

    void _fmpz_poly_evaluate_horner_fmpz(fmpz_t res, fmpz * f, long len,
                                         fmpz_t a)

    void fmpz_poly_evaluate_horner_fmpz(fmpz_t res, fmpz_poly_t f,
                                        fmpz_t a)

    void _fmpz_poly_evaluate_fmpz(fmpz_t res, fmpz * f, long len, fmpz_t a)

    void fmpz_poly_evaluate_fmpz(fmpz_t res, fmpz_poly_t f, fmpz_t a)

    void _fmpz_poly_evaluate_horner_mpq(fmpz_t rnum, fmpz_t rden,
                                        fmpz * f, long len,
                                        fmpz_t anum, fmpz_t aden)

    void fmpz_poly_evaluate_horner_mpq(mpq_t res, fmpz_poly_t f, mpq_t a)

    void _fmpz_poly_evaluate_mpq(fmpz_t rnum, fmpz_t rden,
                                 fmpz * f, long len,
                                 fmpz_t anum, fmpz_t aden)

    void fmpz_poly_evaluate_mpq(mpq_t res, fmpz_poly_t f, mpq_t a)

    mp_limb_t _fmpz_poly_evaluate_mod(fmpz * poly, long len, mp_limb_t a,
                                      mp_limb_t n, mp_limb_t ninv)

    mp_limb_t fmpz_poly_evaluate_mod(fmpz_poly_t poly, mp_limb_t a,
                                     mp_limb_t n)

    void  _fmpz_poly_evaluate_divconquer(fmpz * res, fmpz * poly, long len,
                                         fmpz_t x)

    void  fmpz_poly_evaluate_divconquer(fmpz_t res,
                                        fmpz_poly_t poly, fmpz_t x)

    #  Composition

    void _fmpz_poly_compose_horner(fmpz * res, fmpz * poly1, long len1,
                                   fmpz * poly2, long len2)

    void fmpz_poly_compose_horner(fmpz_poly_t res, fmpz_poly_t poly1,
                                  fmpz_poly_t poly2)

    void _fmpz_poly_compose_divconquer(fmpz * res, fmpz * poly1, long len1,
                                       fmpz * poly2, long len2)

    void fmpz_poly_compose_divconquer(fmpz_poly_t res, fmpz_poly_t poly1,
                                      fmpz_poly_t poly2)

    void _fmpz_poly_compose(fmpz * res, fmpz * poly1, long len1,
                            fmpz * poly2, long len2)

    void fmpz_poly_compose(fmpz_poly_t res, fmpz_poly_t poly1,
                           fmpz_poly_t poly2)

    #  Taylor shift

    void _fmpz_poly_taylor_shift_horner(fmpz * poly, fmpz_t c, long n)

    void fmpz_poly_taylor_shift_horner(fmpz_poly_t g, fmpz_poly_t f,
                                       fmpz_t c)

    void _fmpz_poly_taylor_shift_divconquer(fmpz * poly, fmpz_t c, long n)

    void fmpz_poly_taylor_shift_divconquer(fmpz_poly_t g, fmpz_poly_t f,
                                           fmpz_t c)

    void _fmpz_poly_taylor_shift(fmpz * poly, fmpz_t c, long n)

    void fmpz_poly_taylor_shift(fmpz_poly_t g, fmpz_poly_t f, fmpz_t c)

    #  Power series composition and compositional inverse

    void _fmpz_poly_compose_series_brent_kung(fmpz * res, fmpz * poly1, long len1,
                                              fmpz * poly2, long len2, long n)

    void fmpz_poly_compose_series_brent_kung(fmpz_poly_t res,
                                             fmpz_poly_t poly1, fmpz_poly_t poly2, long n)

    void _fmpz_poly_compose_series_horner(fmpz * res, fmpz * poly1, long len1,
                                          fmpz * poly2, long len2, long n)

    void fmpz_poly_compose_series_horner(fmpz_poly_t res,
                                         fmpz_poly_t poly1, fmpz_poly_t poly2, long n)

    void _fmpz_poly_compose_series(fmpz * res, fmpz * poly1, long len1,
                                   fmpz * poly2, long len2, long n)

    void fmpz_poly_compose_series(fmpz_poly_t res,
                              fmpz_poly_t poly1, fmpz_poly_t poly2, long n)

    void _fmpz_poly_revert_series_lagrange(fmpz * Qinv, fmpz * Q, long n)

    void fmpz_poly_revert_series_lagrange(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)

    void_fmpz_poly_revert_series_lagrange_fast(fmpz * Qinv, fmpz * Q, long n)

    void fmpz_poly_revert_series_lagrange_fast(fmpz_poly_t Qinv,
                                               fmpz_poly_t Q, long n)

    void _fmpz_poly_revert_series_newton(fmpz * Qinv, fmpz * Q, long n)

    void fmpz_poly_revert_series_newton(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)

    void _fmpz_poly_revert_series(fmpz * Qinv, fmpz * Q, long n)

    void fmpz_poly_revert_series(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)

    #  Square root

    int _fmpz_poly_sqrt_classical(fmpz * res, fmpz * poly, long len)

    int fmpz_poly_sqrt_classical(fmpz_poly_t b, fmpz_poly_t a)

    int _fmpz_poly_sqrt(fmpz * res, fmpz * poly, long len)

    int fmpz_poly_sqrt(fmpz_poly_t b, fmpz_poly_t a)

    #  Signature

    void _fmpz_poly_signature(long * r1, long * r2, fmpz * poly, long len)

    void fmpz_poly_signature(long * r1, long * r2, fmpz_poly_t poly)

    #  Input and output

    int fmpz_poly_fprint(FILE * file, fmpz_poly_t poly)

    int _fmpz_poly_fprint_pretty(FILE * file,
                                 fmpz * poly, long len, char * x)

    int fmpz_poly_fprint_pretty(FILE * file,
                                fmpz_poly_t poly, char * x)

    int fmpz_poly_print(fmpz_poly_t poly)

    int fmpz_poly_print_pretty(fmpz_poly_t poly, char * x)

    int fmpz_poly_fread(FILE * file, fmpz_poly_t poly)

    int fmpz_poly_fread_pretty(FILE *file, fmpz_poly_t poly, char **x)

    int fmpz_poly_read(fmpz_poly_t poly)

    int fmpz_poly_read_pretty(fmpz_poly_t poly, char **x)

    void fmpz_poly_debug(fmpz_poly_t poly)

    #  CRT

    void fmpz_poly_get_nmod_poly(nmod_poly_t res, fmpz_poly_t poly)

    void fmpz_poly_set_nmod_poly(fmpz_poly_t res, nmod_poly_t poly)

    void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, nmod_poly_t poly)

    void _fmpz_poly_CRT_ui_precomp(fmpz * res, fmpz * poly1, long len1,
                                   fmpz_t m1, mp_srcptr poly2, long len2, mp_limb_t m2,
                                   mp_limb_t m2inv, fmpz_t m1m2, mp_limb_t c, int sign)

    void _fmpz_poly_CRT_ui(fmpz * res, fmpz * poly1, long len1,
                           fmpz_t m1, mp_srcptr poly2, long len2, mp_limb_t m2,
                           mp_limb_t m2inv, int sign)

    void fmpz_poly_CRT_ui(fmpz_poly_t res, fmpz_poly_t poly1,
                          fmpz_t m1, nmod_poly_t poly2,
                          int sign)


    # Products

    void _fmpz_poly_product_roots_fmpz_vec(fmpz * poly,
                                           fmpz * xs, long n)

    void fmpz_poly_product_roots_fmpz_vec(fmpz_poly_t poly,
                                          fmpz * xs, long n)

    # Newton basis

    void _fmpz_poly_monomial_to_newton(fmpz * poly, fmpz * roots, long n)

    void _fmpz_poly_newton_to_monomial(fmpz * poly, fmpz * roots, long n)


    # Multipoint evaluation and interpolation

    void fmpz_poly_evaluate_fmpz_vec(fmpz * res, fmpz_poly_t f,
                            fmpz * a, long n)

    void fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                        fmpz * xs, fmpz * ys, long n)

    # Hensel lifting

    void fmpz_poly_hensel_build_tree(long * link, fmpz_poly_t *v, fmpz_poly_t *w,
                                     nmod_poly_factor_t fac)

    void fmpz_poly_hensel_lift(fmpz_poly_t Gout, fmpz_poly_t Hout,
                               fmpz_poly_t Aout, fmpz_poly_t Bout,
                               fmpz_poly_t f,
                               fmpz_poly_t g, fmpz_poly_t h,
                               fmpz_poly_t a, fmpz_poly_t b,
                               fmpz_t p, fmpz_t p1)

    void _fmpz_poly_hensel_lift_without_inverse(fmpz *G, fmpz *H,
                                                fmpz *f, long lenF,
                                                fmpz *g, long lenG, fmpz *h, long lenH,
                                                fmpz *a, long lenA, fmpz *b, long lenB,
                                                fmpz_t p, fmpz_t p1)

    void fmpz_poly_hensel_lift_without_inverse(fmpz_poly_t Gout, fmpz_poly_t Hout,
                                               fmpz_poly_t f, fmpz_poly_t g, fmpz_poly_t h,
                                               fmpz_poly_t a, fmpz_poly_t b,
                                               fmpz_t p, fmpz_t p1)

    void _fmpz_poly_hensel_lift_only_inverse(fmpz *A, fmpz *B,
                                             fmpz *G, long lenG, fmpz *H, long lenH,
                                             fmpz *a, long lenA, fmpz *b, long lenB,
                                             fmpz_t p, fmpz_t p1)

    void fmpz_poly_hensel_lift_only_inverse(fmpz_poly_t Aout, fmpz_poly_t Bout,
                                            fmpz_poly_t G, fmpz_poly_t H,
                                            fmpz_poly_t a, fmpz_poly_t b,
                                            fmpz_t p, fmpz_t p1)

    void fmpz_poly_hensel_lift_tree_recursive(long *link,
                                              fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, long j, long inv,
                                              fmpz_t p0, fmpz_t p1)

    void fmpz_poly_hensel_lift_tree(long *link, fmpz_poly_t *v, fmpz_poly_t *w,
                                    fmpz_poly_t f, long r, fmpz_t p, long e0, long e1, long inv)

    long _fmpz_poly_hensel_start_lift(fmpz_poly_factor_t lifted_fac, long *link,
                                      fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f,
                                      nmod_poly_factor_t local_fac, long target_exp)

    long _fmpz_poly_hensel_continue_lift(fmpz_poly_factor_t lifted_fac,
                                         long *link, fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f,
                                         long prev, long curr, long N, fmpz_t p)

    void fmpz_poly_hensel_lift_once(fmpz_poly_factor_t lifted_fac,
                                    fmpz_poly_t f,
                                    nmod_poly_factor_t local_fac, long N)


    # Deprecated

    void fmpz_poly_scalar_mul_mpz(fmpz_poly_t poly1, fmpz_poly_t poly2,  mpz_t x)

    void fmpz_poly_scalar_divexact_mpz(fmpz_poly_t poly1, fmpz_poly_t poly2,  mpz_t x)

    void fmpz_poly_scalar_fdiv_mpz(fmpz_poly_t poly1, fmpz_poly_t poly2,  mpz_t x)

    void fmpz_poly_set_coeff_mpz(fmpz_poly_t poly, long n, mpz_t x)

    void fmpz_poly_get_coeff_mpz(mpz_t x,  fmpz_poly_t poly, long n)
