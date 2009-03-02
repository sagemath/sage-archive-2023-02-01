include "fmpz.pxi"
include "../ntl/decl.pxi"

cdef extern from "FLINT/fmpz_poly.h":

    ctypedef void* fmpz_poly_t

    void fmpz_poly_init(fmpz_poly_t poly)
    void fmpz_poly_init2(fmpz_poly_t poly, unsigned long alloc, \
            unsigned long limbs)
    void fmpz_poly_realloc(fmpz_poly_t poly, unsigned long alloc)

    void fmpz_poly_fit_length(fmpz_poly_t poly, unsigned long alloc)
    void fmpz_poly_resize_limbs(fmpz_poly_t poly, unsigned long limbs)
    void fmpz_poly_fit_limbs(fmpz_poly_t poly, unsigned long limbs)
    unsigned long fmpz_poly_limbs(fmpz_poly_t poly)

    void fmpz_poly_clear(fmpz_poly_t poly)

    long fmpz_poly_degree(fmpz_poly_t poly)
    unsigned long fmpz_poly_length(fmpz_poly_t poly)


    void fmpz_poly_set_length(fmpz_poly_t poly, unsigned long length)
    void fmpz_poly_truncate(fmpz_poly_t poly, unsigned long length)

    void fmpz_poly_set(fmpz_poly_t result, fmpz_poly_t poly)
    void fmpz_poly_set_coeff_si(fmpz_poly_t poly, unsigned long n, long x)
    void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, unsigned long n, \
            unsigned long x)
    void fmpz_poly_set_coeff_mpz(fmpz_poly_t poly, unsigned long n, mpz_t x)
    void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, unsigned long n, fmpz_t x)

    void fmpz_poly_get_coeff_mpz(mpz_t x, fmpz_poly_t poly, unsigned long n)
    void fmpz_poly_get_coeff_mpz_read_only(mpz_t x, fmpz_poly_t poly, unsigned long n)

    void fmpz_poly_add(fmpz_poly_t output, fmpz_poly_t input1, \
            fmpz_poly_t input2)
    void fmpz_poly_sub(fmpz_poly_t output, fmpz_poly_t input1, \
            fmpz_poly_t input2)
    void fmpz_poly_neg(fmpz_poly_t output, fmpz_poly_t input1)

    void fmpz_poly_mul(fmpz_poly_t output, fmpz_poly_t input1, \
            fmpz_poly_t input2)
    void fmpz_poly_mul_trunc_n(fmpz_poly_t output, fmpz_poly_t input1, \
            fmpz_poly_t input2, unsigned long trunc)
    void fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, fmpz_poly_t input1, \
            fmpz_poly_t input2, unsigned long trunc)

#    void fmpz_poly_scalar_mul(fmpz_poly_t output, fmpz_poly_t input, fmpz_t x)
    void fmpz_poly_scalar_mul_ui(fmpz_poly_t output, fmpz_poly_t input, \
            unsigned long x)
    void fmpz_poly_scalar_mul_si(fmpz_poly_t output, fmpz_poly_t input, long x)

    void fmpz_poly_scalar_mul_mpz(fmpz_poly_t output, fmpz_poly_t poly,
            mpz_t x)

    void fmpz_poly_scalar_div_exact_ui(fmpz_poly_t output, fmpz_poly_t poly, \
            unsigned long x)
    void fmpz_poly_scalar_div_exact_si(fmpz_poly_t output, fmpz_poly_t poly, \
            long x)
    void fmpz_poly_scalar_div_exact_fmpz(fmpz_poly_t output, fmpz_poly_t poly, \
            fmpz_t x)

    void fmpz_poly_scalar_div_mpz( fmpz_poly_t output, fmpz_poly_t poly, mpz_t x)

    void fmpz_poly_div(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)
    void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, \
            fmpz_poly_t B)

    void fmpz_poly_pseudo_div(fmpz_poly_t Q, unsigned long *d, fmpz_poly_t A, \
            fmpz_poly_t B)
    void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R, \
            unsigned long *d, fmpz_poly_t A, fmpz_poly_t B)

    int fmpz_poly_equal(fmpz_poly_t poly1, fmpz_poly_t poly2)

    bint fmpz_poly_from_string(fmpz_poly_t poly, char* s)
    char* fmpz_poly_to_string(fmpz_poly_t poly)
    void fmpz_poly_print(fmpz_poly_t poly)
    bint fmpz_poly_read(fmpz_poly_t poly)

    void fmpz_poly_power(fmpz_poly_t output, fmpz_poly_t poly, \
            unsigned long exp)
    void fmpz_poly_power_trunc_n(fmpz_poly_t output, fmpz_poly_t poly, \
            unsigned long exp, unsigned long n)

    void fmpz_poly_content(fmpz_t c, fmpz_poly_t poly)
    void fmpz_poly_primitive_part(fmpz_poly_t prim, fmpz_poly_t poly)

    void fmpz_poly_gcd(fmpz_poly_t res, fmpz_poly_t poly1, \
            fmpz_poly_t poly2)

    void fmpz_poly_xgcd(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t, fmpz_poly_t a,\
            fmpz_poly_t b)

    unsigned long fmpz_poly_resultant_bound(fmpz_poly_t a, fmpz_poly_t b)
    void fmpz_poly_resultant(fmpz_t r, fmpz_poly_t a, fmpz_poly_t b)

    unsigned long fmpz_poly_max_limbs(fmpz_poly_t poly)
