cdef extern from "FLINT/fmpz_poly.h":

    ctypedef void* fmpz_poly_t

    void fmpz_poly_init(fmpz_poly_t poly)
    void fmpz_poly_init2(fmpz_poly_t poly, unsigned long alloc, unsigned long limbs)
    void fmpz_poly_realloc(fmpz_poly_t poly, unsigned long alloc)

    void fmpz_poly_fit_length(fmpz_poly_t poly, unsigned long alloc)
    void fmpz_poly_resize_limbs(fmpz_poly_t poly, unsigned long limbs)
    void fmpz_poly_fit_limbs(fmpz_poly_t poly, unsigned long limbs)

    void fmpz_poly_clear(fmpz_poly_t poly)

    long fmpz_poly_degree(fmpz_poly_t poly)
    unsigned long fmpz_poly_length(fmpz_poly_t poly)


    void fmpz_poly_set_length(fmpz_poly_t poly, unsigned long length)
    void fmpz_poly_truncate(fmpz_poly_t poly, unsigned long length)

    void fmpz_poly_set_coeff_si(fmpz_poly_t poly, unsigned long n, long x)
    void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, unsigned long n, unsigned long x)

    void fmpz_poly_get_coeff_mpz(mpz_t x, fmpz_poly_t poly, unsigned long n)

    void fmpz_poly_mul(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2)
    void fmpz_poly_mul_trunc_n(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2, unsigned long trunc)
    void fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2, unsigned long trunc)

    void fmpz_poly_div_naive(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)

    bint fmpz_poly_from_string(fmpz_poly_t poly, char* s)
    char* fmpz_poly_to_string(fmpz_poly_t poly)
    void fmpz_poly_print(fmpz_poly_t poly)
    bint fmpz_poly_read(fmpz_poly_t poly)

    void fmpz_poly_divrem_naive(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)
    void fmpz_poly_div_karatsuba_recursive(fmpz_poly_t Q, fmpz_poly_t DQ, fmpz_poly_t A, fmpz_poly_t B)
    void fmpz_poly_divrem_karatsuba(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)
    void fmpz_poly_div_karatsuba(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)
    void fmpz_poly_newton_invert_basecase(fmpz_poly_t Q_inv, fmpz_poly_t Q, unsigned long n)
    void fmpz_poly_newton_invert(fmpz_poly_t Q_inv, fmpz_poly_t Q, unsigned long n)
    void fmpz_poly_div_series(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B, unsigned long n)
    void fmpz_poly_div_newton(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)

    void fmpz_poly_power(fmpz_poly_t output, fmpz_poly_t poly, unsigned long exp)


