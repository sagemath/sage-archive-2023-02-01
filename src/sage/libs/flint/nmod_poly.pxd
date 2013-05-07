include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"

from sage.libs.flint.flint cimport *

from flint import *

cdef extern from "flint/nmod_poly.h":
    ctypedef struct nmod_poly_struct:
        unsigned long *coeffs
        unsigned long alloc
        unsigned long length
        unsigned long p
        double p_inv

    ctypedef nmod_poly_struct* nmod_poly_t

    ctypedef struct nmod_poly_factor_struct:
        nmod_poly_t *factors
        unsigned long *exponents
        unsigned long alloc
        unsigned long num_factors

    ctypedef nmod_poly_factor_struct* nmod_poly_factor_t

    cdef void nmod_poly_init(nmod_poly_t poly, unsigned long p)
    cdef void nmod_poly_init_preinv(nmod_poly_t poly, unsigned long p, double p_inv)
    cdef void nmod_poly_init2(nmod_poly_t poly, unsigned long p, unsigned long alloc)
    cdef void nmod_poly_init2_preinv(nmod_poly_t poly, unsigned long p, double p_inv, unsigned long alloc)
    cdef void nmod_poly_clear(nmod_poly_t poly)

    cdef void nmod_poly_realloc(nmod_poly_t poly, unsigned long alloc)
    # _bits_ only applies to newly allocated coefficients, not existing ones...

    # this non-inlined version REQUIRES that alloc > poly->alloc
    void __nmod_poly_fit_length(nmod_poly_t poly, unsigned long alloc)

    # this is arranged so that the initial comparison (very frequent) is inlined,
    # but the actual allocation (infrequent) is not
    cdef void nmod_poly_fit_length(nmod_poly_t poly, unsigned long alloc)

    # ------------------------------------------------------
    # Setting/retrieving coefficients

    cdef unsigned long nmod_poly_get_coeff_ui(nmod_poly_t poly, unsigned long n)

    cdef unsigned long _nmod_poly_get_coeff_ui(nmod_poly_t poly, unsigned long n)

    cdef void nmod_poly_set_coeff_ui(nmod_poly_t poly, unsigned long n, unsigned long c)

    cdef void _nmod_poly_set_coeff_ui(nmod_poly_t poly, unsigned long n, unsigned long c)

    # ------------------------------------------------------
    # String conversions and I/O

    cdef int nmod_poly_from_string(nmod_poly_t poly, char* s)
    cdef char* nmod_poly_to_string(nmod_poly_t poly)
    cdef void nmod_poly_print(nmod_poly_t poly)
    cdef void nmod_poly_fprint(nmod_poly_t poly, FILE* f)
    cdef int nmod_poly_read(nmod_poly_t poly)
    cdef int nmod_poly_fread(nmod_poly_t poly, FILE* f)

    # ------------------------------------------------------
    # Length and degree

    cdef void _nmod_poly_normalise(nmod_poly_t poly)
    cdef int _nmod_poly_normalised(nmod_poly_t poly)
    cdef void nmod_poly_truncate(nmod_poly_t poly, unsigned long length)

    cdef unsigned long nmod_poly_length(nmod_poly_t poly)

    cdef long nmod_poly_degree(nmod_poly_t poly)

    cdef unsigned long nmod_poly_modulus(nmod_poly_t poly)

    cdef double nmod_poly_precomputed_inverse(nmod_poly_t poly)

    # ------------------------------------------------------
    # Assignment

    cdef void _nmod_poly_set(nmod_poly_t res, nmod_poly_t poly)
    cdef void nmod_poly_set(nmod_poly_t res, nmod_poly_t poly)

    cdef void nmod_poly_zero(nmod_poly_t poly)

    cdef void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)

    #
    # Subpolynomials
    #

    cdef void _nmod_poly_attach(nmod_poly_t output, nmod_poly_t input)

    cdef void nmod_poly_attach(nmod_poly_t output, nmod_poly_t input)

    #
    # Attach input shifted right by n to output
    #

    cdef void _nmod_poly_attach_shift(nmod_poly_t output, nmod_poly_t input, unsigned long n)

    cdef void nmod_poly_attach_shift(nmod_poly_t output, nmod_poly_t input, unsigned long n)

    #
    # Attach input to first n coefficients of input
    #

    cdef void _nmod_poly_attach_truncate(nmod_poly_t output,  nmod_poly_t input, unsigned long n)

    cdef void nmod_poly_attach_truncate(nmod_poly_t output, nmod_poly_t input, unsigned long n)

    #
    # Comparison functions
    #

    cdef int nmod_poly_equal(nmod_poly_t poly1, nmod_poly_t poly2)

    cdef int nmod_poly_is_one(nmod_poly_t poly1)

    #
    # Reversal
    #

    cdef void _nmod_poly_reverse(nmod_poly_t output, nmod_poly_t input, unsigned long length)
    cdef void nmod_poly_reverse(nmod_poly_t output, nmod_poly_t input, unsigned long length)

    #
    # Monic polys
    #

    cdef void nmod_poly_make_monic(nmod_poly_t output, nmod_poly_t pol)

    #
    # Addition and subtraction
    #

    cdef void nmod_poly_add(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_add_without_mod(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_sub(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void _nmod_poly_sub(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_neg(nmod_poly_t res, nmod_poly_t poly)

    #
    # Shifting functions
    #

    cdef void nmod_poly_left_shift(nmod_poly_t res, nmod_poly_t poly, unsigned long k)
    cdef void nmod_poly_right_shift(nmod_poly_t res, nmod_poly_t poly, unsigned long k)

    #
    # Polynomial multiplication
    #
    # All multiplication functions require that the modulus be no more than FLINT_BITS-1 bits
    #

    cdef void nmod_poly_mul(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_sqr(nmod_poly_t res, nmod_poly_t poly)

    # Requires that poly1 bits + poly2 bits + log_length is not greater than 2*FLINT_BITS

    cdef void nmod_poly_mul_KS(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits_input)
    cdef void _nmod_poly_mul_KS(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits_input)

    cdef void nmod_poly_mul_KS_trunc(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits_input, unsigned long trunc)
    cdef void _nmod_poly_mul_KS_trunc(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits_input, unsigned long trunc)

    cdef void _nmod_poly_mul_classical(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void __nmod_poly_mul_classical_mod_last(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits)
    cdef void __nmod_poly_mul_classical_mod_throughout(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits)
    cdef void nmod_poly_mul_classical(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void _nmod_poly_sqr_classical(nmod_poly_t res, nmod_poly_t poly)
    cdef void nmod_poly_sqr_classical(nmod_poly_t res, nmod_poly_t poly)

    cdef void _nmod_poly_mul_classical_trunc(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long trunc)
    cdef void __nmod_poly_mul_classical_trunc_mod_last(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void __nmod_poly_mul_classical_trunc_mod_throughout(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void nmod_poly_mul_classical_trunc(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long trunc)

    cdef void _nmod_poly_mul_classical_trunc_left(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long trunc)
    cdef void __nmod_poly_mul_classical_trunc_left_mod_last(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void __nmod_poly_mul_classical_trunc_left_mod_throughout(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void nmod_poly_mul_classical_trunc_left(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long trunc)

    cdef void nmod_poly_mullow(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long trunc)
    cdef void nmod_poly_mulhigh(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, unsigned long trunc)

    #
    # Bit packing functions
    #

    cdef unsigned long nmod_poly_bits(nmod_poly_t poly)
    cdef void _nmod_poly_bit_pack_mpn(mp_limb_t * res, nmod_poly_t poly, unsigned long bits, unsigned long length)
    cdef void _nmod_poly_bit_unpack_mpn(nmod_poly_t poly, mp_limb_t *mpn, unsigned long length, unsigned long bits)

    cdef void print_binary(unsigned long n, unsigned long len)
    cdef void print_binary2(unsigned long n, unsigned long len, unsigned long space_bit)

    #
    # Scalar multiplication
    #

    cdef void _nmod_poly_scalar_mul_nmod(nmod_poly_t res, nmod_poly_t poly, unsigned long scalar)
    cdef void nmod_poly_scalar_mul_nmod(nmod_poly_t res, nmod_poly_t poly, unsigned long scalar)
    cdef void __nmod_poly_scalar_mul_without_mod(nmod_poly_t res, nmod_poly_t poly, unsigned long scalar)

    #
    # Division
    #

    cdef void nmod_poly_divrem_classical(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)
    cdef void __nmod_poly_divrem_classical_mod_last(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_div_classical(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B)
    cdef void __nmod_poly_div_classical_mod_last(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_div_divconquer_recursive(nmod_poly_t Q, nmod_poly_t BQ, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_divrem_divconquer(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_div_divconquer(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B)

    #
    # Newton Inversion
    #

    cdef void nmod_poly_newton_invert_basecase(nmod_poly_t Q_inv, nmod_poly_t Q, unsigned long n)
    cdef void nmod_poly_newton_invert(nmod_poly_t Q_inv, nmod_poly_t Q, unsigned long n)

    #
    # Newton Division
    #

    cdef void nmod_poly_div_series(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B, unsigned long n)
    cdef void nmod_poly_div_newton(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_divrem_newton(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)

    cdef void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)

    cdef void nmod_poly_div(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B)

    #
    # Resultant
    #

    cdef unsigned long nmod_poly_resultant_euclidean(nmod_poly_t a, nmod_poly_t b)

    cdef unsigned long nmod_poly_resultant(nmod_poly_t a, nmod_poly_t b)

    #
    # GCD
    #

    cdef void nmod_poly_gcd(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef int nmod_poly_gcd_invert(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_xgcd(nmod_poly_t res, nmod_poly_t s, nmod_poly_t t, nmod_poly_t poly1, nmod_poly_t poly2)



    # Composition / evaluation

    cdef unsigned long nmod_poly_evaluate_nmod(nmod_poly_t, unsigned long)
    cdef void nmod_poly_compose(nmod_poly_t, nmod_poly_t, nmod_poly_t)

    # Factorization

    cdef bint nmod_poly_is_irreducible(nmod_poly_t p)

    ctypedef struct nmod_poly_factors_struct:
        unsigned long num_factors
        unsigned long* exponents
        nmod_poly_t* factors

    ctypedef nmod_poly_factors_struct* nmod_poly_factor_t

    cdef void nmod_poly_factor_init(nmod_poly_factor_t)
    cdef void nmod_poly_factor_clear(nmod_poly_factor_t)
    cdef unsigned long nmod_poly_factor(nmod_poly_factor_t, nmod_poly_t)
    cdef void nmod_poly_factor_squarefree(nmod_poly_factor_t, nmod_poly_t)
    cdef void nmod_poly_factor_berlekamp(nmod_poly_factor_t factors, nmod_poly_t f)

    cdef void nmod_poly_factor_add(nmod_poly_factor_t fac, nmod_poly_t poly)
    cdef void nmod_poly_factor_concat(nmod_poly_factor_t res, nmod_poly_factor_t fac)
    cdef void nmod_poly_factor_print(nmod_poly_factor_t fac)
    cdef void nmod_poly_factor_pow(nmod_poly_factor_t fac, unsigned long exp)

    #
    # Differentiation
    #

    cdef void nmod_poly_derivative(nmod_poly_t res, nmod_poly_t poly)

    #
    # Arithmetic modulo a polynomial
    #

    cdef void nmod_poly_mulmod(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, nmod_poly_t f)
    cdef void nmod_poly_powmod(nmod_poly_t res,nmod_poly_t pol, long exp, nmod_poly_t f)
