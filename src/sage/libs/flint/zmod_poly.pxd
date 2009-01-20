include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"

from sage.libs.flint.flint cimport *

from flint import *

cdef extern from "FLINT/zmod_poly.h":
    ctypedef struct zmod_poly_struct:
        unsigned long *coeffs
        unsigned long alloc
        unsigned long length
        unsigned long p
        double p_inv

    ctypedef zmod_poly_struct* zmod_poly_t

    cdef void zmod_poly_init(zmod_poly_t poly, unsigned long p)
    cdef void zmod_poly_init_precomp(zmod_poly_t poly, unsigned long p, double p_inv)
    cdef void zmod_poly_init2(zmod_poly_t poly, unsigned long p, unsigned long alloc)
    cdef void zmod_poly_init2_precomp(zmod_poly_t poly, unsigned long p, double p_inv, unsigned long alloc)
    cdef void zmod_poly_clear(zmod_poly_t poly)

    cdef void zmod_poly_realloc(zmod_poly_t poly, unsigned long alloc)
    # _bits_ only applies to newly allocated coefficients, not existing ones...

    # this non-inlined version REQUIRES that alloc > poly->alloc
    void __zmod_poly_fit_length(zmod_poly_t poly, unsigned long alloc)

    # this is arranged so that the initial comparison (very frequent) is inlined,
    # but the actual allocation (infrequent) is not
    cdef void zmod_poly_fit_length(zmod_poly_t poly, unsigned long alloc)

    # ------------------------------------------------------
    # Setting/retrieving coefficients

    cdef unsigned long zmod_poly_get_coeff_ui(zmod_poly_t poly, unsigned long n)

    cdef unsigned long _zmod_poly_get_coeff_ui(zmod_poly_t poly, unsigned long n)

    cdef void zmod_poly_set_coeff_ui(zmod_poly_t poly, unsigned long n, unsigned long c)

    cdef void _zmod_poly_set_coeff_ui(zmod_poly_t poly, unsigned long n, unsigned long c)

    # ------------------------------------------------------
    # String conversions and I/O

    cdef int zmod_poly_from_string(zmod_poly_t poly, char* s)
    cdef char* zmod_poly_to_string(zmod_poly_t poly)
    cdef void zmod_poly_print(zmod_poly_t poly)
    cdef void zmod_poly_fprint(zmod_poly_t poly, FILE* f)
    cdef int zmod_poly_read(zmod_poly_t poly)
    cdef int zmod_poly_fread(zmod_poly_t poly, FILE* f)

    # ------------------------------------------------------
    # Length and degree

    cdef void __zmod_poly_normalise(zmod_poly_t poly)
    cdef int __zmod_poly_normalised(zmod_poly_t poly)
    cdef void zmod_poly_truncate(zmod_poly_t poly, unsigned long length)

    cdef unsigned long zmod_poly_length(zmod_poly_t poly)

    cdef long zmod_poly_degree(zmod_poly_t poly)

    cdef unsigned long zmod_poly_modulus(zmod_poly_t poly)

    cdef double zmod_poly_precomputed_inverse(zmod_poly_t poly)

    # ------------------------------------------------------
    # Assignment

    cdef void _zmod_poly_set(zmod_poly_t res, zmod_poly_t poly)
    cdef void zmod_poly_set(zmod_poly_t res, zmod_poly_t poly)

    cdef void zmod_poly_zero(zmod_poly_t poly)

    cdef void zmod_poly_swap(zmod_poly_t poly1, zmod_poly_t poly2)

    #
    # Subpolynomials
    #

    cdef void _zmod_poly_attach(zmod_poly_t output, zmod_poly_t input)

    cdef void zmod_poly_attach(zmod_poly_t output, zmod_poly_t input)

    #
    # Attach input shifted right by n to output
    #

    cdef void _zmod_poly_attach_shift(zmod_poly_t output, zmod_poly_t input, unsigned long n)

    cdef void zmod_poly_attach_shift(zmod_poly_t output, zmod_poly_t input, unsigned long n)

    #
    # Attach input to first n coefficients of input
    #

    cdef void _zmod_poly_attach_truncate(zmod_poly_t output,  zmod_poly_t input, unsigned long n)

    cdef void zmod_poly_attach_truncate(zmod_poly_t output, zmod_poly_t input, unsigned long n)

    #
    # Comparison functions
    #

    cdef int zmod_poly_equal(zmod_poly_t poly1, zmod_poly_t poly2)

    cdef int zmod_poly_is_one(zmod_poly_t poly1)

    #
    # Reversal
    #

    cdef void _zmod_poly_reverse(zmod_poly_t output, zmod_poly_t input, unsigned long length)
    cdef void zmod_poly_reverse(zmod_poly_t output, zmod_poly_t input, unsigned long length)

    #
    # Monic polys
    #

    cdef void zmod_poly_make_monic(zmod_poly_t output, zmod_poly_t pol)

    #
    # Addition and subtraction
    #

    cdef void zmod_poly_add(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void zmod_poly_add_without_mod(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void zmod_poly_sub(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void _zmod_poly_sub(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void zmod_poly_neg(zmod_poly_t res, zmod_poly_t poly)

    #
    # Shifting functions
    #

    cdef void zmod_poly_left_shift(zmod_poly_t res, zmod_poly_t poly, unsigned long k)
    cdef void zmod_poly_right_shift(zmod_poly_t res, zmod_poly_t poly, unsigned long k)

    #
    # Polynomial multiplication
    #
    # All multiplication functions require that the modulus be no more than FLINT_BITS-1 bits
    #

    cdef void zmod_poly_mul(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void zmod_poly_sqr(zmod_poly_t res, zmod_poly_t poly)

    # Requires that poly1 bits + poly2 bits + log_length is not greater than 2*FLINT_BITS

    cdef void zmod_poly_mul_KS(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits_input)
    cdef void _zmod_poly_mul_KS(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits_input)

    cdef void zmod_poly_mul_KS_trunc(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits_input, unsigned long trunc)
    cdef void _zmod_poly_mul_KS_trunc(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits_input, unsigned long trunc)

    cdef void _zmod_poly_mul_classical(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void __zmod_poly_mul_classical_mod_last(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits)
    cdef void __zmod_poly_mul_classical_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits)
    cdef void zmod_poly_mul_classical(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void _zmod_poly_sqr_classical(zmod_poly_t res, zmod_poly_t poly)
    cdef void zmod_poly_sqr_classical(zmod_poly_t res, zmod_poly_t poly)

    cdef void _zmod_poly_mul_classical_trunc(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
    cdef void __zmod_poly_mul_classical_trunc_mod_last(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void __zmod_poly_mul_classical_trunc_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void zmod_poly_mul_classical_trunc(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)

    cdef void _zmod_poly_mul_classical_trunc_left(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
    cdef void __zmod_poly_mul_classical_trunc_left_mod_last(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void __zmod_poly_mul_classical_trunc_left_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
    cdef void zmod_poly_mul_classical_trunc_left(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)

    cdef void zmod_poly_mul_trunc_n(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
    cdef void zmod_poly_mul_trunc_left_n(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)

    #
    # Bit packing functions
    #

    cdef unsigned long zmod_poly_bits(zmod_poly_t poly)
    cdef void _zmod_poly_bit_pack_mpn(mp_limb_t * res, zmod_poly_t poly, unsigned long bits, unsigned long length)
    cdef void _zmod_poly_bit_unpack_mpn(zmod_poly_t poly, mp_limb_t *mpn, unsigned long length, unsigned long bits)

    cdef void print_binary(unsigned long n, unsigned long len)
    cdef void print_binary2(unsigned long n, unsigned long len, unsigned long space_bit)

    #
    # Scalar multiplication
    #

    cdef void _zmod_poly_scalar_mul(zmod_poly_t res, zmod_poly_t poly, unsigned long scalar)
    cdef void zmod_poly_scalar_mul(zmod_poly_t res, zmod_poly_t poly, unsigned long scalar)
    cdef void __zmod_poly_scalar_mul_without_mod(zmod_poly_t res, zmod_poly_t poly, unsigned long scalar)

    #
    # Division
    #

    cdef void zmod_poly_divrem_classical(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)
    cdef void __zmod_poly_divrem_classical_mod_last(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)
    cdef void zmod_poly_div_classical(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)
    cdef void __zmod_poly_div_classical_mod_last(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)
    cdef void zmod_poly_div_divconquer_recursive(zmod_poly_t Q, zmod_poly_t BQ, zmod_poly_t A, zmod_poly_t B)
    cdef void zmod_poly_divrem_divconquer(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)
    cdef void zmod_poly_div_divconquer(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)

    #
    # Newton Inversion
    #

    cdef void zmod_poly_newton_invert_basecase(zmod_poly_t Q_inv, zmod_poly_t Q, unsigned long n)
    cdef void zmod_poly_newton_invert(zmod_poly_t Q_inv, zmod_poly_t Q, unsigned long n)

    #
    # Newton Division
    #

    cdef void zmod_poly_div_series(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B, unsigned long n)
    cdef void zmod_poly_div_newton(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)
    cdef void zmod_poly_divrem_newton(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)

    cdef void zmod_poly_divrem(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)

    cdef void zmod_poly_div(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)

    #
    # Resultant
    #

    cdef unsigned long zmod_poly_resultant_euclidean(zmod_poly_t a, zmod_poly_t b)

    cdef unsigned long zmod_poly_resultant(zmod_poly_t a, zmod_poly_t b)

    #
    # GCD
    #

    cdef void zmod_poly_gcd(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef int zmod_poly_gcd_invert(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
    cdef void zmod_poly_xgcd(zmod_poly_t res, zmod_poly_t s, zmod_poly_t t, zmod_poly_t poly1, zmod_poly_t poly2)


    # FLINT 1.1 will have:
    # cdef void zmod_poly_powmod(zmod_poly_t res, zmod_poly_t pol, long exp, zmod_poly_t f)
