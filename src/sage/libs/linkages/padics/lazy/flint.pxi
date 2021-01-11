from sage.libs.flint.types cimport flint_rand_t
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
cdef extern from "sage/libs/flint/flint_wrap.h":
    cdef ulong fmpz_bits(fmpz_t f)
    cdef int fmpz_tstbit(fmpz_t f, ulong i)

cdef extern from "sage/libs/linkages/padics/lazy/flint_helper.c":
    cdef void flint_randinit(flint_rand_t state)
    cdef void flint_randclear(flint_rand_t state)
    cdef fmpz* get_coeff(fmpz_poly_t poly, slong i)
    cdef void get_slice(fmpz_poly_t slice, fmpz_poly_t poly, slong start, slong length)
    cdef void iadd_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void isub_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void iadd_shifted(fmpz_poly_t poly, fmpz_poly_t summand, slong shift)
    cdef void reduce_coeff(fmpz_poly_t poly, slong i, fmpz_t modulus)

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.pow_computer_flint cimport PowComputer_flint

from sage.ext.stdsage cimport PY_NEW


cdef flint_rand_t flint_randstate 
flint_randinit(flint_randstate)

# Operations on digits (intended to be small elements in the exact subring)

cdef inline void digit_init(cdigit a):
    fmpz_init(a)

cdef inline void digit_clear(cdigit a):
    fmpz_clear(a)

cdef inline void digit_set(cdigit a, cdigit b):
    fmpz_set(a, b)

cdef inline void digit_set_ui(cdigit a, slong b):
    fmpz_set_ui(a, b)

cdef inline bint digit_equal(cdigit a, cdigit b):
    return fmpz_equal(a, b)

cdef inline bint digit_is_zero(cdigit a):
    return fmpz_is_zero(a)

cdef inline bint digit_equal_ui(cdigit a, slong b):
    return fmpz_equal_ui(a, b)

cdef inline void digit_set_sage(cdigit a, Integer elt):
    fmpz_set_mpz(a, elt.value)

cdef inline Integer digit_get_sage(cdigit a):
    cdef Integer elt = PY_NEW(Integer)
    fmpz_get_mpz(elt.value, a)
    return elt

cdef inline void digit_random(cdigit res, PowComputer_class prime_pow):
    fmpz_randm(res, flint_randstate, (<PowComputer_flint>prime_pow).fprime)

cdef inline void digit_add(cdigit res, cdigit a, cdigit b):
    fmpz_add(res, a, b)

cdef inline void digit_sub(cdigit res, cdigit a, cdigit b):
    fmpz_sub(res, a, b)

cdef inline void digit_mul(cdigit res, cdigit a, cdigit b):
    fmpz_mul(res, a, b)

cdef inline void digit_mod(cdigit res, cdigit a, PowComputer_class prime_pow):
    fmpz_mod(res, a, (<PowComputer_flint>prime_pow).fprime)

cdef inline void digit_quorem(cdigit quo, cdigit rem, cdigit a, PowComputer_class prime_pow):
    fmpz_tdiv_qr(quo, rem, a, (<PowComputer_flint>prime_pow).fprime)

cdef inline void digit_inv(cdigit res, cdigit a, PowComputer_class prime_pow):
    cdef cdigit gcd
    fmpz_init(gcd)
    fmpz_gcdinv(gcd, res, a, (<PowComputer_flint>prime_pow).fprime)
    fmpz_clear(gcd)

cdef bint digit_sqrt(cdigit_ptr ans, cdigit_ptr x, PowComputer_class prime_pow):
    # Need to do something better
    cdef Integer zx = digit_get_sage(x)
    cdef Integer zp = digit_get_sage((<PowComputer_flint>prime_pow).fprime)
    try:
        k = GF(zp)
        zans = ZZ(k(zx).sqrt(extend=False))
    except ValueError:
        return 1
    digit_set_sage(ans, zans)
    return 0


# Operations on elements (represented as series of digits)

cdef inline void element_init(celement x):
    fmpz_poly_init(x)

cdef inline void element_clear(celement x):
    fmpz_poly_clear(x)

cdef inline void element_set(celement x, celement y):
    fmpz_poly_set(x, y)

cdef inline void element_set_coeff(celement x, cdigit a, slong i):
    fmpz_poly_set_coeff_fmpz(x, i, a)

cdef inline void element_set_coeff_ui(celement x, slong a, slong i):
    fmpz_poly_set_coeff_ui(x, i, a)

cdef inline void element_set_coeff_sage(celement x, Integer a, slong i):
    fmpz_poly_set_coeff_mpz(x, i, a.value)

cdef inline cdigit_ptr element_get_coeff(celement x, slong i):
    return get_coeff(x, i)

cdef inline Integer element_get_sage(celement x, PowComputer_class prime_pow):
    cdef fmpz_t value
    fmpz_init(value)
    fmpz_poly_evaluate_fmpz(value, x, (<PowComputer_flint>prime_pow).fprime)
    cdef Integer ans = digit_get_sage(value)
    fmpz_clear(value)
    return ans

cdef inline Integer element_get_coeff_sage(celement x, slong i):
    return digit_get_sage(get_coeff(x, i))

cdef inline void element_get_slice(celement res, celement x, slong start, slong length):
    get_slice(res, x, start, length)

cdef inline void element_iadd_coeff(celement x, cdigit a, slong i):
    iadd_coeff(x, a, i)

cdef inline void element_isub_coeff(celement x, cdigit a, slong i):
    isub_coeff(x, a, i)

cdef inline void element_iadd_slice(celement x, celement slice, slong start):
    iadd_shifted(x, slice, start)

cdef inline void element_isub_slice(celement x, celement slice, slong start):
    raise NotImplementedError

cdef inline void element_scalarmul(celement res, celement x, cdigit a):
    fmpz_poly_scalar_mul_fmpz(res, x, a)

cdef inline void element_mul(celement res, celement x, celement y):
    fmpz_poly_mul(res, x, y)

cdef inline void element_reduce_coeff(celement x, slong i, PowComputer_class prime_pow):
    reduce_coeff(x, i, (<PowComputer_flint>prime_pow).fprime)

cdef inline void element_shift_right(celement x):
    fmpz_poly_shift_right(x, x, 1)
