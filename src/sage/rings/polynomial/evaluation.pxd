from sage.libs.flint.fmpz_poly cimport fmpz_poly_t
from sage.libs.ntl.types cimport ZZX_c
from sage.libs.mpfr cimport mpfr_t, mpfr_rnd_t
from sage.libs.mpfi cimport mpfi_t

cdef fmpz_poly_eval_mpfr(mpfr_t res, const fmpz_poly_t poly, const mpfr_t a, mpfr_rnd_t rnd)
cdef fmpz_poly_eval_mpfi(mpfi_t res, const fmpz_poly_t poly, const mpfi_t a)

cdef ZZX_eval_mpfr(mpfr_t res, ZZX_c poly, const mpfr_t a, const mpfr_rnd_t rnd)
cdef ZZX_eval_mpfi(mpfi_t res, ZZX_c poly, const mpfi_t a)

