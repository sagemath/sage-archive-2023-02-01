include "sage/ext/cdefs.pxi"
include "sage/ext/interrupt.pxi"

from sage.libs.flint.padic cimport padic_ctx_t
from sage.libs.flint.fmpz cimport fmpz_t
from sage.libs.flint.fmpz_poly cimport fmpz_poly_t
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class PowComputer_flint(PowComputer_class):
    cdef padic_ctx_t ctx
    cdef fmpz_t fprime
    cdef fmpz_t ftmp
    cdef fmpz_t ftmp2
    cdef mpz_t top_power

    cdef fmpz_t* pow_fmpz_t_tmp(self, unsigned long n) except NULL
    cdef unsigned long capdiv(self, unsigned long n)

cdef class PowComputer_flint_1step(PowComputer_flint):
    cdef bint is_monic
    cdef fmpz_t q
    cdef fmpz_t _an
    cdef fmpz_t* _inv_an
    cdef fmpz_poly_t modulus
    cdef fmpz_poly_t tmp_poly
    cdef fmpz_poly_t tmp_poly2
    cdef fmpz_poly_t* _moduli
    cdef fmpz_poly_t* get_modulus(self, unsigned long n)
    cdef fmpz_poly_t* get_modulus_capdiv(self, unsigned long n)
    cdef fmpz_t* get_inv_an(self, unsigned long k)
    cdef fmpz_t* get_inv_an_capdiv(self, unsigned long k)

cdef class PowComputer_flint_unram(PowComputer_flint_1step):
    pass

cdef class PowComputer_flint_eis(PowComputer_flint_1step):
    pass
