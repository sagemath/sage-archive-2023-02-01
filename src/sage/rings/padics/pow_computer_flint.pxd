include "sage/ext/cdefs.pxi"
include "sage/libs/flint/fmpz.pxi"
include "sage/libs/flint/fmpz_poly.pxi"
include "sage/libs/flint/padic.pxi"

#from sage.libs.flint.ntl_interface cimport *
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class PowComputer_flint_base(PowComputer_class):
    cdef padic_ctx_t ctx
    cdef fmpz_t fprime
    cdef fmpz_t ftmp
    cdef fmpz_t ftmp2
    cdef mpz_t top_power

    cdef fmpz_t* pow_fmpz_t_tmp(self, unsigned long n)
