from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.flint.fmpz_poly cimport fmpz_poly_t
from sage.rings.padics.pow_computer_flint cimport PowComputer_flint_unram
cdef class PowComputer_(PowComputer_flint_unram):
    pass
ctypedef fmpz_poly_t celement

include "FM_template_header.pxi"

cdef class qAdicFixedModElement(FMElement):
    pass
