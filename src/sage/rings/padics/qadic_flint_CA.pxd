from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.flint.fmpz_poly cimport fmpz_poly_t
from sage.rings.padics.pow_computer_flint cimport PowComputer_flint_unram
from sage.rings.padics.qadic_flint_CR cimport CRElement

cdef class PowComputer_(PowComputer_flint_unram):
    pass
ctypedef fmpz_poly_t celement

include "CA_template_header.pxi"

cdef class qAdicCappedAbsoluteElement(CAElement):
    pass

cdef class qAdicCoercion_Zq_Qq(RingHomomorphism_coercion):
    cdef CRElement _zero
    cdef Morphism _section

cdef class qAdicConvert_Qq_Zq(Morphism):
    cdef CAElement _zero
