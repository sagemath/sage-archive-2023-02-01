#from sage.libs.flint.fmpz_poly cimport fmpz_poly_t
include "sage/libs/flint/fmpz_poly.pxi"

include "sage/ext/cdefs.pxi"

ctypedef fmpz_poly_t celement

include "CR_template_header.pxi"

cdef class qAdicCappedRelativeElement(CRElement):
    pass
