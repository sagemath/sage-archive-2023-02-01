include "sage/ext/cdefs.pxi"
include "sage/libs/flint/fmpz_poly.pxi"

from sage.libs.pari.gen cimport gen as pari_gen
ctypedef fmpz_poly_t celement

include "CR_template_header.pxi"

cdef class qAdicCappedRelativeElement(CRElement):
    pass
