include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"

from flint cimport *
include "fmpz_poly.pxi"

from sage.structure.sage_object cimport SageObject

cdef class Fmpz_poly(SageObject):
    cdef fmpz_poly_t poly
