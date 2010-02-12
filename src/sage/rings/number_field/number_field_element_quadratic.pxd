include "../../ext/cdefs.pxi"
include "../../libs/ntl/decl.pxi"

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.structure.element cimport Element, FieldElement, RingElement, ModuleElement
from sage.structure.parent_base cimport ParentWithBase
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ

from number_field_element cimport NumberFieldElement, NumberFieldElement_absolute

cdef class NumberFieldElement_quadratic(NumberFieldElement_absolute):
    # (a + b sqrt(D)) / denom
    cdef mpz_t a, b, denom
    cdef Integer D
    cdef NumberFieldElement conjugate_c(self)
    cdef bint is_sqrt_disc(self)

    cdef void _randomize(self, num_bound, den_bound, distribution)

cdef class OrderElement_quadratic(NumberFieldElement_quadratic):
    pass
