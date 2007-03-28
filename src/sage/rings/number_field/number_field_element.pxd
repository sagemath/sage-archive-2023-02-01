include "../../ext/cdefs.pxi"
include "../../libs/ntl/decl.pxi"

import sage.structure.element
cimport sage.structure.element
from sage.rings.integer cimport Integer
from sage.rings.polynomial_element cimport Polynomial
from sage.structure.element cimport FieldElement, RingElement, ModuleElement
from sage.structure.parent_base cimport ParentWithBase

cdef class NumberFieldElement(FieldElement):
    cdef ntl_c_ZZX __numerator
    cdef ntl_c_ZZ __denominator
    cdef object __multiplicative_order
    cdef object __pari
    cdef object __matrix

    cdef NumberFieldElement _new(self)

    cdef void _parent_poly_c_(self, ntl_c_ZZX *num, ntl_c_ZZ *den)
    cdef void _invert_c_(self, ntl_c_ZZX *num, ntl_c_ZZ *den)
    cdef void _reduce_c_(self)
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef ModuleElement _neg_c_impl(self)
