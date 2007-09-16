include "../../ext/cdefs.pxi"
include "../../libs/ntl/decl.pxi"

import sage.structure.element
cimport sage.structure.element
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.element cimport FieldElement, RingElement, ModuleElement
from sage.structure.parent_base cimport ParentWithBase

cdef class NumberFieldElement(FieldElement):
    cdef ZZX_c __numerator
    cdef ZZ_c __denominator
    cdef object __multiplicative_order
    cdef object __pari
    cdef object __matrix

    cdef NumberFieldElement _new(self)

    cdef void _parent_poly_c_(self, ZZX_c *num, ZZ_c *den)
    cdef void _invert_c_(self, ZZX_c *num, ZZ_c *den)
    cdef void _reduce_c_(self)
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef ModuleElement _neg_c_impl(self)

cdef class NumberFieldElement_absolute(NumberFieldElement):
    pass

cdef class NumberFieldElement_quadratic(NumberFieldElement_absolute):
    # (a + b sqrt(disc)) / denom
    cdef mpz_t a, b, denom
    cdef Integer disc
    cdef NumberFieldElement conjugate_c(self)

cdef class NumberFieldElement_relative(NumberFieldElement):
    pass

# TODO: cyclotomic and/or quadratic classes? (Both for differing implementations and speed).

cdef class OrderElement_absolute(NumberFieldElement_absolute):
    cdef object _order

cdef class OrderElement_relative(NumberFieldElement_relative):
    cdef object _order
