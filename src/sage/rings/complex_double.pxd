include '../gsl/gsl_complex.pxi'
include '../libs/pari/decl.pxi'

cimport sage.structure.element
cimport sage.rings.ring
cimport sage.libs.pari.gen

import sage.structure.coerce

import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport RingElement


cdef class ComplexDoubleField_class(sage.rings.ring.Field):
	pass

cdef class ComplexDoubleElement(sage.structure.element.FieldElement):
	cdef gsl_complex _complex

	cdef int cmp(ComplexDoubleElement left,ComplexDoubleElement right)

	cdef GEN _gen(self)
	cdef RingElement _add_sibling_cdef(self, RingElement right)

