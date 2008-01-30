include '../gsl/gsl_complex.pxi'
include '../libs/pari/decl.pxi'

cimport sage.structure.element
cimport sage.rings.ring
cimport sage.libs.pari.gen

import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport RingElement, ModuleElement


cdef class ComplexDoubleField_class(sage.rings.ring.Field):
	pass

cdef class ComplexDoubleElement(sage.structure.element.FieldElement):
	cdef gsl_complex _complex
	cdef GEN _gen(self)
	cdef ComplexDoubleElement _new_c(self, gsl_complex x)
	cdef _new_from_gen_c(self, GEN g, pari_sp sp)

cdef public api ComplexDoubleElement new_ComplexDoubleElement()