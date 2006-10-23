include '../gsl/gsl_complex.pxi'
include '../libs/pari/decl.pxi'

cimport sage.structure.element
cimport sage.rings.ring
cimport sage.libs.pari.gen

import sage.structure.coerce
#import complex_number

import integer_ring
import infinity

cdef class ComplexDoubleField_class(sage.rings.ring.Field):
	pass

cdef class ComplexDoubleElement(sage.structure.element.FieldElement):
	cdef gsl_complex _complex


	cdef int cmp(ComplexDoubleElement left,ComplexDoubleElement right)

	cdef GEN _gen(self)


#cdef new_element(gsl_complex x)

#cdef  new_from_gen(GEN g, pari_sp sp)

#cdef GEN complex_gen(x)

#cdef _add_(ComplexDoubleElement self,ComplexDoubleElement right)

