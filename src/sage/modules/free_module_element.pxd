#include  '../ext/cdefs.pxi'
#include  '../ext/interrupt.pxi'

cimport sage.structure.element
import  sage.structure.element
cimport sage.matrix.matrix_generic
cdef class FreeModuleElement(sage.structure.element.ModuleElement):
	cdef object __entries
	cdef object __zero
	cdef FreeModuleElement _scalar_multiply(self,scalar)


	cdef FreeModuleElement _scalar_multiply_coerce(self,scalar)

	cdef FreeModuleElement _matrix_multiply(self,sage.matrix.matrix_generic.Matrix A)

