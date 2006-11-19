cimport free_module_element
import  free_module_element
include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
include '../gsl/gsl.pxi'


cdef class RealDoubleVectorSpace_element(free_module_element.FreeModuleElement):
	cdef gsl_vector * v
	cdef _new_c(self, gsl_vector* v)
	cdef gsl_vector* gsl_vector_copy(self)

