cimport free_module_element
import  free_module_element
include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
include '../gsl/gsl.pxi'


cdef class ComplexDoubleVectorSpace_element(free_module_element.FreeModuleElement):
	cdef gsl_vector_complex * v
	cdef _add_(ComplexDoubleVectorSpace_element left, ComplexDoubleVectorSpace_element right)
