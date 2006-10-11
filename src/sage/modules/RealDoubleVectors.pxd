cimport sage.modules.free_module_element
import sage.modules.free_module_element
#ctypedef int size_t
include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
#include '../gsl/gsl.pxi'


cdef class RealDoubleVectorSpace_element(sage.modules.free_module_element.FreeModuleElement):
	cdef double* vec
	cdef size_t n
	cdef size_t stride
#	ddef void __dealloc__(self)