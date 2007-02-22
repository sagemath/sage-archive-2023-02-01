# choose: dense or sparse
cimport matrix_dense
#from sage.structure.element cimport ModuleElement
include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../gsl/gsl.pxi'
include '../ext/python.pxi'
import matrix
cimport matrix
cdef class Matrix_real_double_dense(matrix_dense.Matrix_dense):
    cdef gsl_matrix *_matrix
    cdef gsl_matrix *_LU
    cdef gsl_permutation *_p
    cdef int _signum
    cdef int _LU_valid
    cdef _c_compute_LU(self)
#    cdef ModuleElement _add_c_impl(self, ModuleElement right )
