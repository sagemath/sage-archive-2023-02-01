from sage.modules.real_double_vector cimport RealDoubleVectorSpaceElement
from matrix_real_double_dense cimport Matrix_real_double_dense
from sage.modules.free_module_element import vector
import sage.rings.real_double
include '../ext/interrupt.pxi'
include '../ext/stdsage.pxi'
include '../ext/cdefs.pxi'
include '../ext/python.pxi'
include '../gsl/gsl.pxi'


def solve_matrix_real_double_dense(mat, vec):
       cdef RealDoubleVectorSpaceElement _vec
       cdef Matrix_real_double_dense _mat
       _vec = vec
       _mat = mat
       cdef gsl_vector * result_vector
       cdef gsl_vector *b
       cdef gsl_permutation *p
       cdef gsl_matrix *LU
       if _mat._LU_valid == 0:
           _mat._c_compute_LU()
       result_vector = <gsl_vector *> gsl_vector_alloc(_mat._matrix.size1)
       cdef int result_solve
       if not isinstance(vec, RealDoubleVectorSpaceElement):
           raise TypeError, "vector must be in RDF"
       p = (<gsl_permutation*> _mat._p)
       LU =  <gsl_matrix *> _mat._LU
       b = <gsl_vector*> ( <RealDoubleVectorSpaceElement> _vec).v

       result_solve = gsl_linalg_LU_solve(LU,p,b,result_vector) #check error_code
       list = []
       for i from 0<=i< _mat._matrix.size1:
           list.append(gsl_vector_get(result_vector, i))
       gsl_vector_free(result_vector)
       return vector(sage.rings.real_double.RDF,list) #todo: don't go through python

