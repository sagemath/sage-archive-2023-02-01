from sage.modules.complex_double_vector cimport ComplexDoubleVectorSpaceElement
from matrix_complex_double_dense cimport Matrix_complex_double_dense
from sage.modules.free_module_element import vector
import sage.rings.complex_double
from matrix_real_double_dense cimport Matrix_real_double_dense
from sage.modules.real_double_vector cimport RealDoubleVectorSpaceElement
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
#       if _mat._LU_valid == 0:
       if _mat.fetch('LU_valid')!=True:
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


def solve_matrix_complex_double_dense(mat, vec):
       cdef ComplexDoubleVectorSpaceElement _vec
       cdef Matrix_complex_double_dense _mat
       _vec = vec
       _mat = mat
       cdef gsl_vector_complex * result_vector
       cdef gsl_vector_complex *b
       cdef gsl_permutation *p
       cdef gsl_matrix_complex *LU
       cdef gsl_complex z
       if _mat.fetch('LU_valid')!=True:
           _mat._c_compute_LU()
       result_vector = <gsl_vector_complex *> gsl_vector_complex_alloc(_mat._matrix.size1)
       cdef int result_solve
       if not isinstance(vec, ComplexDoubleVectorSpaceElement):
           raise TypeError, "vector must be in CDF"
       p = (<gsl_permutation*> _mat._p)
       LU =  <gsl_matrix_complex *> _mat._LU
       b = <gsl_vector_complex*> ( <ComplexDoubleVectorSpaceElement> _vec).v

       result_solve = gsl_linalg_complex_LU_solve(LU,p,b,result_vector) #check error_code
       list = []
       for i from 0<=i< _mat._matrix.size1:
              z=gsl_vector_complex_get(result_vector,i)
              list.append(sage.rings.complex_double.CDF(GSL_REAL(z),GSL_IMAG(z)))
       gsl_vector_complex_free(result_vector)
       return vector(sage.rings.complex_double.CDF,list) #todo: don't go through python

