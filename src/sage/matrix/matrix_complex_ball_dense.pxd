from sage.libs.arb.acb_mat cimport acb_mat_t
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.structure.parent cimport Parent
from sage.structure.sage_object cimport SageObject

cdef void matrix_to_acb_mat(acb_mat_t target, source)
cdef Matrix_generic_dense acb_mat_to_matrix(
    acb_mat_t source, Parent CIF)

cdef class Acb_mat(SageObject):
    cdef acb_mat_t value
    cdef unsigned long _precision_
    cpdef  Matrix_generic_dense _matrix_(self)
