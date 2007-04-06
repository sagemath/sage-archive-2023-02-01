ctypedef unsigned long ulong

include "../../ext/cdefs.pxi"

cimport sage.matrix.matrix_dense
cimport sage.matrix.matrix_integer_dense
cimport sage.rings.integer
cimport sage.structure.element

from sage.rings.integer cimport Integer
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.structure.element cimport RingElement

cdef class Matrix_generic_dense(Matrix_dense):
    cdef char _initialized
    cdef Matrix_integer_dense _value_matrix
    cdef ulong *_relprecs
    cdef ulong **_relprec_mat
    cdef object _valaddeds
    cdef object _padic_values

    cdef void _comp_valaddeds(self)
    cdef void _adjust_prec_info_global(self, RingElement absolute, RingElement relative)
    cdef void _adjust_prec_info_global_local(self, RingElement absolute, object relative)
    cdef void _adjust_prec_info_local_global(self, object absolute, RingElement relative)
    cdef void _adjust_prec_info_local_local(self, object absolute, object relative)
    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object value)
    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j)

