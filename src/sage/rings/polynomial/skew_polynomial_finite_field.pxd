from sage.rings.integer cimport Integer
from skew_polynomial_element cimport SkewPolynomial_generic_dense
from sage.matrix.matrix_dense cimport Matrix_dense
from polynomial_element cimport Polynomial
from sage.structure.element cimport RingElement

cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_generic_dense):
    
    cdef SkewPolynomial_finite_field_dense _rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lmonic(self)
    cdef void _inplace_rmonic(self)
    cdef void _inplace_rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef Py_ssize_t _val_inplace_unit(self)
    cdef SkewPolynomial_finite_field_dense _rquo_inplace_rem(self, SkewPolynomial_finite_field_dense other)

    cdef Matrix_dense _matmul_c(self)

    cpdef _leftpow_(self, exp, modulus=*)
    cpdef _rightpow_(self, exp, modulus=*)
