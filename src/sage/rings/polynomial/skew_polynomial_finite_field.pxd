from sage.rings.integer cimport Integer

from skew_polynomial_element cimport SkewPolynomial_generic_dense
from sage.matrix.matrix_dense cimport Matrix_dense
from skew_polynomial_element cimport CenterSkewPolynomial_generic_dense

from polynomial_element cimport Polynomial
from sage.structure.element cimport RingElement


cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_generic_dense):
    # cache
    cdef list _conjugates
    cdef Polynomial _norm
    cdef _norm_factor
    cdef _optbound
    cdef dict _rdivisors
    cdef dict _types
    cdef _factorization

    cdef inline void _init_cache(self)

    # Karatsuba
    #cpdef RingElement _mul_karatsuba(self, RingElement other, cutoff=*)
    cpdef SkewPolynomial_finite_field_dense _mul_central(self, SkewPolynomial_finite_field_dense right)
    cpdef RingElement _mul_(self, RingElement right)
    cpdef rquo_rem_karatsuba(self, RingElement other, cutoff=*)

    cdef SkewPolynomial_finite_field_dense _rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef SkewPolynomial_finite_field_dense _coeff_llcm(self, SkewPolynomial_finite_field_dense other)

    # Inplace functions
    cdef void _inplace_lrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_pow_mod(self, Integer n, SkewPolynomial_finite_field_dense mod)
    cdef void _inplace_lmonic(self)
    cdef void _inplace_rmonic(self)
    cdef void _inplace_rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef Py_ssize_t _val_inplace_unit(self)
    cdef SkewPolynomial_finite_field_dense _rquo_inplace_rem(self, SkewPolynomial_finite_field_dense other)

    # Specific functions
    cdef Matrix_dense _matphir_c(self)
    cdef Matrix_dense _matmul_c(self)

    # Finding divisors
    cdef SkewPolynomial_finite_field_dense _rdivisor_c(P, CenterSkewPolynomial_generic_dense N)

    # Finding factorizations
    cdef _factor_c(self)
    cdef _factor_uniform_c(self)


cdef class SkewPolynomial_finite_field_karatsuba:
    cdef _parent
    cdef Py_ssize_t _order
    cdef Py_ssize_t _cutoff
    cdef RingElement _zero
    cdef _twist
    cdef char _algo_matrix
    cdef RingElement _t
    cdef Matrix_dense _T, _Tinv

    cdef list mul_step (self, list x, list y)
    cdef list mul_step_matrix(self, list x, list y)
    cdef list mul_iter(self, list x, list y, char flag)
    cdef list _twinv
    cdef list div_step(self, list a, Py_ssize_t ia, Py_ssize_t da, list b, Py_ssize_t ib, Py_ssize_t db)
    cdef list div_iter(self, list a, Py_ssize_t ia, Py_ssize_t da, list b, Py_ssize_t ib, Py_ssize_t db)
