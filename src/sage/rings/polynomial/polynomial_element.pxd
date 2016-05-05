from sage.structure.element import Element, CommutativeAlgebraElement
from sage.structure.element cimport Element, CommutativeAlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from polynomial_compiled import CompiledPolynomialFunction
from polynomial_compiled cimport CompiledPolynomialFunction

cdef class Polynomial(CommutativeAlgebraElement):
    cpdef ModuleElement _neg_(self)
    cdef char _is_gen
    cdef CompiledPolynomialFunction _compiled
    cpdef Polynomial truncate(self, long n)
    cpdef Polynomial inverse_series_trunc(self, long prec)
    cdef long _hash_c(self) except -1
    cpdef constant_coefficient(self)
    cpdef Polynomial _new_constant_poly(self, a, Parent P)

    cpdef bint is_zero(self)
    cpdef bint is_one(self)

    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n)

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long n)

    cdef get_unsafe(self, Py_ssize_t i)

cdef class Polynomial_generic_dense(Polynomial):
    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P)
    cdef list __coeffs
    cdef int __normalize(self) except -1

cdef class Polynomial_generic_dense_inexact(Polynomial_generic_dense):
    pass

cpdef is_Polynomial(f)
