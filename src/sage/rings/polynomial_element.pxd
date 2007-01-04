include "../ext/interrupt.pxi"
include "../ext/cdefs.pxi"
include '../ext/stdsage.pxi'


from sage.structure.element import Element, IntegralDomainElement, CommutativeAlgebraElement
from sage.structure.element cimport Element, IntegralDomainElement, CommutativeAlgebraElement


cdef class Polynomial(CommutativeAlgebraElement):
    cdef size_t degree
    cdef char _is_gen

cdef class Polynomial_generic_dense(Polynomial):
    cdef object __coeffs # a python list

cdef class Polynomial_generic_sparse(Polynomial):
    cdef object __coeffs # a python dict (for now)
