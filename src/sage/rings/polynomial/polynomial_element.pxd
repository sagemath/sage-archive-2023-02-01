include "../../ext/interrupt.pxi"
include "../../ext/cdefs.pxi"
include '../../ext/stdsage.pxi'


from sage.structure.element import Element, CommutativeAlgebraElement
from sage.structure.element cimport Element, CommutativeAlgebraElement, ModuleElement
from polynomial_compiled import CompiledPolynomialFunction
from polynomial_compiled cimport CompiledPolynomialFunction

cdef class Polynomial(CommutativeAlgebraElement):
    cpdef ModuleElement _neg_(self)
    cdef char _is_gen
    cdef CompiledPolynomialFunction _compiled
    cpdef Polynomial truncate(self, long n)
    cdef long _hash_c(self)

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long n)

cdef class Polynomial_generic_dense(Polynomial):
    cdef list __coeffs
    cdef void __normalize(self)
#    cdef _dict_to_list(self, x, zero)

#cdef class Polynomial_generic_sparse(Polynomial):
#    cdef object __coeffs # a python dict (for now)
