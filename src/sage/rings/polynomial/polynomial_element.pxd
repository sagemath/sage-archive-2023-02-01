include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"
include 'sage/ext/stdsage.pxi'


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
    cdef long _hash_c(self) except -1
    cpdef constant_coefficient(self)
    cpdef Polynomial _new_constant_poly(self, a, Parent P)

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long n)

cdef class Polynomial_generic_dense(Polynomial):
    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P)
    cdef list __coeffs
    cdef void __normalize(self)
#    cdef _dict_to_list(self, x, zero)

cpdef is_Polynomial(f)

#cdef class Polynomial_generic_sparse(Polynomial):
#    cdef object __coeffs # a python dict (for now)
