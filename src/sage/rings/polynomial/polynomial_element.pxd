include "../../ext/interrupt.pxi"
include "../../ext/cdefs.pxi"
include '../../ext/stdsage.pxi'


from sage.structure.element import Element, CommutativeAlgebraElement
from sage.structure.element cimport Element, CommutativeAlgebraElement, ModuleElement
from polynomial_compiled import CompiledPolynomialFunction
from polynomial_compiled cimport CompiledPolynomialFunction

cdef class Polynomial(CommutativeAlgebraElement):
    cdef ModuleElement _neg_c_impl(self)
    cdef char _is_gen
    cdef CompiledPolynomialFunction _compiled

cdef class Polynomial_generic_dense(Polynomial):
    cdef object __coeffs # a python list
    cdef void __normalize(self)

#cdef class Polynomial_generic_sparse(Polynomial):
#    cdef object __coeffs # a python dict (for now)
