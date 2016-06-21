#include "../../ext/interrupt.pxi"
include "cysignals/signals.pxi"
include "../../ext/cdefs.pxi"
include '../../ext/stdsage.pxi'


from sage.rings.integer cimport Integer

from sage.structure.element import Element, AlgebraElement
from sage.structure.element cimport Element, AlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from polynomial_compiled import CompiledPolynomialFunction
from polynomial_compiled cimport CompiledPolynomialFunction

from sage.rings.morphism cimport RingHomomorphism
from sage.structure.element cimport RingElement

from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense


cdef class CenterSkewPolynomial_generic_dense(Polynomial_generic_dense):
    pass


cdef class SkewPolynomial(AlgebraElement):
    cdef _is_gen

    cdef long _hash_c(self)
    cdef list _list_c(self)
    cdef SkewPolynomial _new_c(self,list coeffs,Parent P,char check=*)
    cpdef SkewPolynomial _new_constant_poly(self,RingElement a,Parent P,char check=*)
    cpdef ModuleElement _neg_(self)



cdef class SkewPolynomial_generic_dense(SkewPolynomial):
    cdef list __coeffs
    cdef void __normalize(self)
    cpdef _pow_(self,right,modulus=*,side=*)

    # Inplace functions
    cdef void _inplace_add(self, SkewPolynomial_generic_dense right)
    cdef void _inplace_sub(self, SkewPolynomial_generic_dense right)
    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right)
    cdef void _inplace_lmul(self, SkewPolynomial_generic_dense right)
    cdef void _inplace_pow(self, Py_ssize_t n)


cdef class SkewPolynomialBaseringInjection(RingHomomorphism):
    cdef RingElement _an_element
    cdef object _new_constant_poly_
