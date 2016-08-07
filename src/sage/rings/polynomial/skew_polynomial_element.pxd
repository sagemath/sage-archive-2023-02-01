from sage.structure.element cimport AlgebraElement
from sage.structure.parent cimport Parent
from sage.rings.morphism cimport Morphism
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

cdef class SkewPolynomial(AlgebraElement):
    cdef _is_gen

    cdef long _hash_c(self)
    cdef SkewPolynomial _new_c(self,list coeffs,Parent P,char check=*)
    cpdef SkewPolynomial _new_constant_poly(self,RingElement a,Parent P,char check=*)
    cpdef _neg_(self)
    cpdef list list(self)

    cdef list _coeffs
    cdef void __normalize(self)
    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right)
    cdef void _inplace_pow(self, Py_ssize_t n)

cdef class SkewPolynomial_generic_dense(SkewPolynomial):
    cpdef power_right_mod(self,right,modulus)
    cpdef power_left_mod(self,right,modulus)

cdef class SkewPolynomialBaseringInjection(Morphism):
    cdef RingElement _an_element
    cdef object _new_constant_poly_
