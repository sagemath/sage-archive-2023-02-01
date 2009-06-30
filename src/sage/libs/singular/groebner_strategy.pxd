from sage.libs.singular.decl cimport skStrategy

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular, MPolynomial_libsingular
from sage.structure.sage_object cimport SageObject


cdef class GroebnerStrategy(SageObject):
    cdef skStrategy *_strat
    cdef MPolynomialRing_libsingular _parent
    cdef object _ideal

    cpdef MPolynomial_libsingular normal_form(self, MPolynomial_libsingular p)
