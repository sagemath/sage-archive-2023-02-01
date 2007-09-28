from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_dense_mod_n(Polynomial):
    cdef object __poly

cdef class Polynomial_dense_modn_ntl_ZZ(Polynomial_dense_mod_n):
    pass

cdef class Polynomial_dense_modn_ntl_zz(Polynomial_dense_mod_n):
    pass

cdef class Polynomial_dense_mod_p(Polynomial_dense_mod_n):
    pass
