from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense_inexact

cdef class PowComputer_relative(PowComputer_class):
    cdef public object poly_ring
    cdef public object base_ring
    cdef public Polynomial_generic_dense_inexact modulus
    cdef Polynomial_generic_dense_inexact powhelper_oneunit
    cdef Polynomial_generic_dense_inexact powhelper_teichdiff
    cdef Polynomial_generic_dense_inexact powhelper_cconv_out
    cdef unsigned long capdiv(self, unsigned long n)

cdef class PowComputer_relative_unram(PowComputer_relative):
    pass

cdef class PowComputer_relative_eis(PowComputer_relative):
    cdef public Polynomial_generic_dense_inexact pxe
    cdef Polynomial_generic_dense_inexact poly_clist
