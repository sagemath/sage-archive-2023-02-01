from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class PowComputer_relative(PowComputer_class):
    cdef public object poly_ring
    cdef public object base_ring
    cdef public object modulus
    cdef unsigned long capdiv(self, unsigned long n)

cdef class PowComputer_relative_unram(PowComputer_relative):
    pass

cdef class PowComputer_relative_eis(PowComputer_relative):
    cdef public object pxe
