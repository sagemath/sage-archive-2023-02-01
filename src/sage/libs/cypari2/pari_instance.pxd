from .types cimport *
cimport cython

from .gen cimport gen

cpdef long prec_bits_to_words(unsigned long prec_in_bits)
cpdef long prec_words_to_bits(long prec_in_words)
cpdef long default_bitprec()

cdef class PariInstance_auto:
    pass

cdef class PariInstance(PariInstance_auto):
    cdef readonly gen PARI_ZERO, PARI_ONE, PARI_TWO
    cpdef gen zero(self)
    cpdef gen one(self)
    cdef gen _empty_vector(self, long n)

cdef long get_var(v) except -2

# TODO: this should not be needed
cdef PariInstance _pari_instance
