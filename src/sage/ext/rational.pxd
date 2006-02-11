include "cdefs.pxi"
cimport element
import element

cdef class Rational(element.FieldElement):
    cdef mpq_t value
    cdef object _parent

    cdef cmp(Rational self, Rational x)
    cdef void set_from_mpq(Rational self, mpq_t value)

