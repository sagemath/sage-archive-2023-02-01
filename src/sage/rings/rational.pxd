include "../ext/cdefs.pxi"
cimport sage.structure.element
import  sage.structure.element

cdef class Rational(sage.structure.element.FieldElement):
    cdef mpq_t value
    cdef object _parent

    cdef cmp(Rational self, Rational x)
    cdef void set_from_mpq(Rational self, mpq_t value)

