include "../ext/cdefs.pxi"
cimport sage.structure.element
import  sage.structure.element

cimport integer

cdef class Rational(sage.structure.element.FieldElement):
    cdef mpq_t value

    cdef void set_from_mpq(Rational self, mpq_t value)
    cdef _lshift(self, unsigned long int exp)
    cdef _rshift(self, unsigned long int exp)

    cdef integer.Integer _integer_c(self)
