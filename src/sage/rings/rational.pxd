include "../ext/cdefs.pxi"
cimport sage.structure.element
import  sage.structure.element

cimport integer

cpdef rational_power_parts(a, b, factor_limit=?)


cdef class Rational(sage.structure.element.FieldElement):
    cdef mpq_t value

    cdef __set_value(self, x, unsigned int base)
    cdef void set_from_mpq(Rational self, mpq_t value)
    cdef _lshift(self, long int exp)
    cdef _rshift(self, long int exp)

    cdef _val_unit(self, integer.Integer p)
