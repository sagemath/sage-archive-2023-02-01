include "../ext/cdefs.pxi"
import sage.ext.element
cimport sage.ext.element

cdef class IntegerMod(sage.ext.element.CommutativeRingElement):
    cdef public object _parent
    cdef mpz_t value
#    cdef int cmp(self, IntegerMod x)
    cdef void set_from_mpz(IntegerMod self, mpz_t value)
    cdef mpz_t* get_value(IntegerMod self)
    cdef mpz_t* mpz_modulus(IntegerMod self)
    cdef object _pari

