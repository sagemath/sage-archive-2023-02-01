# I would like to use macros to define this to be just under sqrt(INT_MAX)
cdef enum:
    INTEGER_MOD_INT_LIMIT = 46340

include "../ext/cdefs.pxi"
import sage.ext.element
cimport sage.ext.element

cdef class IntegerMod(sage.ext.element.CommutativeRingElement):
    cdef public object _parent
    cdef mpz_t value
    cdef void set_from_mpz(IntegerMod self, mpz_t value)
    cdef mpz_t* get_value(IntegerMod self)
    cdef mpz_t* mpz_modulus(IntegerMod self)


cdef class IntegerMod_int(IntegerMod):
    cdef int imodulus  # one extra word of storage is better than mpz -> int all over the place, maybe if integer_mod_ring was in pyrex we could eventually store it there
    cdef int ivalue
    cdef void set_from_mpz(IntegerMod_int self, mpz_t value)
    cdef void set_from_int(IntegerMod_int self, int value)
    cdef int get_int_value(IntegerMod_int self)

