# I would like to use macros to define this to be just under sqrt(INT_MAX)
cdef enum:
    INTEGER_MOD_INT_LIMIT = 46340

include "../ext/cdefs.pxi"
import sage.ext.element
cimport sage.ext.element

#import sage.ext.integer
#cimport sage.ext.integer

cdef class IntegerMod_abstract(sage.ext.element.CommutativeRingElement):
    cdef public object _parent


cdef class IntegerMod_gmp(IntegerMod_abstract):
    cdef mpz_t value
#    cdef void set_from_Integer(IntegerMod_gmp self, sage.ext.integer.Integer value)
    cdef void set_from_mpz(IntegerMod_gmp self, mpz_t value)
    cdef mpz_t* get_value(IntegerMod_gmp self)
    cdef mpz_t* mpz_modulus(IntegerMod_gmp self)


cdef class IntegerMod_int(IntegerMod_abstract):
    cdef int imodulus  # one extra word of storage is better than mpz -> int all over the place, maybe if integer_mod_ring was in pyrex we could eventually store it there
    cdef int ivalue
#    cdef void set_from_Integer(IntegerMod_int self, sage.ext.integer.Integer value)
    cdef void set_from_mpz(IntegerMod_int self, mpz_t value)
    cdef void set_from_int(IntegerMod_int self, int value)
    cdef int get_int_value(IntegerMod_int self)

