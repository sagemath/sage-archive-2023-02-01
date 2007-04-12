include "../ext/cdefs.pxi"

cdef extern from "stdint.h":
    ctypedef int int_fast32_t
    ctypedef int int_fast64_t
    int_fast32_t INTEGER_MOD_INT32_LIMIT
    int_fast64_t INTEGER_MOD_INT64_LIMIT

cimport sage.structure.element
from sage.rings.integer cimport Integer

cdef class NativeIntStruct:
    cdef Integer sageInteger
    cdef int_fast32_t int32
    cdef int_fast64_t int64
    cdef object table # how much faster is unsafe access to a list than a c array?
    cdef lookup(NativeIntStruct self, Py_ssize_t value)


cdef class IntegerMod_abstract(sage.structure.element.CommutativeRingElement):
    cdef NativeIntStruct __modulus
    cdef _new_c_from_long(self, long value)
    cdef void set_from_mpz(self, mpz_t value)
    cdef void set_from_long(self, long value)

cdef class IntegerMod_gmp(IntegerMod_abstract):
    cdef mpz_t value
    cdef mpz_t* get_value(IntegerMod_gmp self)
    cdef IntegerMod_gmp _new_c(self)

cdef class IntegerMod_int(IntegerMod_abstract):
    cdef public int_fast32_t ivalue
    cdef void set_from_int(IntegerMod_int self, int_fast32_t value)
    cdef int_fast32_t get_int_value(IntegerMod_int self)
    cdef IntegerMod_int _new_c(self, int_fast32_t value)

cdef class IntegerMod_int64(IntegerMod_abstract):
    cdef int_fast64_t ivalue
    cdef void set_from_int(IntegerMod_int64 self, int_fast64_t value)
    cdef int_fast64_t get_int_value(IntegerMod_int64 self)
    cdef IntegerMod_int64 _new_c(self, int_fast64_t value)
