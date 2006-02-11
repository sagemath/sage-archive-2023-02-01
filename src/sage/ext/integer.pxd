include "cdefs.pxi"
cimport element
import element

cdef class Integer(element.EuclideanDomainElement):
    cdef mpz_t value
    cdef int cmp(self, Integer x)
    cdef void set_from_mpz(Integer self, mpz_t value)
    cdef mpz_t* get_value(self)
    cdef object _pari

