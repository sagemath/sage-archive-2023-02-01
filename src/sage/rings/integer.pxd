include "../ext/cdefs.pxi"
include "../libs/ntl/decl.pxi"

import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport EuclideanDomainElement, RingElement, ModuleElement

cdef class Integer(EuclideanDomainElement):
    cdef mpz_t value

    cdef void _to_ZZ(self, ntl_c_ZZ *z)
    cdef void set_from_mpz(self, mpz_t value)
    cdef mpz_t* get_value(self)
    #cdef object _pari

    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef ModuleElement _neg_c_impl(self)

    cdef _lshift(self, long int n)
    cdef _rshift(Integer self, long int n)
    cdef _and(Integer self, Integer other)
    cdef _or(Integer self, Integer other)

    cdef RingElement _valuation(Integer self, Integer p)
    cdef object _val_unit(Integer self, Integer p)
    cdef Integer _divide_knowing_divisible_by(Integer self, Integer right)
    cdef bint _is_power_of(Integer self, Integer n)

cdef extern int set_mpz(Integer self, mpz_t value)
