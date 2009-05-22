include "../ext/cdefs.pxi"
include "../libs/ntl/decl.pxi"

import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport EuclideanDomainElement, RingElement, ModuleElement, Element
from sage.categories.morphism cimport Morphism

cdef class Integer(EuclideanDomainElement):
    cdef mpz_t value

    cdef void _to_ZZ(self, ZZ_c *z)
    cdef void set_from_mpz(self, mpz_t value)
    cdef mpz_t* get_value(self)
    cdef hash_c(self)

    cdef _pari_c(self)

    cdef _shift_helper(Integer self, y, int sign)
    cdef _lshift(self, long int n)
    cpdef _rshift_(Integer self, long int n)
    cdef _and(Integer self, Integer other)
    cdef _or(Integer self, Integer other)
    cdef _xor(Integer self, Integer other)

    cpdef size_t _exact_log_log2_iter(self,Integer m)
    cpdef size_t _exact_log_mpfi_log(self,m)
    cpdef RingElement _valuation(Integer self, Integer p)
    cdef object _val_unit(Integer self, Integer p)
    cdef Integer _divide_knowing_divisible_by(Integer self, Integer right)
    cdef bint _is_power_of(Integer self, Integer n)

    cdef _reduce_set(self, s) # do not use, since integers are immutable.

cdef class int_to_Z(Morphism):
    pass
