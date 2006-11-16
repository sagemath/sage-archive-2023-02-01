include "../ext/cdefs.pxi"

import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport EuclideanDomainElement, RingElement, ModuleElement

cdef class Integer(EuclideanDomainElement):
    cdef mpz_t value

    cdef void set_from_mpz(self, mpz_t value)
    cdef mpz_t* get_value(self)
    cdef object _pari

    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef ModuleElement _neg_c_impl(self)

cdef extern int set_mpz(Integer self, mpz_t value)
