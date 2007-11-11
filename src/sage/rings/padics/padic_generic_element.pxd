include "../../ext/cdefs.pxi"

cimport sage.structure.element
from sage.structure.element cimport Element
from sage.rings.padics.local_generic_element cimport LocalGenericElement
#from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

cdef class pAdicGenericElement(LocalGenericElement):
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef long valuation_c(self)
    cdef val_unit_c(self)
    cdef int _set_from_Integer(self, Integer x, absprec, relprec) except -1
    cdef int _set_from_Rational(self, Rational x, absprec, relprec) except -1

cdef extern void teichmuller_set_c(mpz_t value, mpz_t p, mpz_t ppow)
