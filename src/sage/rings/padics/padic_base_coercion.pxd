include "../../ext/cdefs.pxi"

from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.rings.morphism cimport RingHomomorphism_coercion, RingMap
from sage.rings.padics.pow_computer cimport PowComputer_base

cdef class pAdicCoercion_ZZ_CR(RingHomomorphism_coercion):
    cdef PowComputer_base prime_pow
    cdef RingMap _section
cdef class pAdicConvert_CR_ZZ(RingMap):
    pass
cdef class pAdicCoercion_QQ_CR(RingHomomorphism_coercion):
    cdef PowComputer_base prime_pow
    cdef RingMap _section
cdef class pAdicConvert_CR_QQ(RingMap):
    pass
cdef class pAdicConvert_QQ_CR(Morphism):
    cdef PowComputer_base prime_pow
    cdef RingMap _section

cdef class pAdicCoercion_ZZ_CA(RingHomomorphism_coercion):
    cdef PowComputer_base prime_pow
    cdef RingMap _section
cdef class pAdicConvert_CA_ZZ(RingMap):
    pass
cdef class pAdicConvert_QQ_CA(Morphism):
    cdef PowComputer_base prime_pow
    cdef RingMap _section

cdef class pAdicCoercion_ZZ_FM(RingHomomorphism_coercion):
    cdef PowComputer_base prime_pow
    cdef RingMap _section
cdef class pAdicConvert_FM_ZZ(RingMap):
    pass
cdef class pAdicConvert_QQ_FM(Morphism):
    cdef PowComputer_base prime_pow
    cdef RingMap _section
