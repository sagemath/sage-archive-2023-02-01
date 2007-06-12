include "singular-cdefs.pxi"

from sage.rings.rational cimport Rational
from sage.structure.element cimport Element
from sage.rings.integer cimport Integer
from sage.rings.finite_field_givaro cimport FiniteField_givaro,FiniteField_givaroElement

cdef class Conversion:
    cdef extern Rational si2sa_QQ(self, number (*),ring (*))
    cdef extern FiniteField_givaroElement si2sa_GFqGivaro(self, number *n, ring *_ring, FiniteField_givaro base)
    cdef extern object si2sa_GFqPari(self, number *n, ring *_ring, object base)

    cdef extern number *sa2si_QQ(self, Rational ,ring (*))
    cdef extern number *sa2si_ZZ(self, Integer d, ring *_ring)

    cdef extern number *sa2si_GFqGivaro(self, int exp ,ring (*))
    cdef extern number *sa2si_GFqPari(self, object vector, ring *_ring)

    cdef extern object si2sa(self, number *n, ring *_ring, object base)
    cdef extern number *sa2si(self, Element elem, ring * _ring)
