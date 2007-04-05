include "singular-cdefs.pxi"

from sage.rings.rational cimport Rational

cdef class Conversion:
    cdef extern Rational si2sa_QQ(self, number (*),ring (*))
    cdef extern number *sa2si_QQ(self, Rational ,ring (*))

    cdef extern object si2sa(self, number *n, ring *_ring, object base)
