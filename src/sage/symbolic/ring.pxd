from sage.rings.ring cimport CommutativeRing

cdef class SymbolicRing(CommutativeRing):
    cdef public dict symbols
