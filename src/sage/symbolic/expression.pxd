include "../libs/ginac/decl.pxi"

from sage.structure.element cimport CommutativeRingElement


cdef class Symbol(CommutativeRingElement):
    cdef GSymbol _gobj

cdef class Expression(CommutativeRingElement):
    cdef GEx _gobj

cdef Expression new_Expression_from_GEx(GEx juice)
