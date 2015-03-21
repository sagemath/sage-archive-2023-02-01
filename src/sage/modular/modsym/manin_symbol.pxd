from sage.structure.element cimport Element

cdef class ManinSymbol(Element):
    cdef tuple __t
