include "../libs/ginac/decl.pxi"

cdef class PynacConstant:
    cdef GConstant* pointer
    cdef GConstant object
    cdef object _name
