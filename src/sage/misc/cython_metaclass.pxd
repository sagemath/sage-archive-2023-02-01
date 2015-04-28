cdef extern from "sage/misc/cython_metaclass.h":
    PyMethodDescr_CallSelf(desc, self)

cdef class TypeInitMetaclass(type):
    pass
