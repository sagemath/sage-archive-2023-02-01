cdef extern from "object.h":
    ctypedef class __builtin__.type [object PyHeapTypeObject]:
        pass

cdef class NestedClassMetaclass(type):
    pass
