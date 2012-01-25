cdef class MonoDict:
    cdef __weakref__
    cdef Py_ssize_t _size
    cdef buckets
    cdef dict _refcache
    cdef double threshold
    cdef public MonoDictEraser eraser

cdef class MonoDictEraser:
    cdef object D

cdef class TripleDict(MonoDict):
    cdef get(self, object k1, object k2, object k3)
    cdef set(self, object k1, object k2, object k3, value)

cdef class TripleDictEraser(MonoDictEraser):
    pass
