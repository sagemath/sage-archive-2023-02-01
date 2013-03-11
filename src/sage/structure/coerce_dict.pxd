cdef class MonoDict:
    cdef __weakref__
    cdef Py_ssize_t _size
    cdef buckets
    cdef bint weak_values
    cdef double threshold
    cdef public MonoDictEraser eraser
    cdef get(self, object k)
    cdef set(self, object k, value)

cdef class MonoDictEraser:
    cdef object D

cdef class TripleDict:
    cdef __weakref__
    cdef Py_ssize_t _size
    cdef buckets
    cdef bint weak_values
    cdef double threshold
    cdef public TripleDictEraser eraser
    cdef get(self, object k1, object k2, object k3)
    cdef set(self, object k1, object k2, object k3, value)

cdef class TripleDictEraser:
    cdef object D

