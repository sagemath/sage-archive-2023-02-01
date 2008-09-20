cdef class TripleDict:
    cdef buckets
    cdef double threshold
    cdef get(self, k1, k2, k3)
    cdef set(self, k1, k2, k3, value)

cdef class TripleDictIter:
    cdef TripleDict pairs
    cdef buckets, bucket_iter
