cdef class OrderedPartition:
    cdef public int length
    cdef int** data
    cdef int* sizes

    cdef public int size(OrderedPartition self, int n)

