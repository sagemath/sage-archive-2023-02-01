cdef class TimeSeries:
    cdef double* _values
    cdef Py_ssize_t _length
