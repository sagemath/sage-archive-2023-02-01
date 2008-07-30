
cdef double* to_double_array(v) except NULL:
    cdef double x
    cdef double* w = <double*> safe_malloc(sizeof(double)*len(v))
    cdef Py_ssize_t i = 0
    for x in v:
        w[i] = x
        i += 1
    return w

cdef int* to_int_array(v) except NULL:
    cdef int x
    cdef int* w = <int*> safe_malloc(sizeof(int)*len(v))
    cdef Py_ssize_t i = 0
    for x in v:
        w[i] = x
        i += 1
    return w

cdef void* safe_malloc(int bytes) except NULL:
    """
    malloc the given bytes of memory and check that the malloc
    succeeds -- if not raise a MemoryError.
    """
    cdef void* t = sage_malloc(bytes)
    if not t:
        raise MemoryError, "error allocating memory for Hidden Markov Model"
    return t

