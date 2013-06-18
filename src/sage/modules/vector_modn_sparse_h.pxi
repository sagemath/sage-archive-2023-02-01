cdef struct c_vector_modint:
    int *entries
    int p
    Py_ssize_t *positions
    Py_ssize_t degree
    Py_ssize_t num_nonzero

