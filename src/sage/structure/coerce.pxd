cdef class Coerce:
    cdef cmp_c(self, x, y)
    cdef bin_op_c(self, x, y, op)

