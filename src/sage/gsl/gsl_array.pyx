#include 'gsl.pxi'
cdef class GSLDoubleArray:
    def __init__(self, size_t n, size_t stride = 1, data = None):
        cdef int i

        self.n = n
        self.stride = stride
        self.data = <double *> PyMem_Malloc(sizeof(double)*n)
        if data is not None:
            for i from 0 <= i < n:
                self.data[i] = data[i]
        else:
            for i from 0 <= i < n:
                self.data[i] = 0

    def  __dealloc__(self):
        PyMem_Free(self.data)

    def __len__(self):
        return self.n

    def __getslice__(self, i, j):
        # Todo -- make this actually fast.
        return list(self)[i:j]

    def __repr__(self):
        return str(list(self))

    def __setitem__(self, size_t i, x):
        # just set real for now
        if i < 0 or i >= self.n:
            raise IndexError
        self.data[i] = x

    def __getitem__(self, size_t i):
        if i < 0 or i >= self.n:
            raise IndexError
        return self.data[i]
