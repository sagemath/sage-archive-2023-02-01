"""
GSL arrays
"""

include 'sage/ext/stdsage.pxi'

cdef class GSLDoubleArray:
    r"""
    EXAMPLES::

        sage: a = WaveletTransform(128,'daubechies',4)
        sage: for i in range(1, 11):
        ....:     a[i] = 1
        sage: a[:6:2]
        [0.0, 1.0, 1.0]
    """
    def __init__(self, size_t n, size_t stride = 1, data = None):
        cdef int i

        self.n = n
        self.stride = stride
        self.data = <double *> sage_malloc(sizeof(double)*n)
        if data is not None:
            for i from 0 <= i < n:
                self.data[i] = data[i]
        else:
            for i from 0 <= i < n:
                self.data[i] = 0

    def  __dealloc__(self):
        sage_free(self.data)

    def __len__(self):
        return self.n

    def __repr__(self):
        return str(list(self))

    def __setitem__(self, size_t i, x):
        # just set real for now
        if i < 0 or i >= self.n:
            raise IndexError
        self.data[i] = x

    def __getitem__(self, i):
        if isinstance(i, slice):
            start, stop, step = i.indices(len(self))
            # TODO -- make this actually fast.
            return list(self)[start:stop:step]
        else:
            if i < 0 or i >= self.n:
                raise IndexError
            return self.data[i]
