cdef class MemoryAllocator:
    r"""
    An object with for memory allocation, whose resources are freed upon __dealloc__.
    """
    def __cinit__(self):
        self.n = 0
        self.max_size = 0
        self.pointers = NULL

    cdef enlarge_if_needed(self):
        r"""
        Enlarge the list of local pointer.
        """
        if self.n < self.max_size:
            return
        cdef int    new_max_size = 1 if self.max_size == 0 else 2*self.max_size
        cdef void ** new_pointers = <void **> realloc(self.pointers, new_max_size*sizeof(void*))
        if new_pointers == NULL:
            raise MemoryError()
        self.pointers = new_pointers
        self.max_size = new_max_size

    cdef void * malloc(self, size_t size) except NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        self.enlarge_if_needed()
        cdef void * val = malloc(size)
        if val == NULL:
            raise MemoryError()
        self.pointers[self.n] = val
        self.n += 1
        return val

    def __dealloc__(self):
        r"""
        Free the allocated resources
        """
        cdef int i
        for i in range(self.n):
            free(self.pointers[i])
        free(self.pointers)
