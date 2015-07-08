include 'sage/ext/interrupt.pxi'

cdef class MemoryAllocator:
    r"""
    An object with for memory allocation, whose resources are freed upon __dealloc__.
    """
    def __cinit__(self):
        self.n = 0
        self.max_size = 0
        self.pointers = NULL

    cdef int enlarge_if_needed(self) except -1:
        r"""
        Enlarge the list of local pointer.
        """
        if self.n < self.max_size:
            return
        cdef int     new_max_size = 2*self.max_size if self.max_size else 1
        cdef void ** new_pointers = <void **> check_realloc(self.pointers, new_max_size*sizeof(void*))
        self.pointers = new_pointers
        self.max_size = new_max_size

    cdef void * malloc(self, size_t size) except NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        self.enlarge_if_needed()
        cdef void * val = check_malloc(size)
        self.pointers[self.n] = val
        self.n += 1
        return val

    cdef void * calloc(self, size_t nmemb, size_t size) except NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        self.enlarge_if_needed()
        cdef void * val = check_calloc(nmemb, size)
        self.pointers[self.n] = val
        self.n += 1
        return val

    def __dealloc__(self):
        r"""
        Free the allocated resources
        """
        cdef size_t i
        for i in range(self.n):
            free(self.pointers[i])
        free(self.pointers)
