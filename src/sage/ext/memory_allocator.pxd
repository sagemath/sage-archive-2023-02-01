cimport cython

cdef extern from *:
    int unlikely(int) nogil  #Defined by Cython

@cython.final
cdef class MemoryAllocator:
    cdef size_t n
    cdef size_t size
    cdef void ** pointers
    cdef void * static_pointers[16]  # If n <= 16, store pointers here

    cdef void * malloc(self, size_t size) except? NULL
    cdef void * realloc(self, void* ptr, size_t size) except? NULL
    cdef void * calloc(self, size_t nmemb, size_t size) except? NULL
    cdef void * allocarray(self, size_t nmemb, size_t size) except? NULL
    cdef void * reallocarray(self, void * ptr, size_t nmemb,
                             size_t size) except? NULL
    cdef inline void * aligned_malloc(self, size_t alignment,
                                      size_t size) except? NULL
    cdef inline void * aligned_calloc(self, size_t alignment, size_t nmemb,
                                      size_t size) except? NULL
    cdef inline void * aligned_allocarray(self, size_t alignment, size_t nmemb,
                                          size_t size) except? NULL
    cdef int resize(self, size_t new_size) except -1
    cdef inline int enlarge_if_needed(self) except -1:
        r"""
        Enlarge the list of pointers if needed such that there is at
        least one free entry.
        """
        if unlikely(self.n >= self.size):
            return self.resize(self.size * 2)

    cdef inline size_t find_pointer(self, void * ptr) except 0:
        r"""
        Returns the index of ptr plus 1.
        """
        cdef size_t i = 0
        while (i < self.n and self.pointers[i] != ptr):
            i += 1
        if unlikely(i == self.n):
            raise ValueError("given pointer not found in MemoryAllocator")
        return i + 1
