from cysignals.memory cimport *

cdef class MemoryAllocator:
    r"""
    An object for memory allocation, whose resources are freed upon
    ``__dealloc__``.

    EXAMPLES::

        sage: cython(
        ....: '''
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: cdef MemoryAllocator mem = MemoryAllocator()
        ....: cdef void * ptr
        ....: for n in range(100):
        ....:     ptr = mem.malloc(n)
        ....:     mem.realloc(ptr, 2*n)
        ....:     mem.calloc(n, n)
        ....:     ptr = mem.allocarray(n, n)
        ....:     mem.reallocarray(ptr, n + 1, n)
        ....:     mem.aligned_malloc(32, (n//32 + 1)*32)
        ....:     mem.aligned_calloc(16, n, 16)
        ....:     mem.aligned_allocarray(8, n, 8)
        ....: ''')
    """
    def __cinit__(self):
        """
        EXAMPLES::

            sage: cython(
            ....: '''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: cdef MemoryAllocator mem = MemoryAllocator.__new__(MemoryAllocator)
            ....: mem.malloc(10000)
            ....: print(mem.n)
            ....: print(mem.size)
            ....: ''')
            1
            16
        """
        self.n = 0
        self.size = 16
        self.pointers = self.static_pointers

    cdef int resize(self, size_t new_size) except -1:
        r"""
        Resize the list of pointers to contain ``new_size`` elements.

        It is required that ``new_size`` is at least ``self.n``, but
        this condition is not checked.
        """
        cdef size_t i
        if self.pointers == self.static_pointers:
            # Case 1: allocate pointers for the first time
            self.pointers = <void **>check_allocarray(new_size, sizeof(void*))
            for i in range(self.n):
                self.pointers[i] = self.static_pointers[i]
        else:
            # Case 2: resize pointers
            self.pointers = <void **>check_reallocarray(self.pointers, new_size, sizeof(void*))
        self.size = new_size

    cdef void * malloc(self, size_t size) except? NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        return self.aligned_malloc(1, size)

    cdef void * realloc(self, void * ptr, size_t size) except? NULL:
        r"""
        Re-allocates `ptr` and automatically frees it later.

        TESTS::

            sage: cython('''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: cdef MemoryAllocator mem = MemoryAllocator()
            ....: cdef MemoryAllocator mem2 = MemoryAllocator()
            ....: cdef void * ptr = mem.malloc(20)
            ....: cdef void * ptr2 = mem2.malloc(20)
            ....: try:
            ....:     mem.realloc(ptr2, 21)
            ....: except:
            ....:     print('mem knows not this memory')
            ....: ''')
            mem knows not this memory
        """
        if unlikely(ptr is NULL):
            return self.malloc(size)
        cdef void * val
        cdef size_t i = self.find_pointer(ptr) - 1
        # find_pointer returns one too much to have a unique except value
        val = check_realloc(ptr, size)
        self.pointers[i] = val
        return val

    cdef void * calloc(self, size_t nmemb, size_t size) except? NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        return self.aligned_calloc(1, nmemb, size)

    cdef void * allocarray(self, size_t nmemb, size_t size) except? NULL:
        r"""
        Returns a new pointer and stores it to be automatically freed later.
        """
        return self.aligned_allocarray(1, nmemb, size)

    cdef void * reallocarray(self, void * ptr, size_t nmemb,
                             size_t size) except? NULL:
        r"""
        Re-allocates `ptr` and automatically frees it later.
        """
        if unlikely(ptr is NULL):
            return self.allocarray(nmemb, size)
        cdef void * val
        cdef size_t i = self.find_pointer(ptr) - 1
        # find_pointer returns one too much to have a unique except value
        val = check_reallocarray(ptr, nmemb, size)
        self.pointers[i] = val
        return val

    cdef inline void * aligned_malloc(self, size_t alignment,
                                      size_t size) except? NULL:
        r"""
        Returns new aligned pointer. Stores it to be automatically freed later.

        Alignment must be a power of two.

        .. NOTE::

            If you want to allocate multiple (small) aligned arrays with the
            same alignment and really want to be memory efficient, you can
            allocate one large aligned array instead.

        TESTS::

            sage: cython('''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: cdef MemoryAllocator mem = MemoryAllocator()
            ....: cdef void * ptr
            ....: for i in range(12):
            ....:     ptr = mem.aligned_malloc(2**i, 4048)
            ....:     if <size_t> ptr != (<size_t> ptr) & ~(2**i-1):
            ....:         print('alignment not working')
            ....: ''')
        """
        self.enlarge_if_needed()
        if unlikely(size == 0):
            return NULL
        cdef void * val = check_malloc(size + alignment - 1)
        cdef void * ptr = <void *> ((<size_t> val + alignment-1) & ~(alignment-1))
        # setting ptr to be aligned according to the users wishes
        self.pointers[self.n] = val
        self.n += 1
        return ptr

    cdef inline void * aligned_calloc(self, size_t alignment, size_t nmemb,
                                      size_t size) except? NULL:
        r"""
        Returns new aligned pointer. Stores it to be automatically freed later.

        Alignment must be a power of two.

        .. NOTE::

            If you want to allocate multiple (small) aligned arrays with the
            same alignment and really want to be memory efficient, you can
            allocate one large aligned array instead.

        TESTS::

            sage: cython('''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: cdef MemoryAllocator mem = MemoryAllocator()
            ....: cdef void * ptr
            ....: for i in range(12):
            ....:     ptr = mem.aligned_calloc(2**i, i, 2**i)
            ....:     if <size_t> ptr != (<size_t> ptr) & ~(2**i-1):
            ....:         print('alignment not working')
            ....: ''')
        """
        self.enlarge_if_needed()
        if unlikely(size == 0 or nmemb == 0):
            return NULL
        cdef void * val = check_calloc(nmemb + (alignment+size-2)//size, size)
        cdef void * ptr = <void *> ((<size_t> val + alignment-1) & ~(alignment-1))
        # setting ptr to be aligned according to the users wishes
        self.pointers[self.n] = val
        self.n += 1
        return ptr

    cdef inline void * aligned_allocarray(self, size_t alignment, size_t nmemb,
                                          size_t size) except? NULL:
        r"""
        Returns new aligned pointer. Stores it to be automatically freed later.

        Alignment must be a power of two.

        .. NOTE::

            If you want to allocate multiple (small) aligned arrays with the
            same alignment and really want to be memory efficient, you can
            allocate one large aligned array instead.

        TESTS::

            sage: cython('''
            ....: from sage.ext.memory_allocator cimport MemoryAllocator
            ....: cdef MemoryAllocator mem = MemoryAllocator()
            ....: cdef void * ptr
            ....: for i in range(12):
            ....:     ptr = mem.aligned_allocarray(2**i, i, 2**i)
            ....:     if <size_t> ptr != (<size_t> ptr) & ~(2**i-1):
            ....:         print('alignment not working')
            ....: ''')
        """
        self.enlarge_if_needed()
        if unlikely(size == 0 or nmemb == 0):
            return NULL
        cdef void * val = check_allocarray(nmemb + (alignment+size-2)//size, size)
        cdef void * ptr = <void *> ((<size_t> val + alignment-1) & ~(alignment-1))
        # setting ptr to be aligned according to the users wishes
        self.pointers[self.n] = val
        self.n += 1
        return ptr

    def __dealloc__(self):
        r"""
        Free the allocated resources

        EXAMPLES::

            sage: from sage.ext.memory_allocator import MemoryAllocator
            sage: _ = MemoryAllocator()
        """
        cdef size_t i
        for i in range(self.n):
            sig_free(self.pointers[i])
        if self.pointers != self.static_pointers:
            sig_free(self.pointers)
