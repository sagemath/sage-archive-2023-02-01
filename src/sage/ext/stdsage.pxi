cdef extern from "stdsage.h":
    # Memory management
    void  sage_free(void *p)
    void* sage_realloc(void *p, size_t n)
    void* sage_malloc(size_t)

