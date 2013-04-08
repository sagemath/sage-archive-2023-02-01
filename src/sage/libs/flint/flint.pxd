cdef extern from "flint.h":

    cdef long FLINT_BITS
    cdef long FLINT_D_BITS

    cdef unsigned long FLINT_BIT_COUNT(unsigned long)
    void flint_stack_cleanup()
