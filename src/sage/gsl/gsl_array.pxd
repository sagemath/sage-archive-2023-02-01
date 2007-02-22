include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
include 'gsl.pxi'
cdef class GSLDoubleArray:
    cdef size_t n
    cdef size_t stride
    cdef double * data

#cdef class GSLDoubleComplexArray(GSLDoubleArray):
#    cdef size_t m
