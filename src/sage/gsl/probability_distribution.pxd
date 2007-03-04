include '../ext/cdefs.pxi'
include 'gsl.pxi'
cdef class ProbabilityDistribution:
    pass

cdef class SphericalDistribution(ProbabilityDistribution):
    cdef gsl_rng_type *T
    cdef gsl_rng*r
    cdef Py_ssize_t dimension
    cdef long int seed
    cdef double* vec

cdef class RealDistribution(ProbabilityDistribution):
    cdef gsl_rng_type *T
    cdef gsl_rng *r
    cdef int distribution_type
    cdef double* parameters
    cdef long int seed
    cdef object name
#    cdef double (*generator_1)(gsl_rng*)
#    cdef double (*generator_2)(gsl_rng*,double)
#    cdef _get_random_element_c(self)
