#
# Declare C library functions used in Sage
#

include "python.pxi"

from libc.stdio cimport *
from libc.string cimport strlen, strcpy, memset, memcpy

from libc.math cimport sqrt
# Cython misdeclares these: http://trac.cython.org/cython_trac/ticket/801
cdef extern from "<math.h>":
    double frexp(double x, int* exponent)
    double ldexp(double x, int exponent)


from sage.libs.gmp.all cimport *
cdef extern from "<gmp.h>":
    pass  # Cython bug sometimes includes this in the wrong place
