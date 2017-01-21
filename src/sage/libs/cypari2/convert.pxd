from .types cimport GEN
from .gen cimport Gen

cpdef integer_to_gen(x)
cpdef gen_to_integer(Gen x)
cpdef gen_to_python(Gen z)

cdef GEN gtoi(GEN g0) except NULL
cdef GEN PyLong_AsGEN(x)
cdef PyLong_FromGEN(GEN g)

cdef Gen new_t_POL_from_int_star(int* vals, unsigned long length, long varnum)
cdef Gen new_gen_from_double(double)
cdef Gen new_t_COMPLEX_from_double(double re, double im)
