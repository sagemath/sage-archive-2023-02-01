cdef extern from *:
    ctypedef double* const_double_ptr "const double*"
    cdef int BINARY = 1
    cdef int REAL = -1
    cdef int INTEGER = 0

from sage.numerical.backends.generic_backend cimport GenericBackend

cdef class MIPVariable
cdef class MixedIntegerLinearProgram



cdef class MixedIntegerLinearProgram:
    cdef GenericBackend _backend
    cdef list _mipvariables
    cdef MIPVariable _default_mipvariable
    cdef dict _variables
    cdef int __BINARY
    cdef int __REAL
    cdef int __INTEGER
    cpdef int number_of_constraints(self)

cdef class MIPVariable:
    cdef MixedIntegerLinearProgram _p
    cdef int _dim
    cdef dict _dict
    cdef int _vtype
    cdef char * _name
    cdef bint _hasname
