from sage.structure.sage_object cimport SageObject

cdef class SFunction(SageObject):
    cdef unsigned int _serial
    cdef object _name
    cdef object _ginac_name
    cdef object _latex_name
    cdef int _nargs
    cdef dict __dict__

    # cache hash value
    cdef long _hash_(self)
    cdef bint __hinit
    cdef long __hcache

    # common initialization from __init__ and __setstate__
    cdef _init_(self)

    cdef object _eval_
    cdef object _evalf_
    cdef object _conjugate_
    cdef object _real_part_
    cdef object _imag_part_
    cdef object _derivative_
    cdef object _power_
    cdef object _series_
    cdef object _print_
    cdef object _print_latex_
