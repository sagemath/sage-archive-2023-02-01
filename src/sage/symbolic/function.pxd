
cdef class SFunction:
    cdef unsigned int serial
    cdef object name
    cdef int nargs

    cdef object eval_f
    cdef object evalf_f
    cdef object conjugate_f
    cdef object real_part_f
    cdef object imag_part_f
    cdef object derivative_f
    cdef object power_f
    cdef object series_f
    cdef object print_f
    cdef object print_latex_f
    cdef object latex_name

    # cache hash value
    cdef long _hash_(self)
    cdef bint __hinit
    cdef long __hcache

    # common initialization from __init__ and __setstate__
    cdef _init_(self)
