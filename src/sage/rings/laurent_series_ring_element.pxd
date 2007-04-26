from sage.structure.element cimport AlgebraElement

cdef class LaurentSeries(AlgebraElement):
    cdef object __u
    cdef long __n
