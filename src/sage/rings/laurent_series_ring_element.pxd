from sage.structure.element cimport ModuleElement, RingElement, AlgebraElement

cdef class LaurentSeries(AlgebraElement):
    cdef object __u
    cdef object __n
