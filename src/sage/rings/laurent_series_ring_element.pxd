from sage.structure.element cimport AlgebraElement, ModuleElement

cdef class LaurentSeries(AlgebraElement):
    cdef ModuleElement __u
    cdef long __n
