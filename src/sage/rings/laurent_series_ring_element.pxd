from sage.structure.element cimport AlgebraElement, ModuleElement
from power_series_ring_element cimport PowerSeries

cdef class LaurentSeries(AlgebraElement):
    cdef ModuleElement __u
    cdef long __n
