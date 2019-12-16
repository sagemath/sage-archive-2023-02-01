from sage.structure.element cimport AlgebraElement, ModuleElement
from sage.rings.laurent_series_ring_element cimport LaurentSeries

cdef class PuiseuxSeries(AlgebraElement):
     cpdef ModuleElement _l
     cdef long _e
