from sage.structure.element cimport ModuleElement
from power_series_ring_element cimport PowerSeries

cdef class PowerSeries_poly(PowerSeries):
    cdef ModuleElement __f
