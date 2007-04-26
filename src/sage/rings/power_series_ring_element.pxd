from sage.structure.element cimport AlgebraElement

cdef class PowerSeries(AlgebraElement):
    cdef char __is_gen
    cdef _prec
    cdef common_prec_c(self, PowerSeries other)

cdef class PowerSeries_poly(PowerSeries):
    cdef __f
