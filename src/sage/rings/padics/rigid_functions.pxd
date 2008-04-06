
from sage.structure.sage_object cimport SageObject
from sage.rings.real_mpfi cimport RealIntervalFieldElement
from sage.rings.integer cimport Integer

cdef class RigidAnalyticFunction(SageObject):
    pass

cdef class RigidAnalyticFunction_disc(SageObject):
    cdef _base_ring
    cdef _log_of_radius
    cdef Integer _e
    cdef Integer _p
    cdef _RIF_p
    cdef _SR_p
    cdef RealIntervalFieldElement _meL
    cdef _log_order
    cdef RealIntervalFieldElement _en
    cdef _log_offset
    cdef RealIntervalFieldElement _eC
    cdef _coeffs

    cpdef Integer coeff_val_bound(self, i)
    cpdef RealIntervalFieldElement coeff_val_bound_RIF(self, _i)
    cpdef Integer term_val_bound(self, i, val)
    cpdef RealIntervalFieldElement term_val_bound_RIF(self, _i, val)
    cpdef output_valuation(self, input_valuation)
    cpdef output_valuation_fast(self, input_valuation)
    cpdef valuation_low_point(self, input_valuation)
    cpdef valuation_low_point_pair(self, input_valuation)
    cpdef RealIntervalFieldElement valuation_low_point_fast(self, input_valuation)

cdef class RigidAnalyticFunctionCoeffs(SageObject):
    cdef _base_ring

cdef class ExpCoeffs(RigidAnalyticFunctionCoeffs):
    pass
cdef class LogCoeffs(RigidAnalyticFunctionCoeffs):
    pass
cdef class RigidAnalyticFunction_exp(RigidAnalyticFunction_disc):
    pass
cdef class RigidAnalyticFunction_log(RigidAnalyticFunction_disc):
    pass
