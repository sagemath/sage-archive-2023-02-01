from types import MethodType

include "sage/libs/linkages/padics/fmpz_poly_unram.pxi"
include "sage/libs/linkages/padics/unram_shared.pxi"
include "FM_template.pxi"

cdef class PowComputer_(PowComputer_flint_unram):
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        self._prec_type = 'fixed-mod'
        PowComputer_flint_unram.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

cdef class qAdicFixedModElement(FMElement):
    frobenius = MethodType(frobenius_unram, None, qAdicFixedModElement)
    trace = MethodType(trace_unram, None, qAdicFixedModElement)
    norm = MethodType(norm_unram, None, qAdicFixedModElement)
