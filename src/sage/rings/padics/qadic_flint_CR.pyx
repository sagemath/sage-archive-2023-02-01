include "sage/libs/linkages/padics/fmpz_poly_unram.pxi"
include "CR_template.pxi"

cdef class PowComputer_(PowComputer_flint_unram):
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        self._prec_type = 'capped-rel'
        PowComputer_flint_unram.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

cdef class qAdicCappedRelativeElement(CRElement):
    pass
