include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "FP_template.pxi"

cdef class RelativeRamifiedFloatingPointElement(FPElement):
    def _poly_rep(self):
        return self.unit
