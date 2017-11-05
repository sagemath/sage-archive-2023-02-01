include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "CA_template.pxi"

cdef class RelativeRamifiedCappedAbsoluteElement(CAElement):
    def _poly_rep(self):
        return self.value
