include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "FM_template.pxi"

cdef class RelativeRamifiedFixedModElement(FMElement):
    def _poly_rep(self):
        return self.value
