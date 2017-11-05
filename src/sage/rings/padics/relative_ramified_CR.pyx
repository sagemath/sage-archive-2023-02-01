include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "CR_template.pxi"

cdef class RelativeRamifiedCappedRelativeElement(CRElement):
    def _poly_rep(self):
        return self.unit
