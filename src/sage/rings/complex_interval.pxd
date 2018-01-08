from sage.libs.mpfr.types cimport mpfr_prec_t
from sage.libs.mpfi.types cimport mpfi_t

cimport sage.structure.element
from .real_mpfi cimport RealIntervalFieldElement


cdef class ComplexIntervalFieldElement(sage.structure.element.FieldElement):
    cdef mpfi_t __re
    cdef mpfi_t __im
    cdef mpfr_prec_t _prec

    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cdef RealIntervalFieldElement abs_c(ComplexIntervalFieldElement self)
    cdef RealIntervalFieldElement norm_c(ComplexIntervalFieldElement self)

    cdef inline ComplexIntervalFieldElement _new(self):
        """
        Quickly create a new complex interval with the same parent as
        ``self``.
        """
        cdef type t = type(self)
        return t.__new__(t, self._parent)
