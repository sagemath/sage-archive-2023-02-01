include "../../libs/polybori/decl.pxi"

from sage.structure.parent_base cimport ParentWithBase
from sage.rings.polynomial.multi_polynomial_ring_generic cimport \
                                                MPolynomialRing_generic
from sage.rings.polynomial.multi_polynomial cimport MPolynomial

cdef class BPRing(MPolynomialRing_generic):
    cdef PBRing _R

cdef class BPolynomial(MPolynomial):
    cdef PBPoly *_P
