from sage.structure.element cimport MultiplicativeGroupElement, Element, MonoidElement, Matrix
from sage.groups.libgap_wrapper cimport ElementLibGAP

cpdef is_MatrixGroupElement(x)

cdef class MatrixGroupElement_generic(MultiplicativeGroupElement):
    cdef public Matrix _matrix

    cpdef _act_on_(self, x, bint self_on_left)
    cpdef int _cmp_(self, Element other) except -2
    cpdef list list(self)
    cpdef MonoidElement _mul_(self, MonoidElement other)
    cpdef is_one(self)

cdef class MatrixGroupElement_gap(ElementLibGAP):
    cpdef _act_on_(self, x, bint self_on_left)
    cpdef int _cmp_(self, Element other) except -2
    cpdef list list(self)

