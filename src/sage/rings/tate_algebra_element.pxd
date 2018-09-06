from sage.structure.element cimport CommutativeAlgebraElement
#from sage.rings.integer cimport Integer

cdef class TateAlgebraElement(CommutativeAlgebraElement):
    cdef _poly
    cdef _prec
    cdef list _coeffs
