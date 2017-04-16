from sage.structure.element_wrapper cimport ElementWrapper

cdef class LieAlgebraElementWrapper(ElementWrapper):
    cpdef _add_(self, right)
    cpdef _sub_(self, right)

cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    pass

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    cpdef bracket(self, right)
    cpdef _bracket_(self, right)
    cpdef to_vector(self)
    cpdef dict monomial_coefficients(self, bint copy=*)
    #cpdef lift(self)

