from sage.structure.element_wrapper cimport ElementWrapper
from sage.structure.sage_object cimport SageObject

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

cdef class LieObject(SageObject):
    cdef tuple _word
    cpdef tuple to_word(self)

cdef class LieGenerator(LieObject):
    cdef public str _name

cdef class LieBracket(LieObject):
    cdef public LieObject _left
    cdef public LieObject _right
    cdef long _hash

    cpdef lift(self, dict UEA_gens_dict)

cdef class GradedLieBracket(LieBracket):
    cdef public _grade

cdef class LyndonBracket(GradedLieBracket):
    pass
