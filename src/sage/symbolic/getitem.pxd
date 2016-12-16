from sage.symbolic.expression cimport Expression
from sage.structure.sage_object cimport SageObject

cdef class OperandsWrapper(SageObject):
    cdef Expression _expr
