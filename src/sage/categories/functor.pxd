from sage.structure.sage_object cimport SageObject

cdef class Functor(SageObject):
    cdef object __domain
    cdef object __codomain