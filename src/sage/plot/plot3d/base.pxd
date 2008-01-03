from sage.structure.sage_object cimport SageObject


cdef class Graphics3d(SageObject):
    cdef object texture

cdef class PrimitiveObject(Graphics3d):
    pass
