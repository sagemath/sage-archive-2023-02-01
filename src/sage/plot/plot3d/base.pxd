from sage.structure.sage_object cimport SageObject


cdef class Graphics3d(SageObject):
    cdef object texture

cdef class PrimativeObject(Graphics3d):
    pass