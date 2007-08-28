from sage.structure.sage_object cimport SageObject


cdef class Graphics3d_new(SageObject):
    cdef object texture

cdef class PrimativeObject(Graphics3d_new):
    pass