from sage.structure.sage_object cimport SageObject


cdef class Graphics3d(SageObject):
    cdef object texture
    cdef object _aspect_ratio
    cdef object _frame_aspect_ratio

cdef class PrimitiveObject(Graphics3d):
    pass
