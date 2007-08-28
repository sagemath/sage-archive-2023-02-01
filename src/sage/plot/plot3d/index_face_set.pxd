from sage.plot.plot3d.base cimport PrimativeObject

cdef struct point_c:
  double x
  double y
  double z

cdef struct face_c:
  int n
  int* vertices


cdef class IndexFaceSet(PrimativeObject):
    cdef bint enclosed
    cdef Py_ssize_t fcount, vcount, icount
    cdef alloc(self, vcount, fcount, icount)
    # array of {x,y,z}
    cdef point_c* vs
    # array of array of fcount indices into points, each ending with -1
    cdef int* face_indices
    # pointers into face_indices marking the begining of each face
    cdef face_c* _faces

cdef class FaceIter:
    cdef Py_ssize_t i
    cdef IndexFaceSet set


cdef class VertexIter:
    cdef Py_ssize_t i
    cdef IndexFaceSet set
