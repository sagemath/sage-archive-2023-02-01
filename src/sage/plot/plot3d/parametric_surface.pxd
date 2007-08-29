from index_face_set cimport IndexFaceSet, point_c, face_c

cdef class ParametricSurface(IndexFaceSet):
    cdef object f
    cdef object render_grid
    cdef eval_c(self, point_c *res, double u, double v)

