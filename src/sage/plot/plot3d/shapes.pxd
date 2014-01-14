from base cimport PrimitiveObject
from index_face_set cimport IndexFaceSet, point_c, face_c
from parametric_surface cimport ParametricSurface


cdef class Cone(ParametricSurface):
    cdef double radius
    cdef double height
    cdef bint closed

cdef class Cylinder(ParametricSurface):
    cdef double radius
    cdef double height
    cdef bint closed

cdef class Sphere(ParametricSurface):
    cdef double radius

cdef class Torus(ParametricSurface):
    cdef double R, r
