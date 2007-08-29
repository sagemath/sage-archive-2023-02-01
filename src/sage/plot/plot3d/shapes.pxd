from base cimport PrimativeObject
from index_face_set cimport IndexFaceSet, point_c, face_c
from parametric_surface cimport ParametricSurface


cdef class Cone(ParametricSurface):
    cdef double radius
    cdef double height

cdef class Cylinder(ParametricSurface):
    cdef double radius
    cdef double height

cdef class Sphere(ParametricSurface):
    cdef double radius

cdef class Torus(ParametricSurface):
    cdef double R, r