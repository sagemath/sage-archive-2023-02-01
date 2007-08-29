include "../../ext/stdsage.pxi"

cdef extern from *:
     int sprintf_3d "sprintf" (char*, char*, double, double, double)
     int sprintf_3i "sprintf" (char*, char*, int, int, int)
     int sprintf_4i "sprintf" (char*, char*, int, int, int, int)
     int sprintf_9d "sprintf" (char*, char*, double, double, double, double, double, double, double, double, double)

include "../../ext/python_list.pxi"
include "../../ext/python_string.pxi"


from math import sin, cos, sqrt

from sage.rings.real_double import RDF

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from sage.plot.plot3d.base import Graphics3dGroup

from transform cimport Transformation

# --------------------------------------------------------------------
# Fast routines for generating string representations of the polygons.
# --------------------------------------------------------------------

cdef inline format_tachyon_triangle(point_c P, point_c Q, point_c R):
    cdef char ss[250]
    # PyString_FromFormat doesn't do floats?
    cdef Py_ssize_t r = sprintf_9d(ss,
                                   "TRI V0 %g %g %g V1 %g %g %g V2 %g %g %g",
                                   P.x, P.y, P.z,
                                   Q.x, Q.y, Q.z,
                                   R.x, R.y, R.z )
    return PyString_FromStringAndSize(ss, r)

cdef inline format_obj_vertex(point_c P):
    cdef char ss[100]
    # PyString_FromFormat doesn't do floats?
    cdef Py_ssize_t r = sprintf_3d(ss, "v %g %g %g", P.x, P.y, P.z)
    return PyString_FromStringAndSize(ss, r)

cdef inline format_obj_face(face_c face, int off):
    cdef char ss[100]
    cdef Py_ssize_t r, i
    if face.n == 3:
        r = sprintf_3i(ss, "f %d %d %d", face.vertices[0] + off, face.vertices[1] + off, face.vertices[2] + off)
    elif face.n == 4:
        r = sprintf_4i(ss, "f %d %d %d %d", face.vertices[0] + off, face.vertices[1] + off, face.vertices[2] + off, face.vertices[3] + off)
    else:
        return "f " + " ".join([str(face.vertices[i] + off) for i from 0 <= i < face.n])
    # PyString_FromFormat is almost twice as slow
    return PyString_FromStringAndSize(ss, r)

cdef inline format_obj_face_back(face_c face, int off):
    cdef char ss[100]
    cdef Py_ssize_t r, i
    if face.n == 3:
        r = sprintf_3i(ss, "f %d %d %d", face.vertices[2] + off, face.vertices[1] + off, face.vertices[1] + off)
    elif face.n == 4:
        r = sprintf_4i(ss, "f %d %d %d %d", face.vertices[3] + off, face.vertices[2] + off, face.vertices[1] + off, face.vertices[0] + off)
    else:
        return "f " + " ".join([str(face.vertices[i] + off) for i from face.n > i >= 0])
    return PyString_FromStringAndSize(ss, r)



cdef class IndexFaceSet(PrimativeObject):

    def __new__(self, faces, point_list=None, enclosde=False, **kwds):
        self.vs = <point_c *>NULL
        self.face_indices = <int *>NULL
        self._faces = <face_c *>NULL


    def __init__(self, faces, point_list=None, enclosed=False, **kwds):
        PrimativeObject.__init__(self, **kwds)

        self.enclosed = enclosed

        if point_list is None:
            face_list = faces
            faces = []
            point_list = []
            point_index = {}
            for face in face_list:
                iface = []
                for p in face:
                    try:
                        ix = point_index[p]
                    except KeyError:
                        ix = len(point_list)
                        point_index[p] = ix
                        point_list.append(p)
                    iface.append(ix)
                faces.append(iface)

        cdef Py_ssize_t i
        cdef Py_ssize_t index_len = len(faces)
        for i from 0 <= i < len(faces):
            index_len += len(faces[i])

        self.vcount = len(point_list)
        self.fcount = len(faces)

        self.realloc(self.vcount, self.fcount, index_len)

        for i from 0 <= i < self.vcount:
            self.vs[i].x, self.vs[i].y, self.vs[i].z = point_list[i]

        cdef int cur_pt = 0
        for i from 0 <= i < self.fcount:
            self._faces[i].n = len(faces[i])
            self._faces[i].vertices = &self.face_indices[cur_pt]
            for ix in faces[i]:
                self.face_indices[cur_pt] = ix
                cur_pt += 1

    cdef realloc(self, vcount, fcount, icount):
        if self.vs:
            self.vs = <point_c *> sage_realloc(self.vs, sizeof(point_c) * vcount)
        else:
            self.vs = <point_c *> sage_malloc(sizeof(point_c) * vcount)
        if self._faces:
            self._faces = <face_c *> sage_realloc(self._faces, sizeof(face_c) * fcount)
        else:
            self._faces = <face_c *> sage_malloc(sizeof(face_c) * fcount)
        if self.face_indices:
            self.face_indices = <int *> sage_realloc(self.face_indices, sizeof(int) * icount)
        else:
            self.face_indices = <int *> sage_malloc(sizeof(int) * icount)
        if self.vs == NULL or self.face_indices == NULL or self._faces == NULL:
            raise MemoryError, "Out of memory allocating triangulation for %s" % type(self)

    def _clean_point_list(self):
        # TODO: remove redundant/unused points...
        pass

    def __dealloc__(self):
        if self.vs: sage_free(self.vs)
        if self._faces: sage_free(self._faces)
        if self.face_indices: sage_free(self.face_indices)

    def is_enclosed(self):
        return self.enclosed

    def index_faces(self):
        cdef Py_ssize_t i, j
        return [[self._faces[i].vertices[j] for j from 0 <= j < self._faces[i].n] for i from 0 <= i < self.fcount]

    def faces(self):
        points = self.vertices()
        cdef Py_ssize_t i, j
        return [[points[self._faces[i].vertices[j]] for j from 0 <= j < self._faces[i].n] for i from 0 <= i < self.fcount]

    def vertices(self):
        cdef Py_ssize_t i
        return [(self.vs[i].x, self.vs[i].y, self.vs[i].z) for i from 0 <= i < self.vcount]

    def x3d_geometry(self):
        cdef Py_ssize_t i
        points = ",".join(["%s %s %s"%(self.vs[i].x, self.vs[i].y, self.vs[i].z) for i from 0 <= i < self.vcount])
        coordIndex = ",-1,".join([",".join([str(self._faces[i].vertices[j])
                                            for j from 0 <= j < self._faces[i].n])
                                  for i from 0 <= i < self.fcount])
        return """
<IndexedFaceSet coordIndex='%s,-1'>
  <Coordinate point='%s'/>
</IndexedFaceSet>
"""%(coordIndex, points)


    def tachyon_repr(self, render_params):
        cdef Transformation transform = render_params.transform
        lines = []
        cdef point_c P, Q, R
        cdef face_c face
        cdef Py_ssize_t i, k
        for i from 0 <= i < self.fcount:
            face = self._faces[i]
            if transform is not None:
                transform.transform_point_c(&P, self.vs[face.vertices[0]])
                transform.transform_point_c(&Q, self.vs[face.vertices[1]])
                transform.transform_point_c(&R, self.vs[face.vertices[2]])
            else:
                P = self.vs[face.vertices[0]]
                Q = self.vs[face.vertices[1]]
                R = self.vs[face.vertices[2]]
            PyList_Append(lines, format_tachyon_triangle(P, Q, R))
            PyList_Append(lines, self.texture.id)
            if face.n > 3:
                for k from 3 <= k < face.n:
                    Q = R
                    if transform is not None:
                        transform.transform_point_c(&R, self.vs[face.vertices[k]])
                    else:
                        R = self.vs[face.vertices[k]]
                    PyList_Append(lines, format_tachyon_triangle(P, Q, R))
                    PyList_Append(lines, self.texture.id)

        return lines


    def obj_repr(self, render_params):

        cdef Transformation transform = render_params.transform
        cdef int off = render_params.obj_vertex_offset
        cdef Py_ssize_t i
        cdef point_c res

        if transform is None:
#            points = ["v %s %s %s"%(self.vs[i].x, self.vs[i].y, self.vs[i].z) for i from 0 <= i < self.vcount]
            points = [format_obj_vertex(self.vs[i]) for i from 0 <= i < self.vcount]
        else:
            points = []
            for i from 0 <= i < self.vcount:
                transform.transform_point_c(&res, self.vs[i])
                PyList_Append(points, format_obj_vertex(res))

        faces = [format_obj_face(self._faces[i], off) for i from 0 <= i < self.fcount]
        if not self.enclosed:
            back_faces = [format_obj_face_back(self._faces[i], off) for i from 0 <= i < self.fcount]
        else:
            back_faces = []

        render_params.obj_vertex_offset += self.vcount

        return ["g " + render_params.unique_name('obj'),
                "usemtl " + self.texture.id,
                points,
                faces,
                back_faces]

    def stickers(self, colors, width, hover):
        all = []
        n = self.fcount; ct = len(colors)
        for k in range(len(colors)):
            if colors[k]:
                all.append(self.sticker(range(k,n,ct), width, hover, texture=colors[k]))
        return Graphics3dGroup(all)

    def sticker(self, face_list, width, hover, **kwds):
        if not isinstance(face_list, (list, tuple)):
            face_list = (face_list,)
        faces = self.getFaceList()
        all = []
        for k in face_list:
            all.append(sticker(faces[k], width, hover))
        return IndexFaceSet(all, **kwds)


cdef class FaceIter:
    def __init__(self, face_set):
        self.set = face_set
        self.i = 0
    def __iter__(self):
        return self
    def __next__(self):
        if self.i >= self.set.fcount:
            raise StopIteration
        else:
            return [self.set._faces[self.i].vertices[j] for j from 0 <= j < self.set._faces[self.i].n]
            self.i += 1

cdef class VertexIter:
    def __init__(self, face_set):
        self.set = face_set
        self.i = 0
    def __iter__(self):
        return self
    def __next__(self):
        if self.i >= self.set.vcount:
            raise StopIteration
        else:
            return (self.set.vs[self.i].x, self.set.vs[self.i].y, self.set.vs[self.i].z)
            self.i += 1

def len3d(v):
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def sticker(face, width, hover):
    n = len(face)
    edges = []
    for i from 0 <= i < n:
        edges.append(vector(RDF, [face[i-1][0] - face[i][0],
                                  face[i-1][1] - face[i][1],
                                  face[i-1][2] - face[i][2]]))
    sticker = []
    for i in range(n):
        v = -edges[i]
        w = edges[i-1]
        N = v.cross_product(w)
        lenN = len3d(N)
        dv = v*(width*len3d(w)/lenN)
        dw = w*(width*len3d(v)/lenN)
        sticker.append(tuple(vector(RDF, face[i-1]) + dv + dw + N*(hover/lenN)))
    return sticker


class RectangularGridSurface(IndexFaceSet):
    def __init__(self, **kwds):
        PrimativeObject.__init__(self, **kwds)
    def getFaceList(self):
        #return [[(0,0,0), (0,1,0), (0,1,1), (0,0,1)]]
        grid = self.getGrid()
        faces = []
        for i in range(1, len(grid)):
            line = grid[i]
            last_line = grid[i-1]
            for j in range(1, len(line)):
                face = [line[j], line[j-1], last_line[j-1], last_line[j]]
                # remove extra vertex of degenerate quads.
                if   face[3] == face[0]: face.remove(face[0])
                elif face[0] == face[1]: face.remove(face[1])
                elif face[1] == face[2]: face.remove(face[2])
                elif face[2] == face[3]: face.remove(face[3])
                faces.append(face)
        return faces
