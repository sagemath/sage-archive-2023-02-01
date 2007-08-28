include "../../ext/stdsage.pxi"


from math import sin, cos, sqrt

from sage.rings.real_double import RDF

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from sage.plot.plot3d.base import Graphics3dGroup



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

    def tachyon_str(self, transform):
        faces = self.faces()
        lines = []
        for face in faces:
            if transform is not None:
                face = [transform(p) for p in face]
            lines.append("""
TRI
    V0 %s %s %s
    V1 %s %s %s
    V2 %s %s %s
                """ % (face[0][0], face[0][1], face[0][2],
                       face[1][0], face[1][1], face[1][2],
                       face[2][0], face[2][1], face[2][2]))
            lines.append(self.texture.id)
            if len(face) > 3:
                for k in range(2, len(face)-1):
                    lines.append("""
TRI
    V0 %s %s %s
    V1 %s %s %s
    V2 %s %s %s
                """ % (face[0][0], face[0][1], face[0][2],
                       face[k][0], face[k][1], face[k][2],
                       face[k+1][0], face[k+1][1], face[k+1][2]))
                    lines.append(self.texture.id)
        return "".join(lines)


    def obj_geometry(self, transform, point_offset=0):
        if transform is None:
            points = ["v %s %s %s"%(self.vs[i].x, self.vs[i].y, self.vs[i].z) for i from 0 <= i < self.vcount]
        else:
            points = ["v %s %s %s"%transform(self.vs[i].x, self.vs[i].y, self.vs[i].z) for i from 0 <= i < self.vcount]
        faces = [" ".join([str(self._faces[i].vertices[j]+1) for j from 0 <= j < self._faces[i].n]) for i from 0 <= i < self.fcount]
        if not self.enclosed:
            faces += [" ".join([str(self._faces[i].vertices[j]+1) for j from self._faces[i].n > j >= 0]) for i from 0 <= i < self.fcount]
        return "%s\n\nf %s\n\n" % ("\n".join(points), "\nf ".join(faces))

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
