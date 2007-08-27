from math import sin, cos, sqrt

from sage.rings.real_double import RDF

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from sage.plot.plot3d.base import Graphics3dGroup, PrimativeObject
from sage.plot.plot3d.base cimport PrimativeObject


cdef class IndexFaceSet(PrimativeObject):
    def __init__(self, faces, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.faces = faces

    def getFaceList(self):
        return self.faces

    def x3d_geometry(self):
        faces = self.getFaceList()
        point_index = {}
        point_list = []
        coordIndex = ""
        for face in faces:
            for p in face:
                try:
                    ix = point_index[p]
                except KeyError:
                    ix = len(point_list)
                    point_index[p] = ix
                    point_list.append(p)
                coordIndex += ",%s"%ix
            coordIndex += ",-1"
        points = ",".join(["%s %s %s"%p for p in point_list])
        return """
<IndexedFaceSet coordIndex='%s'>
  <Coordinate point='%s'/>
</IndexedFaceSet>
"""%(coordIndex, points)

    def tachyon_str(self, transform):
        faces = self.getFaceList()
        s = ""
        for face in faces:
            if transform is not None:
                face = [transform(p) for p in face]
            s += """
TRI
    V0 %s %s %s
    V1 %s %s %s
    V2 %s %s %s
                """ % (face[0][0], face[0][1], face[0][2],
                       face[1][0], face[1][1], face[1][2],
                       face[2][0], face[2][1], face[2][2])
            s += "\n" + self.texture.tachyon_str()
            if len(face) > 3:
                for k in range(2, len(face)-1):
                    s += """
TRI
    V0 %s %s %s
    V1 %s %s %s
    V2 %s %s %s
                """ % (face[0][0], face[0][1], face[0][2],
                       face[k][0], face[k][1], face[k][2],
                       face[k+1][0], face[k+1][1], face[k+1][2])
                    s += "\n" + self.texture.tachyon_str()
        return s


    def obj_geometry(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        faces = self.getFaceList()
        point_index = {}
        face_list = []
        start_index = len(point_list)
        for face in faces:
            cur_face = "f"
            for p in face:
                if transform is not None:
                    p = transform(p)
                try:
                    ix = point_index[p]
                except KeyError:
                    ix = len(point_list)
                    point_index[p] = ix
                    point_list.append(p)
                cur_face += " %s"%(ix+1)
            face_list.append(cur_face)
        s = "\n".join(["v %s %s %s"%(p[0], p[1], p[2]) for p in point_list[start_index:]])
        s += "\n"
        s += "\n".join(face_list)
        s += "\n\n"
        return s

    def stickers(self, colors, width, hover):
        faces = self.getFaceList()
        all = []
        n = len(faces); ct = len(colors)
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

def len3d(v):
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def sticker(face, width, hover):
    n = len(face)
    edges = []
    for i in range(n):
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
