r"""
3D Graphics objects and plotting.

AUTHOR:
    -- Robert Bradshaw

TODO:
   -- integrate tachyon
"""
from random import randint

from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.real_double import RDF
from sage.misc.functional import sqrt, atan
from sage.functions.all import *
from texture import *

default_texture = Texture()

class Graphics3d_new(SageObject):

    def __add__(self, other):
        return Graphics3dGroup([self, other])

    def translate(self, *x):
        if isinstance(x[0], (tuple, list)):
            x = x[0]
        return TransformGroup([self], trans=x)


    def scale(self, *x):
        if isinstance(x[0], (tuple, list)):
            x = x[0]
        return TransformGroup([self], scale=x)

    def rotate(self, v, theta):
        vx, vy, vz = v
        return TransformGroup([self], rot=[vx, vy, vz, theta])

    def rotateX(self, theta):
        return self.rotate((1,0,0), theta)

    def rotateY(self, theta):
        return self.rotate((0,1,0), theta)

    def rotateZ(self, theta):
        return self.rotate((0,0,1), theta)


    def x3d(self):
        return """
<X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
<head>
<meta name='title' content='sage3d'/>
</head>
<Scene>
%s
%s
</Scene>
</X3D>
"""%(self.viewpoint().x3d_str(), self.x3d_str())

    def viewpoint(self):
        return Viewpoint(0,0,6)

    def show(self, **kwds):
        f = open("test.x3d", "w")
        f.write(self.x3d())
        f.close()

    def mtl_str(self):
        return "\n\n".join([t.mtl_str() for t in self.mtl_set()]) + "\n"

    def mtl_set():
        return set()

class Viewpoint(Graphics3d_new):

    def __init__(self, *x):
        if isinstance(x[0], (tuple, list)):
            x = tuple(x[0])
        self.pos = x

    def x3d_str(self):
        return "<Viewpoint position='%s %s %s'/>"%self.pos

class Graphics3dGroup(Graphics3d_new):
    def __init__(self, all=[]):
        self.all = all

    def tachyon_str(self, transform):
        return "\n".join([g.tachyon_str(transform) for g in self.all])

    def x3d_str(self):
        return "\n".join([g.x3d_str() for g in self.all])

    def obj_str(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        return "\n\n".join([g.obj_str(transform, point_list) for g in self.all]) + "\n"

    def mtl_set(self):
        return reduce(set.union, [g.mtl_set() for g in self.all])

class TransformGroup(Graphics3dGroup):

    def __init__(self, all=[], rot=None, trans=None, scale=None):
        Graphics3dGroup.__init__(self, all)
        self._rot = rot
        self._trans = trans
        if scale is not None and len(scale) == 1:
            if isinstance(scale, (tuple, list)):
                scale = scale[0]
            scale = (scale, scale, scale)
        self._scale = scale

    def x3d_str(self):
        s = "<Transform"
        if self._rot is not None:
            s += " rotation='%s %s %s %s'"%tuple(self._rot)
        if self._trans is not None:
            s += " translation='%s %s %s'"%tuple(self._trans)
        if self._scale is not None:
            s += " scale='%s %s %s'"%tuple(self._scale)
        s += ">\n"
        s += Graphics3dGroup.x3d_str(self)
        s += "\n</Transform>"
        return s

    def tachyon_str(self, transform):
        if transform is None:
            composite_transform = self.get_transformation()
        else:
            composite_transform = self.get_transformation() * transform
        return "\n".join([g.tachyon_str(composite_transform) for g in self.all])

    def obj_str(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        if transform is None:
            composite_transform = self.get_transformation()
        else:
            composite_transform = self.get_transformation() * transform
        return "\n\n".join([g.obj_str(composite_transform, point_list) for g in self.all])

    def get_transformation(self):
        try:
            return self.T
        except AttributeError:
            self.T = Transformation(self._scale, self._rot, self._trans)
            return self.T

class Transformation:
    def __init__(self, scale=(1,1,1),
                       rot=None,
                       trans=[0,0,0]):

        if scale is None:
            scale = (1,1,1)
        if trans is None:
            trans = [0,0,0]

        # TODO: determine for sure if x3d does scale or rotation first
        m = matrix(RDF, 3, 3,
                  [scale[0], 0, 0, 0, scale[1], 0, 0, 0, scale[2]])

        if rot is not None:
            # rotate about v by theta
            vx, vy, vz, theta = rot
            if vz != 0:
                t = atan(vy/vz) + pi/2
                m = self.rotX(t) * m
                new_y = vy*cos(t) - vz*sin(t)
            else:
                t = 0
                new_y = vy
            # v = [vx, new_y, 0]
            if new_y != 0:
                s = atan(vx/new_y) + pi/2
                m = self.rotZ(s) * m
            else:
                s = 0
            # v = [new_x, 0, 0]
            m = self.rotX(theta) * m
            # now put back to our former reference frame
            m = self.rotZ(-s) * m
            m = self.rotX(-t) * m

        self.matrix = m.augment(matrix(RDF, 3, 1, list(trans))) \
                       .stack(matrix(RDF, 1, 4, [0,0,0,1]))

    def rotX(self, theta):
        return matrix(RDF, 3, 3, [1, 0, 0,
                                  0, cos(theta), -sin(theta),
                                  0, sin(theta), cos(theta)])

    def rotZ(self, theta):
        return matrix(RDF, 3, 3, [cos(theta), -sin(theta), 0,
                                  sin(theta), cos(theta), 0,
                                  0, 0, 1])

    def transform_point(self, x):
        Tx = self.matrix * vector(RDF, [x[0], x[1], x[2], 1])
        return (Tx[0], Tx[1], Tx[2])

    def transform_vector(self, x):
        Tx = self.matrix * vector(RDF, [x[0], x[1], x[2], 0])
        return (Tx[0], Tx[1], Tx[2])

    def __mul__(self, other):
        T = Transformation()
        T.matrix = self.matrix * other.matrix
        return T

    def __call__(self, p):
        res = self.matrix * vector(RDF, [p[0], p[1], p[2], 1])
        return tuple(res[:3])

class PrimativeObject(Graphics3d_new):
    def __init__(self, **kwds):
        try:
            self.texture = kwds['texture']
            if not isinstance(self.texture, Texture_class):
                self.texture = Texture(self.texture)
        except KeyError:
            self.texture = default_texture

    def x3d_str(self):
        return "<Shape>" + self.x3d_geometry() + self.texture.x3d_str() + "</Shape>\n"

    def tachyon_str(self, transform):
        return self.tachyon_geometry(transform) + self.texture.tachyon_str()

    def tachyon_geometry(self, transform):
        return self.triangulation().tachyon_geometry(transform)

    def obj_str(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        return "g obj%s\n\nusemtl "%randint(0,10000000) + self.texture.id + "\n" + self.obj_geometry(transform, point_list)

    def obj_geometry(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        return self.triangulation().obj_geometry(transform, point_list)

    def mtl_set(self):
        return set([self.texture])

class Light(PrimativeObject):
    def __init__(self, intensity=.3, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.intensity = intensity
    def x3d_geometry(self):
        return "<Light intensity='%s'/>"%("%s %s %s"%self.location, self.intensity)

class Box(PrimativeObject):
    def __init__(self, *size, **kwds):
        PrimativeObject.__init__(self, **kwds)
        if isinstance(size[0], (tuple, list)):
            size = size[0]
        self.size = size
    def x3d_geometry(self):
        return "<Box size='%s %s %s'/>"%self.size
    def triangulation(self):
        """
        Returns an IndexFaceSet (which may be either triangles or quadrilaterals).
        """
        x, y, z = self.size
        faces = [[(x, y, z), (-x, y, z), (-x,-y, z), ( x,-y, z)],
                 [(x, y, z), ( x, y,-z), (-x, y,-z), (-x, y, z)],
                 [(x, y, z), ( x,-y, z), ( x,-y,-z), ( x, y,-z)] ]
        faces += [reversed([(-x,-y,-z) for x,y,z in face]) for face in faces]
        return IndexFaceSet(faces)

def ColorCube(size, colors):
    if not isinstance(size, (tuple, list)):
        size = (size, size, size)
    box = Box(size)
    faces = box.triangulation().getFaceList()
    if len(colors) == 3:
        colors = colors * 2
    all = []
    for k in range(6):
        all.append(IndexFaceSet([faces[k]], texture=colors[k]))
    return Graphics3dGroup(all)

class Cone(PrimativeObject):
    def __init__(self, radius, height, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.radius = radius
        self.height = height
    def x3d_geometry(self):
        return "<Cone bottomRadius='%s' height='%s'/>"%(self.radius, self.height)
    def triangulation(self, res=30):
        def f(u, v):
            if u == -1:
                return (0,0,0)
            elif u == 1:
                return (0,0,self.height)
            else:
                return (self.radius*sin(v), self.radius*cos(v), 0)
        twoPi = RDF(2*pi)
        return ParametricSurface(f, [-1,0,1], [twoPi*k/res for k in range(res)] + [RDF(0)])

class Cylinder(PrimativeObject):
    def __init__(self, radius, height, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.radius = radius
        self.height = height
    def x3d_geometry(self):
        return "<Cylinder radius='%s' height='%s'/>"%(self.radius, self.height)

    def tachyon_geometry(self, transform):
        cen = transform.transform_point((0,0,0))
        axis = transform.transform_vector((0,self.height,0))
        radv = transform.transform_vector((self.radius,0,0)) # Tachyon can't do sqashed, just scale uniform
        rad = sqrt(sum([x*x for x in radv]))
        return """
FCylinder
   Center %s %s %s
   Axis %s %s %s
   Rad %s
        """%(cen[0], cen[1], cen[2], axis[0], axis[1], axis[2], rad)

    def triangulation(self, res=30):
        def f(u, v):
            if u == -2:
                return (0, 0, -self.height)
            elif u == -1:
                return (self.radius*sin(v), self.radius*cos(v), -self.height)
            elif u == 1:
                return (self.radius*sin(v), self.radius*cos(v), self.height)
            else: # u == 2:
                return (0, 0, self.height)
        twoPi = RDF(2*pi)
        return ParametricSurface(f, [-2,-1,1,2], [twoPi*k/res for k in range(res)] + [RDF(0)])


class Sphere(PrimativeObject):
    def __init__(self, radius, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.radius = radius

    def x3d_geometry(self):
        return "<Sphere radius='%s'/>"%(self.radius)

    def tachyon_geometry(self, transform):
        cen = transform.transform_point((0,0,0))
        radv = transform.transform_vector((1,0,0)) # Tachyon can't do ellipsoids, just scale uniform
        rad = sqrt(sum([x*x for x in radv]))
        return """
Sphere center %s %s %s
   Rad %s
        """%(cen[0], cen[1], cen[2], rad)

    def triangulation(self, res=30, vres=None):
        if vres is None:
            vres = res
        def f(u, v):
            if u == -10:
                return (0, 0, -self.radius)
            elif u == 10:
                return (0, 0, self.radius)
            else:
                return (self.radius*sin(v) * cos(u), self.radius*cos(v) * cos(u), self.radius * sin(u))
        twoPi = RDF(2*pi)
        return ParametricSurface(f, [-10] + [twoPi*k/vres - twoPi/4 for k in range(1,vres)] + [10], [twoPi*k/res for k in range(res)] + [RDF(0)])

class Text(PrimativeObject):
    def __init__(self, string, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.string = string
    def x3d_geometry(self):
        return "<Text string='%s' solid='true'/>"%self.string

class IndexFaceSet(PrimativeObject):
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
                try:
                    ix = point_index[p]
                except KeyError:
                    ix = len(point_list)
                    point_index[p] = ix
                    point_list.append(p)
                cur_face += " %s"%(ix+1)
            face_list.append(cur_face)
        if transform is not None:
            point_list = [transform(p) for p in point_list]
        s = "\n".join(["v %s %s %s"%(p[0], p[1], p[2]) for p in point_list[start_index:]])
        s += "\n"
        s += "\n".join(face_list)
        s += "\n\n"
        return s

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
                if   face[3] == face[0]: face.remove(face[0])
                elif face[0] == face[1]: face.remove(face[1])
                elif face[1] == face[2]: face.remove(face[2])
                elif face[2] == face[3]: face.remove(face[3])
                faces.append(face)
        return faces


# TODO: I want to be able to override f in subclasses, but also pass in an f
class ParametricSurface(RectangularGridSurface):
    def __init__(self, f, urange, vrange, **kwds):
        PrimativeObject.__init__(self, **kwds)
        if f is not None:
            self.f = f
        self.urange = urange
        self.vrange = vrange
    def getGrid(self):
        return [[self.f(u,v) for u in self.urange] for v in self.vrange]

class Torus(ParametricSurface):
# e.g  show(sum([Torus(1,.03,20,20, color=[1, float(t/30), 0]).rotate((1,1,1),t) for t in range(30)], Sphere(.3)))
    def __init__(self, R=1, r=.3, u_divs=10, v_divs=10, **kwds):
        twoPi = RDF(2*pi)
        urange = [twoPi*k/u_divs for k in range(u_divs)] + [RDF(0)]
        vrange = [twoPi*k/v_divs for k in range(v_divs)] + [RDF(0)]
        ParametricSurface.__init__(self, None, urange, vrange, **kwds)
        self.R = RDF(R)
        self.r = RDF(r)
    def f(self,u,v):
        return (self.R+self.r*sin(v))*sin(u), (self.R+self.r*sin(v))*cos(u), self.r*cos(v)
