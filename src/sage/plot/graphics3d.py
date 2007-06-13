r"""
3D Graphics objects and plotting, esp. for integration with x3d

AUTHOR:
    -- Robert Bradshaw

TODO:
   -- integrate tachyon
"""
from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.real_double import RDF
from sage.misc.functional import sqrt, atan
from sage.functions.all import *

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
        composite_transform = self.get_transformation() * transform
        return "\n".join([g.tachyon_str(composite_transform) for g in self.all])

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


class PrimativeObject(Graphics3d_new):
    def __init__(self, **kwds):
        try:
            self.color = kwds['color']
        except KeyError:
            self.color = [1,0,0]

    def x3d_str(self):
        return "<Shape>" + self.x3d_geometry() + self.x3d_appearance() + "</Shape>\n"

    def x3d_appearance(self):
        return "<Appearance><Material diffuseColor='%s %s %s' shininess='0.1' specularColor='0.7 0.7 0.7'/></Appearance>"%(self.color[0], self.color[1], self.color[2])

    def tachyon_str(self, transform):
        return self.tachyon_geometry(transform) + self.tachyon_appearance()

    def tachyon_appearance(self):
        return """
   Texture Ambient 0.2 Diffuse 0.8 Specular 0.0 Opacity 1.0
      Color 1.0 0.0 0.0
      TexFunc 0"""

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

class Cone(PrimativeObject):
    def __init__(self, radius, height, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.radius = radius
        self.height = height
    def x3d_geometry(self):
        return "<Cone bottomRadius='%s' height='%s'/>"%(self.radius, self.height)

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

class Text(PrimativeObject):
    def __init__(self, string, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.string = string
    def x3d_geometry(self):
        return "<Text string='%s' solid='true'/>"%self.string

class IndexFaceSet(PrimativeObject):
    def __init__(self, faces, **kwds):
        raise NotImplementedException, "No generic IndexFaceSet"
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
                faces.append([line[j], line[j-1], last_line[j-1], last_line[j]])
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
        urange = [twoPi*k/u_divs for k in range(u_divs+1)]
        vrange = [twoPi*k/v_divs for k in range(v_divs+1)]
        ParametricSurface.__init__(self, None, urange, vrange, **kwds)
        self.R = RDF(R)
        self.r = RDF(r)
    def f(self,u,v):
        return (self.R+self.r*sin(v))*sin(u), (self.R+self.r*sin(v))*cos(u), self.r*cos(v)
