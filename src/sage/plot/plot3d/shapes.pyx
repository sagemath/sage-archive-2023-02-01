"""
Basic objects such as Sphere, Box, Cone, etc.

AUTHORS:
    -- Robert Bradshaw 2007-02: inital version
    -- Robert Bradshaw 2007-08: obj/tachon rendering, much updating
    -- Robert Bradshaw 2007-08: cythonization

EXAMPLES:
    sage: from sage.plot.plot3d.shapes import *
    sage: S = Sphere(.5, color='yellow')
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)
    sage: S += Sphere(.1, color='white').translate(.45,-.1,.15) + Sphere(.05, color='black').translate(.51,-.1,.17)
    sage: S += Sphere(.1, color='white').translate(.45, .1,.15) + Sphere(.05, color='black').translate(.51, .1,.17)
    sage: S += Sphere(.1, color='yellow').translate(.5, 0, -.2)
    sage: S.show()
    sage: S.scale(1,1,2).show()

    sage: from sage.plot.plot3d.shapes import *
    sage: Torus(.7, .2, color=(0,.3,0)).show()
"""


#*****************************************************************************
#      Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef extern from "math.h":
    double sqrt(double)
    double sin(double)
    double cos(double)
    double tan(double)
    double asin(double)
    double acos(double)
    double atan(double)

from sage.rings.real_double  import RDF
from sage.modules.free_module_element import vector

from base import Graphics3dGroup


class Box(IndexFaceSet):

    def __init__(self, *size, **kwds):
        if isinstance(size[0], (tuple, list)):
            size = size[0]
        self.size = size
        x, y, z = self.size
        faces = [[(x, y, z), (-x, y, z), (-x,-y, z), ( x,-y, z)],
                 [(x, y, z), ( x, y,-z), (-x, y,-z), (-x, y, z)],
                 [(x, y, z), ( x,-y, z), ( x,-y,-z), ( x, y,-z)] ]
        faces += [list(reversed([(-x,-y,-z) for x,y,z in face])) for face in faces]
        IndexFaceSet.__init__(self, faces, enclosed=True, **kwds)

    def x3d_geometry(self):
        return "<Box size='%s %s %s'/>"%self.size


def ColorCube(size, colors):
    if not isinstance(size, (tuple, list)):
        size = (size, size, size)
    box = Box(size)
    faces = box.face_list()
    if len(colors) == 3:
        colors = colors * 2
    all = []
    for k in range(6):
        all.append(IndexFaceSet([faces[k]], enclosed=True, texture=colors[k]))
    return Graphics3dGroup(all)

cdef class Cone(ParametricSurface):

    def __init__(self, radius, height, closed=True, **kwds):
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius
        self.height = height
        self.closed = closed

    def x3d_geometry(self):
        return "<Cone bottomRadius='%s' height='%s'/>"%(self.radius, self.height)

    def get_grid(self, ds):
        twoPi = 2*RDF.pi()
        v_res = min(max(int(twoPi*self.radius/ds), 5), 37)
        if self.closed:
            urange = [1,0,-1]
        else:
            urange = [1,0]
        vrange = [float(twoPi*k/v_res) for k in range(v_res)] + [0.0]
        return urange, vrange

    cdef eval_c(self, point_c *res, double u, double v):
        if u == -1:
            res.x, res.y, res.z = 0, 0, 0
        elif u == 0:
            res.x = self.radius*sin(v)
            res.y = self.radius*cos(v)
            res.z = 0
        else: # u == 1:
            res.x, res.y, res.z = 0, 0, self.height


cdef class Cylinder(ParametricSurface):

    def __init__(self, radius, height, closed=True, **kwds):
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius
        self.height = height
        self.closed = closed

    def x3d_geometry(self):
        return "<Cylinder radius='%s' height='%s'/>"%(self.radius, self.height)

    def tachyon_repr(self, render_params):
        transform = render_params.transform
        if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])):
            # Tachyon can't do sqashed
            return ParametricSurface.tachyon_repr(self, render_params)

        if transform is None:
            base = (0,0,0)
            top = (0,0,self.height)
            rad = self.radius
        else:
            base = transform.transform_point((0,0,0))
            top = transform.transform_point((0,0,self.height))
            radv = transform.transform_vector((self.radius,0,0))
            rad = sqrt(sum([x*x for x in radv]))
        cyl = """FCylinder
   Base %s %s %s
   Apex %s %s %s
   Rad %s
   %s     """%(base[0], base[1], base[2], top[0], top[1], top[2], rad, self.texture.id)
        if self.closed:
            normal = transform.transform_vector((0,0,1))
            base_cap = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % (base[0], base[1], base[2], normal[0], normal[1], normal[2], rad, self.texture.id)
            top_cap  = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % ( top[0],  top[1],  top[2], normal[0], normal[1], normal[2], rad, self.texture.id)
            return [base_cap, cyl, top_cap]
        else:
            return cyl

    def get_grid(self, ds):
        twoPi = 2*RDF.pi()
        v_res = min(max(int(twoPi*self.radius/ds), 5), 37)
        if self.closed:
            urange = [2,1,-1,-2]
        else:
            urange = [1,-1]
        vrange = [float(twoPi*k/v_res) for k in range(v_res)] + [0.0]
        return urange, vrange

    cdef eval_c(self, point_c *res, double u, double v):
        if u == -2:
            res.x, res.y, res.z = 0, 0, 0
        elif u == -1:
            res.x = self.radius*sin(v)
            res.y = self.radius*cos(v)
            res.z = 0
        elif u == 1:
            res.x = self.radius*sin(v)
            res.y = self.radius*cos(v)
            res.z = self.height
        else: # u == 2:
            res.x, res.y, res.z = 0, 0, self.height


def LineSegment(start, end, thickness=1, radius=None, **kwds):
    """
    Create a line segment, which is drawn as a cylinder from start to
    end with radius radius.

    EXAMPLES:
        sage: from sage.plot.plot3d.shapes import LineSegment, Sphere
        sage: P = (0,0,0.1)
        sage: Q = (0.5,0.6,0.7)
        sage: S = Sphere(.2, color='red').translate(P) + \
                  Sphere(.2, color='blue').translate(Q) + \
                  LineSegment(P, Q, .05, color='black')
        sage: S.show()
        sage: S = Sphere(.1, color='red').translate(P) + \
                  Sphere(.1, color='blue').translate(Q) + \
                  LineSegment(P, Q, .15, color='black')
        sage: S.show()

    AUTHOR:
        -- Robert Bradshaw
    """
    if radius is None:
        radius = thickness/50.0
    start = vector(RDF, start, sparse=False)
    end = vector(RDF, end, sparse=False)
    zaxis = vector(RDF, (0,0,1), sparse=False)
    diff = end - start
    height = sqrt(diff.dot_product(diff))
    cyl = Cylinder(radius, height, **kwds)
    axis = zaxis.cross_product(diff)
    if axis == 0:
        return cyl.translate(start)
    else:
        theta = -acos(diff[2]/height)
        return cyl.rotate(axis, theta).translate(start)

def Arrow(start, end, thickness=1, radius=None, head_radius=None, head_len=None, **kwds):
    if radius is None:
        radius = thickness/50.0
    if head_radius == None:
        head_radius = 3*radius
    if head_len == None:
        head_len = 3*head_radius
    start = vector(RDF, start, sparse=False)
    end = vector(RDF, end, sparse=False)
    zaxis = vector(RDF, (0,0,1), sparse=False)
    diff = end - start
    height = sqrt(diff.dot_product(diff))
    arrow = Cylinder(radius, height-head_len, **kwds) \
          + Cone(head_radius, head_len, **kwds).translate(0, 0, height-head_len)
    axis = zaxis.cross_product(diff)
    if axis == 0:
        return arrow.translate(start)
    else:
        theta = -acos(diff[2]/height)
        return arrow.rotate(axis, theta).translate(start)



cdef class Sphere(ParametricSurface):

    def __init__(self, radius, **kwds):
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius

    def x3d_geometry(self):
        return "<Sphere radius='%s'/>"%(self.radius)

    def tachyon_repr(self, render_params):
        transform = render_params.transform
        if not (transform is None or transform.is_uniform()):
            return ParametricSurface.tachyon_repr(self, render_params)

        if transform is None:
            cen = (0,0,0)
            rad = self.radius
        else:
            cen = transform.transform_point((0,0,0))
            radv = transform.transform_vector((self.radius,0,0))
            rad = sqrt(sum([x*x for x in radv]))
        return "Sphere center %s %s %s Rad %s %s" % (cen[0], cen[1], cen[2], rad, self.texture.id)

    def jmol_repr(self, render_params):
        name = render_params.unique_name('sphere')
        transform = render_params.transform
        if not (transform is None or transform.is_uniform()):
            return ParametricSurface.jmol_repr(self, render_params)

        if transform is None:
            cen = (0,0,0)
            rad = self.radius
        else:
            cen = transform.transform_point((0,0,0))
            radv = transform.transform_vector((self.radius,0,0))
            rad = sqrt(sum([x*x for x in radv]))
        return ["isosurface %s center {%s %s %s} sphere %s\n%s" % (name, cen[0], cen[1], cen[2], rad, self.texture.jmol_str("isosurface"))]

    def get_grid(self, ds):
        pi = RDF.pi()
        twoPi = 2*pi
        u_res = min(max(int(RDF.pi()*self.radius/ds), 6), 20)
        v_res = min(max(int(twoPi*self.radius/ds), 6), 36)
        urange = [-10] + [pi*k/u_res - twoPi/4 for k in range(1,u_res)] + [10]
        vrange = [float(twoPi*k/v_res) for k in range(v_res)] + [0.0]
        return urange, vrange

    cdef eval_c(self, point_c *res, double u, double v):
        if u == -10:
            res.x, res.y, res.z = 0, 0, -self.radius
        elif u == 10:
            res.x, res.y, res.z = 0, 0,  self.radius
        else:
            res.x = self.radius*cos(v) * cos(u)
            res.y = self.radius*sin(v) * cos(u)
            res.z = self.radius * sin(u)


cdef class Torus(ParametricSurface):
# e.g  show(sum([Torus(1,.03,20,20, color=[1, float(t/30), 0]).rotate((1,1,1),t) for t in range(30)], Sphere(.3)))
    """
    INPUT:
        R -- (default: 1) outer radius
        r -- (default: .3) inner radius

    OUTPUT:
        a 3d torus
    """
    def __init__(self, R=1, r=.3, **kwds):
        ParametricSurface.__init__(self, None, **kwds)
        self.R = R
        self.r = r

    def get_grid(self, ds):
        twoPi = 2*RDF.pi()
        u_divs = min(max(int(twoPi*self.R/ds), 6), 37)
        v_divs = min(max(int(twoPi*self.r/ds), 6), 37)
        urange = [-twoPi*k/u_divs for k in range(u_divs)] + [0.0]
        vrange = [ twoPi*k/v_divs for k in range(v_divs)] + [0.0]
        return urange, vrange

    cdef eval_c(self, point_c *res, double u, double v):
        res.x = (self.R+self.r*sin(v))*sin(u)
        res.y = (self.R+self.r*sin(v))*cos(u)
        res.z = self.r*cos(v)


class Text(PrimativeObject):
    def __init__(self, string, **kwds):
        PrimativeObject.__init__(self, **kwds)
        self.string = string
    def x3d_geometry(self):
        return "<Text string='%s' solid='true'/>"%self.string


