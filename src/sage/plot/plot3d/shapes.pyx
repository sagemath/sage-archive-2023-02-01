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

from sage.misc.all import srange

from base import Graphics3dGroup, Graphics3d


class Box(IndexFaceSet):
    """
    EXAMPLES:
        sage: from sage.plot.plot3d.shapes import Box

    A square black box:
        sage: show(Box([1,1,1]))

    A red rectangular box.
        sage: show(Box([2,3,4], color="red"))

    A stack of boxes:
        sage: show(sum([Box([2,3,1], color="red").translate((0,0,6*i)) for i in [0..3]]))

    A sinusoidal stack of multicolored boxes:
        sage: B = sum([Box([2,4,1/4], color=(i/4,i/5,1)).translate((sin(i),0,5-i)) for i in [0..20]])
        sage: show(B, figsize=6)
    """
    def __init__(self, *size, **kwds):
        if isinstance(size[0], (tuple, list)):
            from shapes2 import validate_frame_size
            size = validate_frame_size(size[0])
        self.size = size
        x, y, z = self.size
        faces = [[(x, y, z), (-x, y, z), (-x,-y, z), ( x,-y, z)],
                 [(x, y, z), ( x, y,-z), (-x, y,-z), (-x, y, z)],
                 [(x, y, z), ( x,-y, z), ( x,-y,-z), ( x, y,-z)] ]
        faces += [list(reversed([(-x,-y,-z) for x,y,z in face])) for face in faces]
        IndexFaceSet.__init__(self, faces, enclosed=True, **kwds)

    def bounding_box(self):
        """
        EXAMPLES:
            sage: from sage.plot.plot3d.shapes import Box
            sage: Box([1,2,3]).bounding_box()
            ((-1.0, -2.0, -3.0), (1.0, 2.0, 3.0))
        """
        return tuple([-a for a in self.size]), tuple(self.size)

    def x3d_geometry(self):
        return "<Box size='%s %s %s'/>"%self.size

def ColorCube(size, colors, opacity=1, **kwds):
    """
    Return a cube with given size and sides with given colors.

    INPUT:
        size -- 3-tuple of sizes (same as for box and frame)
        colors -- a list of either 3 or 6 color
        opacity -- (default: 1) opacity of cube sides
        **kwds -- passed to the face constructor

    OUTPUT:
        -- a 3d graphics object

    EXAMPLES:
    A color cube with translucent sides:
        sage: from sage.plot.plot3d.shapes import ColorCube
        sage: c = ColorCube([1,2,3], ['red', 'blue', 'green', 'black', 'white', 'orange'], opacity=0.5)
        sage: c.show()
        sage: list(c.texture_set())[0].opacity
        0.500000000000000

    If you omit the last 3 colors then the first three are repeated.
        sage: c = ColorCube([0.5,0.5,0.5], ['red', 'blue', 'green'])
    """
    if not isinstance(size, (tuple, list)):
        size = (size, size, size)
    box = Box(size)
    faces = box.face_list()
    if len(colors) == 3:
        colors = colors * 2
    all = []

    from texture import Texture
    for k in range(6):
        all.append(IndexFaceSet([faces[k]], enclosed=True,
             texture=Texture(colors[k], opacity=opacity),
             **kwds))
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

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
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

    def bounding_box(self):
        return (-self.radius, -self.radius, 0), (self.radius, self.radius, self.height)

    def x3d_geometry(self):
        return "<Cylinder radius='%s' height='%s'/>"%(self.radius, self.height)

    def tachyon_repr(self, render_params):
        transform = render_params.transform
        if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])):
            # Tachyon can't do sqashed
            return ParametricSurface.tachyon_repr(self, render_params)

        base, top = self.get_endpoints(transform)
        rad = self.get_radius(transform)
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

    def jmol_repr(self, render_params):
        transform = render_params.transform
        base, top = self.get_endpoints(transform)
        rad = self.get_radius(transform)

        cdef double ratio = sqrt(rad*rad / ((base[0]-top[0])**2 + (base[1]-top[1])**2 + (base[2]-top[2])**2))
        #print ratio

        if ratio > .02:
            if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])) or ratio > .05:
                # Jmol can't do sqashed
                return ParametricSurface.jmol_repr(self, render_params)

        name = render_params.unique_name('line')
        return ["""
draw %s width %s {%s %s %s} {%s %s %s}\n%s
""" % (name,
       rad,
       base[0], base[1], base[2],
       top [0], top [1], top [2],
       self.texture.jmol_str("$" + name)) ]

    def get_endpoints(self, transform=None):
        if transform is None:
            return (0,0,0), (0,0,self.height)
        else:
            return transform.transform_point((0,0,0)), transform.transform_point((0,0,self.height))

    def get_radius(self, transform=None):
        if transform is None:
            return self.radius
        else:
            radv = transform.transform_vector((self.radius,0,0))
            return sqrt(sum([x*x for x in radv]))


    def get_grid(self, ds):
        twoPi = 2*RDF.pi()
        v_res = min(max(int(twoPi*self.radius/ds), 5), 37)
        if self.closed:
            urange = [2,1,-1,-2]
        else:
            urange = [1,-1]
        vrange = [float(twoPi*k/v_res) for k in range(v_res)] + [0.0]
        return urange, vrange

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
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
    end   = vector(RDF, end, sparse=False)
    zaxis = vector(RDF, (0,0,1), sparse=False)
    diff  = end - start
    height= sqrt(diff.dot_product(diff))
    cyl   = Cylinder(radius, height, **kwds)
    axis  = zaxis.cross_product(diff)
    if axis == 0:
        if diff[2] < 0:
            return cyl.translate(end)
        else:
            return cyl.translate(start)
    else:
        theta = -acos(diff[2]/height)
        return cyl.rotate(axis, theta).translate(start)

def arrow3d(start, end, thickness=1, radius=None, head_radius=None, head_len=None, **kwds):
    """
    Create a 3d arrow.

    INPUT:
        start -- (x,y,z) point; the starting point of the arrow
        end -- (x,y,z) point; the end point
        thickness -- (default: 1); how thick the arrow is
        radius -- (default: thickness/50.0) the radius of the arrow
        head_radius -- (default: 3*radius); radius of arrow head
        head_len -- (default: 3*head_radius); len of arrow head

    EXAMPLES:
    The default arrow:
        sage: arrow3d((0,0,0), (1,1,1), 1)

    A fat arrow:
        sage: arrow3d((0,0,0), (1,1,1), radius=0.1)

    A green arrow:
        sage: arrow3d((0,0,0), (1,1,1), color='green')

    A fat arrow head:
        sage: arrow3d((2,1,0), (1,1,1), color='green', head_radius=0.3, aspect_ratio=[1,1,1])

    Many arrow arranged in a circle (flying spears?):
        sage: sum([arrow3d((cos(t),sin(t),0),(cos(t),sin(t),1)) for t in [0,0.3,..,2*pi]])

    TESTS:
    If the arrow is too long, the shaft and part of the head is cut off.
        sage: a = arrow3d((0,0,0), (0,0,0.5), head_len=1)
        sage: len(a.all)
        1
        sage: type(a.all[0])
        <type 'sage.plot.plot3d.shapes.Cone'>

    Arrows are always constructed pointing up in the z direction from
    the origin, and then rotated/translated into place. This works for
    every arrow direction except the -z direction.  We take care of the
    anomoly by testing to see if the arrow should point in the -z
    direction, and if it should, just scaling the constructed arrow by
    -1 (i.e., every point is sent to its negative). The scaled arrow
    then points downwards. The doctest just tests that the scale of -1
    is applied to the arrow.

        sage: a = arrow3d((0,0,0), (0,0,-1))
        sage: a.all[0].get_transformation().transform_point((0,0,1))
        (0.0, 0.0, -1.0)
    """
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
    length = sqrt(diff.dot_product(diff))
    if length <= head_len:
        arrow = Cone(head_radius*length/head_len, length, **kwds)
    else:
        arrow = Cylinder(radius, length-head_len, **kwds) \
                + Cone(head_radius, head_len, **kwds).translate(0, 0, length-head_len)
    axis = zaxis.cross_product(diff)
    if axis == 0:
        if diff[2] >= 0:
            return arrow.translate(start)
        else:
            return arrow.scale(-1).translate(start)
    else:
        theta = -acos(diff[2]/length)
        return arrow.rotate(axis, theta).translate(start)



cdef class Sphere(ParametricSurface):

    def __init__(self, radius, **kwds):
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius

    def bounding_box(self):
        """
        Return the bounding box that contains this sphere.

        EXAMPLES:
            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(3).bounding_box()
            ((-3.0, -3.0, -3.0), (3.0, 3.0, 3.0))
        """
        return ((-self.radius, -self.radius, -self.radius),
                (self.radius, self.radius, self.radius))

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
        if rad < 0.5:
            res = "resolution %s" % min(int(7/rad), 100)
        else:
            res = ""
        return ["isosurface %s %s center {%s %s %s} sphere %s\n%s" % (name, res, cen[0], cen[1], cen[2], rad, self.texture.jmol_str("isosurface"))]

    def get_grid(self, ds):
        pi = RDF.pi()
        twoPi = 2*pi
        u_res = min(max(int(RDF.pi()*self.radius/ds), 6), 20)
        v_res = min(max(int(twoPi*self.radius/ds), 6), 36)
        urange = [-10] + [pi*k/u_res - twoPi/4 for k in range(1,u_res)] + [10]
        vrange = [float(twoPi*k/v_res) for k in range(v_res)] + [0.0]
        return urange, vrange

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
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

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
        res.x = (self.R+self.r*sin(v))*sin(u)
        res.y = (self.R+self.r*sin(v))*cos(u)
        res.z = self.r*cos(v)


class Text(PrimitiveObject):
    def __init__(self, string, **kwds):
        PrimitiveObject.__init__(self, **kwds)
        self.string = string

    def x3d_geometry(self):
        return "<Text string='%s' solid='true'/>"%self.string

    def obj_repr(self, render_params):
        return ''

    def tachyon_repr(self, render_params):
        return ''
        # Text in Tachyon not implemented yet.
        # I have no idea what the code below is supposed to do.
##         transform = render_params.transform
##         if not (transform is None or transform.is_uniform()):
##             return ParametricSurface.tachyon_repr(self, render_params)

##         if transform is None:
##             cen = (0,0,0)
##             rad = self.radius
##         else:
##             cen = transform.transform_point((0,0,0))
##             radv = transform.transform_vector((self.radius,0,0))
##             rad = sqrt(sum([x*x for x in radv]))
##         return "Sphere center %s %s %s Rad %s %s" % (cen[0], cen[1], cen[2], rad, self.texture.id)

    def jmol_repr(self, render_params):
        cen = render_params.transform.transform_point((0,0,0))
        render_params.atom_list.append(cen)
        atom_no = len(render_params.atom_list)
        return ['select atomno = %s' % atom_no,
                self.get_texture().jmol_str("atom"),
                'label "%s"' % self.string] #.replace('\n', '|')]

    def bounding_box(self):
        return (0,0,0), (0,0,0)

