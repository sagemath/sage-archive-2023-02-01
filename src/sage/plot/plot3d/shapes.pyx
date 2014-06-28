"""
Basic objects such as Sphere, Box, Cone, etc.

AUTHORS:
    -- Robert Bradshaw 2007-02: initial version
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
    double M_PI

from sage.rings.real_double  import RDF
from sage.modules.free_module_element import vector

from sage.misc.all import srange
from sage.plot.misc import rename_keyword

from base import Graphics3dGroup, Graphics3d

# Helper function to check that Box input is right
def validate_frame_size(size):
    """
    Checks that the input is an iterable of length 3 with all
    elements nonnegative and coercible to floats.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import validate_frame_size
        sage: validate_frame_size([3,2,1])
        [3.0, 2.0, 1.0]

    TESTS::

        sage: from sage.plot.plot3d.shapes import validate_frame_size
        sage: validate_frame_size([3,2,-1])
        Traceback (most recent call last):
        ...
        ValueError: each box dimension must be nonnegative
        sage: validate_frame_size([sqrt(-1),3,2])
        Traceback (most recent call last):
        ...
        TypeError: each box dimension must coerce to a float
    """
    if not isinstance(size, (list, tuple)):
        raise TypeError("size must be a list or tuple")
    if len(size) != 3:
        raise TypeError("size must be of length 3")
    try:
        size = [float(x) for x in size]
    except TypeError:
        raise TypeError("each box dimension must coerce to a float")
    for x in size:
        if x < 0:
            raise ValueError("each box dimension must be nonnegative")
    return size


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
        """
        EXAMPLES:
            sage: from sage.plot.plot3d.shapes import Box
            sage: Box(10, 1, 1) + Box(1, 10, 1) + Box(1, 1, 10)
        """
        if isinstance(size[0], (tuple, list)):
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
        """
        EXAMPLES:
            sage: from sage.plot.plot3d.shapes import Box
            sage: Box([1,2,1/4]).x3d_geometry()
            "<Box size='1.0 2.0 0.25'/>"
        """
        return "<Box size='%s %s %s'/>" % tuple(self.size)

def ColorCube(size, colors, opacity=1, **kwds):
    """
    Return a cube with given size and sides with given colors.

    INPUT:
        size -- 3-tuple of sizes (same as for box and frame)
        colors -- a list of either 3 or 6 colors
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

    If you omit the last 3 colors then the first three are repeated (with
    repeated colors on opposing faces):
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
    """
    A cone, with base in the xy-plane pointing up the z-axis.

    INPUT:

        - ``radius`` - positive real number

        - ``height`` - positive real number

        - ``closed`` - whether or not to include the base (default True)

        - ``**kwds`` -- passed to the ParametricSurface constructor

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Cone
        sage: c = Cone(3/2, 1, color='red') + Cone(1, 2, color='yellow').translate(3, 0, 0)
        sage: c.show(aspect_ratio=1)

    We may omit the base::

        sage: Cone(1, 1, closed=False)

    A spiky plot of the sine function::

        sage: sum(Cone(.1, sin(n), color='yellow').translate(n, sin(n), 0) for n in [0..10, step=.1])

    A Christmas tree::

        sage: T = sum(Cone(exp(-n/5), 4/3*exp(-n/5), color=(0, .5, 0)).translate(0, 0, -3*exp(-n/5)) for n in [1..7])
        sage: T += Cone(1/8, 1, color='brown').translate(0, 0, -3)
        sage: T.show(aspect_ratio=1, frame=False)
    """
    def __init__(self, radius, height, closed=True, **kwds):
        """
        TESTS:
            sage: from sage.plot.plot3d.shapes import Cone
            sage: c = Cone(1/2, 1, opacity=.5)
        """
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius
        self.height = height
        self.closed = closed

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: Cone(1, 3).x3d_geometry()
            "<Cone bottomRadius='1.0' height='3.0'/>"

        """
        return "<Cone bottomRadius='%s' height='%s'/>"%(self.radius, self.height)

    def get_grid(self, ds):
        """
        Returns the grid on which to evaluate this parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: Cone(1, 3, closed=True).get_grid(100)
            ([1, 0, -1], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: Cone(1, 3, closed=False).get_grid(100)
            ([1, 0], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: len(Cone(1, 3).get_grid(.001)[1])
            38
        """
        cdef int k, t_res = min(max(int(2*M_PI*self.radius/ds), 5), 37)
        if self.closed:
            urange = [1,0,-1]
        else:
            urange = [1,0]
        vrange = [2*M_PI*k/t_res for k from 0 <= k < t_res] + [0.0]
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
    """
    A cone, with base in the xy-plane pointing up the z-axis.

    INPUT:

        - ``radius`` - positive real number

        - ``height`` - positive real number

        - ``closed`` - whether or not to include the ends (default True)

        - ``**kwds`` -- passed to the ParametricSurface constructor

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Cylinder
        sage: c = Cylinder(3/2, 1, color='red') + Cylinder(1, 2, color='yellow').translate(3, 0, 0)
        sage: c.show(aspect_ratio=1)

    We may omit the base::

        sage: Cylinder(1, 1, closed=False)

    Some gears::

        sage: G = Cylinder(1, .5) + Cylinder(.25, 3).translate(0, 0, -3)
        sage: G += sum(Cylinder(.2, 1).translate(cos(2*pi*n/9), sin(2*pi*n/9), 0) for n in [1..9])
        sage: G += G.translate(2.3, 0, -.5)
        sage: G += G.translate(3.5, 2, -1)
        sage: G.show(aspect_ratio=1, frame=False)
    """
    def __init__(self, radius, height, closed=True, **kwds):
        """
        TESTS:
            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 1, color='red')
        """
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius
        self.height = height
        self.closed = closed

    def bounding_box(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 2).bounding_box()
            ((-1.0, -1.0, 0), (1.0, 1.0, 2.0))

        """
        return (-self.radius, -self.radius, 0), (self.radius, self.radius, self.height)

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 2).x3d_geometry()
            "<Cylinder radius='1.0' height='2.0'/>"

        """
        return "<Cylinder radius='%s' height='%s'/>"%(self.radius, self.height)

    def tachyon_repr(self, render_params):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: C = Cylinder(1/2, 4, closed=False)
            sage: C.tachyon_repr(C.default_render_params())
            'FCylinder\n   Base 0 0 0\n   Apex 0 0 4.0\n   Rad 0.5\n   texture...     '
            sage: C = Cylinder(1, 2)
            sage: C.tachyon_repr(C.default_render_params())
                ['Ring Center 0 0 0 Normal 0 0 1 Inner 0 Outer 1.0 texture...',
                 'FCylinder\n   Base 0 0 0\n   Apex 0 0 2.0\n   Rad 1.0\n   texture...     ',
                 'Ring Center 0 0 2.0 Normal 0 0 1 Inner 0 Outer 1.0 texture...']
        """
        transform = render_params.transform
        if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])):
            # Tachyon can't do squashed
            return ParametricSurface.tachyon_repr(self, render_params)

        base, top = self.get_endpoints(transform)
        rad = self.get_radius(transform)
        cyl = """FCylinder
   Base %s %s %s
   Apex %s %s %s
   Rad %s
   %s     """%(base[0], base[1], base[2], top[0], top[1], top[2], rad, self.texture.id)
        if self.closed:
            normal = (0,0,1)
            if transform is not None:
                normal = transform.transform_vector(normal)
            base_cap = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % (base[0], base[1], base[2], normal[0], normal[1], normal[2], rad, self.texture.id)
            top_cap  = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % ( top[0],  top[1],  top[2], normal[0], normal[1], normal[2], rad, self.texture.id)
            return [base_cap, cyl, top_cap]
        else:
            return cyl

    def jmol_repr(self, render_params):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder

        For thin cylinders, lines are used::
            sage: C = Cylinder(.1, 4)
            sage: C.jmol_repr(C.default_render_params())
            ['\ndraw line_1 width 0.1 {0 0 0} {0 0 4.0}\ncolor $line_1  [102,102,255]\n']

        For anything larger, we use a pmesh::
            sage: C = Cylinder(3, 1, closed=False)
            sage: C.jmol_repr(C.testing_render_params())
            ['pmesh obj_1 "obj_1.pmesh"\ncolor pmesh  [102,102,255]']
        """
        transform = render_params.transform
        base, top = self.get_endpoints(transform)
        rad = self.get_radius(transform)

        cdef double ratio = sqrt(rad*rad / ((base[0]-top[0])**2 + (base[1]-top[1])**2 + (base[2]-top[2])**2))

        if ratio > .02:
            if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])) or ratio > .05:
                # Jmol can't do squashed
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
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: from sage.plot.plot3d.transform import Transformation
            sage: Cylinder(1, 5).get_endpoints()
            ((0, 0, 0), (0, 0, 5.0))
            sage: Cylinder(1, 5).get_endpoints(Transformation(trans=(1,2,3), scale=(2,2,2)))
            ((1.0, 2.0, 3.0), (1.0, 2.0, 13.0))
        """
        if transform is None:
            return (0,0,0), (0,0,self.height)
        else:
            return transform.transform_point((0,0,0)), transform.transform_point((0,0,self.height))

    def get_radius(self, transform=None):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: from sage.plot.plot3d.transform import Transformation
            sage: Cylinder(3, 1).get_radius()
            3.0
            sage: Cylinder(3, 1).get_radius(Transformation(trans=(1,2,3), scale=(2,2,2)))
            6.0
        """
        if transform is None:
            return self.radius
        else:
            radv = transform.transform_vector((self.radius,0,0))
            return sqrt(sum([x*x for x in radv]))

    def get_grid(self, ds):
        """
        Returns the grid on which to evaluate this parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 3, closed=True).get_grid(100)
            ([2, 1, -1, -2], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: Cylinder(1, 3, closed=False).get_grid(100)
            ([1, -1], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: len(Cylinder(1, 3).get_grid(.001)[1])
            38
        """
        cdef int k, v_res = min(max(int(2*M_PI*self.radius/ds), 5), 37)
        if self.closed:
            urange = [2,1,-1,-2]
        else:
            urange = [1,-1]
        vrange = [2*M_PI*k/v_res for k from 0 <= k < v_res] + [0.0]
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

@rename_keyword(deprecation=7154, deprecated_option='thickness', thickness='width')
def arrow3d(start, end, width=1, radius=None, head_radius=None, head_len=None, **kwds):
    """
    Create a 3d arrow.

    INPUT:
        start -- (x,y,z) point; the starting point of the arrow
        end -- (x,y,z) point; the end point
        width -- (default: 1); how wide the arrow is
        radius -- (default: width/50.0) the radius of the arrow
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

    Change the width of the arrow. (Note: for an arrow that scales with zoom, please consider
    the 'line3d' function with the option 'arrow_head=True'):
        sage: arrow3d((0,0,0), (1,1,1), width=1)

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
    anomaly by testing to see if the arrow should point in the -z
    direction, and if it should, just scaling the constructed arrow by
    -1 (i.e., every point is sent to its negative). The scaled arrow
    then points downwards. The doctest just tests that the scale of -1
    is applied to the arrow.

        sage: a = arrow3d((0,0,0), (0,0,-1))
        sage: a.all[0].get_transformation().transform_point((0,0,1))
        (0.0, 0.0, -1.0)

    The thickness option is now deprecated.  It has been replaced by the width option.

        sage: arrow3d((0,0,0), (1,1,1), thickness=1)
        doctest:...: DeprecationWarning: use the option 'width' instead of 'thickness'
        See http://trac.sagemath.org/7154 for details.
        <BLANKLINE>
    """
    if radius is None:
        radius = width/50.0
    if head_radius is None:
        head_radius = 3*radius
    if head_len is None:
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
    """
    This class represents a sphere centered at the origin.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Sphere
        sage: Sphere(3)

    Plot with aspect_ratio=1 to see it unsquashed::

        sage: S = Sphere(3, color='blue') + Sphere(2, color='red').translate(0,3,0)
        sage: S.show(aspect_ratio=1)

    Scale to get an ellipsoid::
        sage: S = Sphere(1).scale(1,2,1/2)
        sage: S.show(aspect_ratio=1)

    """
    def __init__(self, radius, **kwds):
        """
        TESTS:
            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(3)
        """
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius

    def bounding_box(self):
        """
        Return the bounding box that contains this sphere.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(3).bounding_box()
            ((-3.0, -3.0, -3.0), (3.0, 3.0, 3.0))
        """
        return ((-self.radius, -self.radius, -self.radius),
                (self.radius, self.radius, self.radius))

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(12).x3d_geometry()
            "<Sphere radius='12.0'/>"
        """
        return "<Sphere radius='%s'/>"%(self.radius)

    def tachyon_repr(self, render_params):
        """
        Tachyon can natively handle spheres. Ellipsoids rendering is done
        as a parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: S = Sphere(2)
            sage: S.tachyon_repr(S.default_render_params())
            'Sphere center 0 0 0 Rad 2.0 texture...'
            sage: S.translate(1, 2, 3).scale(3).tachyon_repr(S.default_render_params())
            [['Sphere center 3.0 6.0 9.0 Rad 6.0 texture...']]
            sage: S.scale(1,1/2,1/4).tachyon_repr(S.default_render_params())
            [['TRI V0 0 0 -0.5 V1 0.308116 0.0271646 -0.493844 V2 0.312869 0 -0.493844',
              'texture...',
               ...
              'TRI V0 0.308116 -0.0271646 0.493844 V1 0.312869 0 0.493844 V2 0 0 0.5',
              'texture...']]
        """
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
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere

        Jmol has native code for handling spheres::

            sage: S = Sphere(2)
            sage: S.jmol_repr(S.default_render_params())
            ['isosurface sphere_1  center {0 0 0} sphere 2.0\ncolor isosurface  [102,102,255]']
            sage: S.translate(10, 100, 1000).jmol_repr(S.default_render_params())
            [['isosurface sphere_1  center {10.0 100.0 1000.0} sphere 2.0\ncolor isosurface  [102,102,255]']]

        It can't natively handle ellipsoids::

            sage: Sphere(1).scale(2, 3, 4).jmol_repr(S.testing_render_params())
            [['pmesh obj_2 "obj_2.pmesh"\ncolor pmesh  [102,102,255]']]

        Small spheres need extra hints to render well::

            sage: Sphere(.01).jmol_repr(S.default_render_params())
            ['isosurface sphere_1 resolution 100 center {0 0 0} sphere 0.01\ncolor isosurface  [102,102,255]']
        """
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

    def get_grid(self, double ds):
        """
        Returns the the range of variables to be evaluated on to render as a
        parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(1).get_grid(100)
            ([-10.0, ..., 0.0, ..., 10.0],
             [0.0, ..., 3.141592653589793, ..., 0.0])
        """
        cdef int K, u_res, v_res
        u_res = min(max(int(M_PI*self.radius/ds), 6), 20)
        v_res = min(max(int(2*M_PI * self.radius/ds), 6), 36)
        urange = [-10.0] + [M_PI * k/u_res - M_PI/2 for k in range(1, u_res)] + [10.0]
        vrange = [2*M_PI * k/v_res for k in range(v_res)] + [0.0]
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

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Torus
        sage: Torus(1, .2).show(aspect_ratio=1)
        sage: Torus(1, .7, color='red').show(aspect_ratio=1)

    A rubberband ball::

        sage: show(sum([Torus(1, .03, color=(1, t/30.0, 0)).rotate((1,1,1),t) for t in range(30)]))

    Mmm... doughnuts::

        sage: D = Torus(1, .4, color=(.5, .3, .2)) + Torus(1, .3, color='yellow').translate(0, 0, .15)
        sage: G = sum(D.translate(RDF.random_element(-.2, .2), RDF.random_element(-.2, .2), .8*t) for t in range(10))
        sage: G.show(aspect_ratio=1, frame=False)
    """
    def __init__(self, R=1, r=.3, **kwds):
        """
        TESTS:
            sage: from sage.plot.plot3d.shapes import Torus
            sage: T = Torus(1, .5)
        """
        ParametricSurface.__init__(self, None, **kwds)
        self.R = R
        self.r = r

    def get_grid(self, ds):
        """
        Returns the the range of variables to be evaluated on to render as a
        parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: Torus(2, 1).get_grid(100)
            ([0.0, -1.047..., -3.141592653589793, ..., 0.0],
             [0.0, 1.047..., 3.141592653589793, ..., 0.0])
        """
        cdef int k, u_divs, v_divs
        u_divs = min(max(int(4*M_PI * self.R/ds), 6), 37)
        v_divs = min(max(int(4*M_PI * self.r/ds), 6), 37)
        urange = [0.0] + [-2*M_PI * k/u_divs for k in range(1, u_divs)] + [0.0]
        vrange = [ 2*M_PI * k/v_divs for k in range(v_divs)] + [0.0]
        return urange, vrange

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
        res.x = (self.R+self.r*sin(v))*sin(u)
        res.y = (self.R+self.r*sin(v))*cos(u)
        res.z = self.r*cos(v)


class Text(PrimitiveObject):
    """
    A text label attached to a point in 3d space. It always starts at the
    origin, translate it to move it elsewhere.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Text
        sage: Text("Just a lonely label.")
        sage: pts = [(RealField(10)^3).random_element() for k in range(20)]
        sage: sum(Text(str(P)).translate(P) for P in pts)

    """
    def __init__(self, string, **kwds):
        """
        TEST:
            sage: from sage.plot.plot3d.shapes import Text
            sage: T = Text("Hi")
        """
        PrimitiveObject.__init__(self, **kwds)
        self.string = string

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").x3d_geometry()
            "<Text string='Hi' solid='true'/>"
        """
        return "<Text string='%s' solid='true'/>"%self.string

    def obj_repr(self, render_params):
        """
        The obj file format doesn't support text strings::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").obj_repr(None)
            ''
        """
        return ''

    def tachyon_repr(self, render_params):
        """
        Strings are not yet supported in Tachyon, so we ignore them for now::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").tachyon_repr(None)
            ''
        """
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
        """
        Labels in jmol must be attached to atoms.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Text
            sage: T = Text("Hi")
            sage: T.jmol_repr(T.testing_render_params())
            ['select atomno = 1', 'color atom  [102,102,255]', 'label "Hi"']
            sage: T = Text("Hi").translate(-1, 0, 0) + Text("Bye").translate(1, 0, 0)
            sage: T.jmol_repr(T.testing_render_params())
            [[['select atomno = 1', 'color atom  [102,102,255]', 'label "Hi"']],
             [['select atomno = 2', 'color atom  [102,102,255]', 'label "Bye"']]]
        """
        cen = (0,0,0)
        if render_params.transform is not None:
            cen = render_params.transform.transform_point(cen)
        render_params.atom_list.append(cen)
        atom_no = len(render_params.atom_list)
        return ['select atomno = %s' % atom_no,
                self.get_texture().jmol_str("atom"),
                'label "%s"' % self.string] #.replace('\n', '|')]

    def bounding_box(self):
        """
        Text labels have no extent::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").bounding_box()
            ((0, 0, 0), (0, 0, 0))
        """
        return (0,0,0), (0,0,0)

