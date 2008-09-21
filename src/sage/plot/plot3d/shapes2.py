"""
Lines, Frames, Spheres, Points, Dots, and Text
"""

import math
import shapes

from base import PrimitiveObject, point_list_bounding_box

from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.misc.misc import srange

from texture import Texture

TACHYON_PIXEL = 1/200.0

from shapes import Text, Sphere

def line3d(points, thickness=1, radius=None, arrow_head=False, **kwds):
    r"""
    Draw a 3d line joining a sequence of points.

    One may specify either a thickness or radius. If a thickness
    is specified, this line will have a constant diameter regardless
    of scaling and zooming. If a radius is specified, it will
    behave as a series of cylinders.

    INPUT:
        points -- a list of at least 2 points
        thickness -- (default: 1)
        radius -- (default: None)
        arrow_head -- (default: False)
        color -- a word that describes a color
        rgbcolor -- (r,g,b) with r, g, b between 0 and 1 that describes a color
        opacity -- (default: 1) if less than 1 then is transparent

    EXAMPLES:
    A line in 3-space:
        sage: line3d([(1,2,3), (1,0,-2), (3,1,4), (2,1,-2)])

    The same line but red:
        sage: line3d([(1,2,3), (1,0,-2), (3,1,4), (2,1,-2)], color='red')

    A transparent thick green line and a little blue line:
        sage: line3d([(0,0,0), (1,1,1), (1,0,2)], opacity=0.5, radius=0.1, \
                     color='green') + line3d([(0,1,0), (1,0,2)])

    A Dodecahedral complex of 5 tetrahedrons (a more elaborate examples from Peter Jipsen):
        sage: def tetra(col):
        ...       return line3d([(0,0,1), (2*sqrt(2.)/3,0,-1./3), (-sqrt(2.)/3, sqrt(6.)/3,-1./3),\
        ...              (-sqrt(2.)/3,-sqrt(6.)/3,-1./3), (0,0,1), (-sqrt(2.)/3, sqrt(6.)/3,-1./3),\
        ...              (-sqrt(2.)/3,-sqrt(6.)/3,-1./3), (2*sqrt(2.)/3,0,-1./3)],\
        ...              color=col, thickness=10, aspect_ratio=[1,1,1])
        ...
        sage: v  = (sqrt(5.)/2-5/6, 5/6*sqrt(3.)-sqrt(15.)/2, sqrt(5.)/3)
        sage: t  = acos(sqrt(5.)/3)/2
        sage: t1 = tetra('blue').rotateZ(t)
        sage: t2 = tetra('red').rotateZ(t).rotate(v,2*pi/5)
        sage: t3 = tetra('green').rotateZ(t).rotate(v,4*pi/5)
        sage: t4 = tetra('yellow').rotateZ(t).rotate(v,6*pi/5)
        sage: t5 = tetra('orange').rotateZ(t).rotate(v,8*pi/5)
        sage: show(t1+t2+t3+t4+t5, frame=False)


    """
    if len(points) < 2:
        raise ValueError, "there must be at least 2 points"
    for i in range(len(points)):
        x, y, z = points[i]
        points[i] = float(x), float(y), float(z)
    if radius is None:
        L = Line(points, thickness=thickness, arrow_head=arrow_head, **kwds)
        L._set_extra_kwds(kwds)
        return L
    else:
        v = []
        texture = Texture(kwds)
        if kwds.has_key('texture'):
            kwds = kwds.copy()
            del kwds['texture']
        for i in range(len(points) - 1):
            line = shapes.arrow3d if i == len(points)-2 and arrow_head else shapes.LineSegment
            v.append(line(points[i], points[i+1], texture=texture, radius=radius, **kwds))
        w = sum(v)
        w._set_extra_kwds(kwds)
        return w


def frame3d(lower_left, upper_right, **kwds):
    """
    Draw a frame in 3d.
    """
    x0,y0,z0 = lower_left
    x1,y1,z1 = upper_right
    L1 = line3d([(x0,y0,z0), (x0,y1,z0), (x1,y1,z0), (x1,y0,z0),  (x0,y0,z0), # top square
                 (x0,y0,z1), (x0,y1,z1), (x1,y1,z1), (x1,y0,z1),  (x0,y0,z1)],  # bottom square
                **kwds)
    # 3 additional lines joining top to bottom
    v2 = line3d([(x0,y1,z0), (x0,y1,z1)], **kwds)
    v3 = line3d([(x1,y0,z0), (x1,y0,z1)], **kwds)
    v4 = line3d([(x1,y1,z0), (x1,y1,z1)], **kwds)
    F  = L1 + v2 + v3 + v4
    F._set_extra_kwds(kwds)
    return F

def frame_labels(lower_left, upper_right,
                 label_lower_left, label_upper_right, eps = 1,
                 **kwds):
    x0,y0,z0 = lower_left
    x1,y1,z1 = upper_right
    lx0,ly0,lz0 = label_lower_left
    lx1,ly1,lz1 = label_upper_right

    from math import log
    log10 = log(10)
    nd = lambda a: int(log(a)/log10)
    def fmt_string(a):
        b = a/2.0
        if b >= 1:
            return "%.0f"
        n = max(0, 2 - nd(a/2.0))
        return "%%.%sf"%n

    fmt = fmt_string(lx1 - lx0)

    color = (0.3,0.3,0.3)
    T =  Text(fmt%lx0, color=color).translate((x0,y0-eps,z0))
    T += Text(fmt%avg(lx0,lx1), color=color).translate((avg(x0,x1),y0-eps,z0))
    T += Text(fmt%lx1, color=color).translate((x1,y0-eps,z0))

    fmt = fmt_string(ly1 - ly0)
    T += Text(fmt%ly0, color=color).translate((x1+eps,y0,z0))
    T += Text(fmt%avg(ly0,ly1), color=color).translate((x1+eps,avg(y0,y1),z0))
    T += Text(fmt%ly1, color=color).translate((x1+eps,y1,z0))

    fmt = fmt_string(lz1 - lz0)
    T += Text(fmt%lz0, color=color).translate((x0-eps,y0,z0))
    T += Text(fmt%avg(lz0,lz1), color=color).translate((x0-eps,y0,avg(z0,z1)))
    T += Text(fmt%lz1, color=color).translate((x0-eps,y0,z1))
    return T

def validate_frame_size(size):
    # Verify that the input is valid
    if not isinstance(size, (list, tuple)):
        raise TypeError, "size must be a list or tuple"
    if len(size) != 3:
        raise TypeError, "size must be of length 3"
    try:
        size = [float(x) for x in size]
    except TypeError:
        raise TypeError, "each box dimension must coerce to a float"
    for x in size:
        if x < 0:
            raise ValueError, "each box dimension must be nonnegative"
    return size


def ruler(start, end, ticks=4, sub_ticks=4, absolute=False, snap=False, **kwds):
    start = vector(RDF, start)
    end   = vector(RDF, end)
    dir = end - start
    dist = math.sqrt(dir.dot_product(dir))
    dir /= dist

    one_tick = dist/ticks * 1.414
    unit = 10 ** math.floor(math.log(dist/ticks, 10))
    if unit * 5 < one_tick:
        unit *= 5
    elif unit * 2 < one_tick:
        unit *= 2

    if dir[0]:
        tick = dir.cross_product(vector(RDF, (0,0,-dist/30)))
    elif dir[1]:
        tick = dir.cross_product(vector(RDF, (0,0,dist/30)))
    else:
        tick = vector(RDF, (dist/30,0,0))

    if snap:
        for i in range(3):
            start[i] = unit * math.floor(start[i]/unit + 1e-5)
            end[i] = unit * math.ceil(end[i]/unit - 1e-5)

    if absolute:
        if dir[0]*dir[1] or dir[1]*dir[2] or dir[0]*dir[2]:
            raise ValueError, "Absolute rulers only valid for axis-aligned paths"
        m = max(dir[0], dir[1], dir[2])
        if dir[0] == m:
            off = start[0]
        elif dir[1] == m:
            off = start[1]
        else:
            off = start[2]
        first_tick = unit * math.ceil(off/unit - 1e-5) - off
    else:
        off = 0
        first_tick = 0

    ruler = shapes.LineSegment(start, end, **kwds)
    for k in range(1, int(sub_ticks * first_tick/unit)):
        P = start + dir*(k*unit/sub_ticks)
        ruler += shapes.LineSegment(P, P + tick/2, **kwds)
    for d in srange(first_tick, dist + unit/(sub_ticks+1), unit):
        P = start + dir*d
        ruler += shapes.LineSegment(P, P + tick, **kwds)
        ruler += shapes.Text(str(d+off), **kwds).translate(P - tick)
        if dist - d < unit:
            sub_ticks = int(sub_ticks * (dist - d)/unit)
        for k in range(1, sub_ticks):
            P += dir * (unit/sub_ticks)
            ruler += shapes.LineSegment(P, P + tick/2, **kwds)
    return ruler

def ruler_frame(lower_left, upper_right, ticks=4, sub_ticks=4, **kwds):
    return ruler(lower_left, (upper_right[0], lower_left[1], lower_left[2]), ticks=ticks, sub_ticks=sub_ticks, absolute=True) \
         + ruler(lower_left, (lower_left[0], upper_right[1], lower_left[2]), ticks=ticks, sub_ticks=sub_ticks, absolute=True) \
         + ruler(lower_left, (lower_left[0], lower_left[1], upper_right[2]), ticks=ticks, sub_ticks=sub_ticks, absolute=True)


def avg(a,b):
    return (a+b)/2.0



###########################


def sphere(center=(0,0,0), size=1, **kwds):
    r"""
    Return a plot of a sphere of radius size centered at $(x,y,z)$.

    INPUT:
       (x,y,z) -- center (default: (0,0,0)
       size -- the radius (default: 1)

    EXAMPLES:
    A simple sphere:
       sage: sphere()

    Two spheres touching:
       sage: sphere(center=(-1,0,0)) + sphere(center=(1,0,0), aspect_ratio=[1,1,1])

    Spheres of radii 1 and 2 one stuck into the other:
       sage: sphere(color='orange') + sphere(color=(0,0,0.3), \
                    center=(0,0,-2),size=2,opacity=0.9)

    We draw a transparent sphere on a saddle.
       sage: u,v = var('u v')
       sage: saddle = plot3d(u^2 - v^2, (u,-2,2), (v,-2,2))
       sage: sphere((0,0,1), color='red', opacity=0.5, aspect_ratio=[1,1,1]) + saddle
    """
    G = Sphere(size, texture=Texture(kwds), **kwds)
    H = G.translate(center)
    H._set_extra_kwds(kwds)
    return H

def text3d(txt, (x,y,z), **kwds):
    r"""
    Display 3d text.

    INPUT:
        txt -- some text
        (x,y,z) -- position
        **kwds -- standard 3d graphics options

    NOTE: There is no way to change the font size or opacity yet.

    EXAMPLES:
    We write the word SAGE in red at position (1,2,3):
        sage: text3d("SAGE", (1,2,3), color=(0.5,0,0))

    We draw a multicolor spiral of numbers:
        sage: sum([text3d('%.1f'%n, (cos(n),sin(n),n), color=(n/2,1-n/2,0)) \
                    for n in [0,0.2,..,8]])

    Another example
        sage: text3d("Sage is really neat!!",(2,12,1))

    And in 3d in two places:
        sage: text3d("Sage is...",(2,12,1), rgbcolor=(1,0,0)) + text3d("quite powerful!!",(4,10,0), rgbcolor=(0,0,1))


    """
    if not kwds.has_key('color') and not kwds.has_key('rgbcolor'):
        kwds['color'] = (0,0,0)
    G = Text(txt, **kwds).translate((x,y,z))
    G._set_extra_kwds(kwds)

    return G


class Point(PrimitiveObject):
    """
    Create a position in 3-space, represented by a sphere of fixed size.

    INPUT:
        center -- point (3-tuple)
        size -- (default: 1)
    """
    def __init__(self, center, size=1, **kwds):
        PrimitiveObject.__init__(self, **kwds)
        self.loc = (float(center[0]), float(center[1]), float(center[2]))
        self.size = size
        self._set_extra_kwds(kwds)

    def bounding_box(self):
        return self.loc, self.loc

    def tachyon_repr(self, render_params):
        transform = render_params.transform
        if transform is None:
            cen = self.loc
        else:
            cen = transform.transform_point(self.loc)
        return "Sphere center %s %s %s Rad %s %s" % (cen[0], cen[1], cen[2], self.size * TACHYON_PIXEL, self.texture.id)

    def obj_repr(self, render_params):
        T = render_params.transform
        if T is None:
            import transform
            T = transform.Transformation()
        render_params.push_transform(~T)
        S = shapes.Sphere(self.size / 200.0).translate(T(self.loc))
        cmds = S.obj_repr(render_params)
        render_params.pop_transform()
        return cmds

    def jmol_repr(self, render_params):
        name = render_params.unique_name('point')
        transform = render_params.transform
        cen = self.loc if transform is None else transform(self.loc)
        return ["draw %s DIAMETER %s {%s %s %s}\n%s" % (name, int(self.size), cen[0], cen[1], cen[2], self.texture.jmol_str('$' + name))]


class Line(PrimitiveObject):
    r"""
    Draw a 3d line joining a sequence of points.

    This line has a fixed diameter unaffected by transformations and zooming.
    It may be smoothed if \code{corner_cutoff < 1}.

    INPUT:
        points        -- list of points to pass through
        thickness     -- diameter of the line
        corner_cutoff -- threshold for smoothing (see the corners()
                         method) this is the minimum cosine between
                         adjacent segments to smooth
        arrow_head    -- if True make this curve into an arrow

    EXAMPLES:
        sage: from sage.plot.plot3d.shapes2 import Line
        sage: Line([(i*math.sin(i), i*math.cos(i), i/3) for i in range(30)], \
                   arrow_head=True)

    Smooth angles less than 90 degrees:
        sage: Line([(0,0,0),(1,0,0),(2,1,0),(0,1,0)], corner_cutoff=0)
    """
    def __init__(self, points, thickness=5, corner_cutoff=.5, arrow_head=False, **kwds):
        if len(points) < 2:
            raise ValueError, "there must be at least 2 points"
        PrimitiveObject.__init__(self, **kwds)
        self.points = points
        self.thickness = thickness
        self.corner_cutoff = corner_cutoff
        self.arrow_head = arrow_head

    def bounding_box(self):
        try:
            return self.__bounding_box
        except AttributeError:
            self.__bounding_box = point_list_bounding_box(self.points)
        return self.__bounding_box


    def tachyon_repr(self, render_params):
        T = render_params.transform
        cmds = []
        px, py, pz = self.points[0] if T is None else T(self.points[0])
        radius = self.thickness * TACHYON_PIXEL
        for P in self.points[1:]:
            x, y, z = P if T is None else T(P)
            if self.arrow_head and P is self.points[-1]:
                A = shapes.arrow3d((px, py, pz), (x, y, z), radius = radius, texture = self.texture)
                render_params.push_transform(~T)
                cmds.append(A.tachyon_repr(render_params))
                render_params.pop_transform()
            else:
                cmds.append("FCylinder base %s %s %s apex %s %s %s rad %s %s" % (px, py, pz,
                                                                                 x, y, z,
                                                                                 radius,
                                                                                 self.texture.id))
            px, py, pz = x, y, z
        return cmds

    def obj_repr(self, render_params):
        T = render_params.transform
        if T is None:
            import transform
            T = transform.Transformation()
        render_params.push_transform(~T)
        L = line3d([T(P) for P in self.points], radius=self.thickness / 200.0, arrow_head=self.arrow_head, texture=self.texture)
        cmds = L.obj_repr(render_params)
        render_params.pop_transform()
        return cmds

    def jmol_repr(self, render_params):
        T = render_params.transform
        corners = self.corners(max_len=255) # hardcoded limit in jmol
        last_corner = corners[-1]
        corners = set(corners)
        cmds = []
        cmd = None
        for P in self.points:
            TP = P if T is None else T(P)
            if P in corners:
                if cmd:
                    cmds.append(cmd + " {%s %s %s} " % TP)
                    cmds.append(self.texture.jmol_str('$'+name))
                type = 'arrow' if self.arrow_head and P is last_corner else 'curve'
                name = render_params.unique_name('line')
                cmd = "draw %s diameter %s %s {%s %s %s} " % (name, int(self.thickness), type, TP[0], TP[1], TP[2])
            else:
                cmd += " {%s %s %s} " % TP
        cmds.append(cmd)
        cmds.append(self.texture.jmol_str('$'+name))
        return cmds

    def corners(self, corner_cutoff=None, max_len=None):
        """
        Figures out where the curve turns too sharply to pretend it's smooth.

        INPUT:
            Maximum cosine of angle between adjacent line segments
            before adding a corner

        OUTPUT:
            List of points at which to start a new line.
            This always includes the first point, and never the last.

        EXAMPLES:

              Every point:
              sage: from sage.plot.plot3d.shapes2 import Line
              sage: Line([(0,0,0),(1,0,0),(2,1,0),(0,1,0)], corner_cutoff=1).corners()
              [(0, 0, 0), (1, 0, 0), (2, 1, 0)]

              Greater than 90 degrees:
              sage: Line([(0,0,0),(1,0,0),(2,1,0),(0,1,0)], corner_cutoff=0).corners()
              [(0, 0, 0), (2, 1, 0)]

              No corners:
              sage: Line([(0,0,0),(1,0,0),(2,1,0),(0,1,0)], corner_cutoff=-1).corners()
              (0, 0, 0)

              An intermediate value:
              sage: Line([(0,0,0),(1,0,0),(2,1,0),(0,1,0)], corner_cutoff=.5).corners()
              [(0, 0, 0), (2, 1, 0)]
        """
        if corner_cutoff is None:
            corner_cutoff = self.corner_cutoff
        if corner_cutoff >= 1:
            if max_len:
                self.points[:-1][::max_len-1]
            else:
                return self.points[:-1]
        elif corner_cutoff <= -1:
            return self.points[0]
        else:
            if not max_len:
                max_len = len(self.points)+1
            count = 2
            # ... -- prev -- cur -- next -- ...
            cur  = self.points[0]
            next = self.points[1]
            next_dir = [next[i] - cur[i] for i in range(3)]
            corners = [cur]
            cur, prev_dir = next, next_dir
            for next in self.points[2:]:
                if next == cur:
                    corners.append(cur)
                    cur = next
                    count = 1
                    continue
                next_dir = [next[i] - cur[i] for i in range(3)]
                cos_angle = dot(prev_dir, next_dir) / math.sqrt(dot(prev_dir, prev_dir) * dot(next_dir, next_dir))
                if cos_angle <= corner_cutoff or count > max_len-1:
                    corners.append(cur)
                    count = 1
                cur, prev_dir = next, next_dir
                count += 1
            return corners



def point3d(v, size=5, **kwds):
    """
    Plot a point or list of points in 3d space.

    INPUT:
        v -- a point or list of points
        size -- (default: 5) size of the point (or points)
        color -- a word that describes a color
        rgbcolor -- (r,g,b) with r, g, b between 0 and 1 that describes a color
        opacity -- (default: 1) if less than 1 then is transparent

    EXAMPLES:
        sage: sum([point3d((i,i^2,i^3), size=5) for i in range(10)])
    """
    if isinstance(v,(list,tuple)) and len(v) == 3 and not isinstance(v[0],(list,tuple)):
        return Point(v, size, **kwds)
    else:
        A = sum([Point(z, size, **kwds) for z in v])
        A._set_extra_kwds(kwds)
        return A


def dot((x0,y0,z0), (x1,y1,z1)):
    return x0*x1 + y0*y1 + z0*z1
