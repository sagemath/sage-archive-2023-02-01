import math
import shapes

from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.misc.misc import srange

from texture import Texture

from shapes import Text, Sphere

def line3d(points, **kwds):
    """
    Draw a 3d line joining a sequence of points.
    """
    if len(points) < 2:
        raise ValueError, "there must be at least 2 points"
    v = []
    texture = Texture(kwds)
    for i in range(len(points) - 1):
        v.append(shapes.LineSegment(points[i], points[i+1], texture=texture, **kwds))
    w = sum(v)
    w._set_extra_kwds(kwds)
    return w


def frame3d(lower_left, upper_right, **kwds):
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
        n = max(0, 1 - nd(a/2.0))
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


def sphere((x,y,z)=(0,0,0), r=1, **kwds):
    """
    Return a plot of a sphere of radius $r$ centered at $(x,y,z)$.

    INPUT:
       (x,y,z) -- center (default: (0,0,0)
       r -- radius (default: 1)

    EXAMPLES:
    We draw a transparent sphere on a saddle.
       sage: u,v = var('u v')
       sage: sphere((0,0,1), color='red', opacity=0.5, aspect_ratio=[1,1,1]) + plot3d(u^2 - v^2, (u,-2,2), (v,-2,2))
    """
    G = Sphere(r, texture=Texture(kwds), **kwds)
    H = G.translate((x,y,z))
    H._set_extra_kwds(kwds)
    return H

def text3d(txt, (x,y,z), **kwds):
    """
    Display 3d text.

    INPUT:
        txt -- some text
        (x,y,z) -- position
        **kwds -- standard 3d graphics options

    This function called implicitly when you use the text command with a 3d position.

    NOTE: There is no way to change the font size or opacity yet.

    EXAMPLES:
    We write the word SAGE in red at position (1,2,3):
        sage: text("SAGE", (1,2,3), color=(0.5,0,0))

    We draw a multicolore spiral of numbers:
        sage: sum([text('%.1f'%n, (cos(n),sin(n),n), color=(n/2,1-n/2,0)) for n in [0,0.2,..,8]])
    """
    if not kwds.has_key('color') and not kwds.has_key('rgbcolor'):
        kwds['color'] = (0,0,0)
    G = Text(txt, **kwds).translate((x,y,z))
    G._set_extra_kwds(kwds)

    return G
