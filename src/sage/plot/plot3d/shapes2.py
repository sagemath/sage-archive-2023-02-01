import math
import shapes

from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.misc.misc import srange

def line3d(points, coerce=True, color="lightblue", **kwds):
    if len(points) < 2:
        raise ValueError, "there must be at least 2 points"
    v = []
    for i in range(len(points) - 1):
        v.append(shapes.LineSegment(points[i], points[i+1], color=color, **kwds))
    return sum(v)

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
    return L1 + v2 + v3 + v4

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

