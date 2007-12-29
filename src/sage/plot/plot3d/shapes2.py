import shapes

def line3d(points, coerce=True, **kwds):
    v = []
    for i in range(len(points) - 1):
        v.append(shapes.LineSegment(points[i], points[i+1], **kwds))
    return sum(v)


def frame3d(size, **kwds):
    x,y,z = validate_frame_size(size)
    L1 = line3d([(x,y,z), (x,-y,z), (-x,-y,z), (-x,y,z),  (x,y,z), # top square
                 (x,y,-z), (x,-y,-z), (-x,-y,-z), (-x,y,-z),  (x,y,-z)],  # bottom square
                coerce=False, **kwds)
    # 3 additional lines joining top to bottom
    v2 = line3d([(x,-y,z), (x,-y,-z)], coerce=False, **kwds)
    v3 = line3d([(-x,y,z), (-x,y,-z)], coerce=False, **kwds)
    v4 = line3d([(-x,-y,z), (-x,-y,-z)], coerce=False, **kwds)
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


