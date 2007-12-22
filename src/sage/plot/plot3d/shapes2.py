import shapes

def line3d(points, coerce=True, **kwds):
    v = []
    for i in range(len(points) - 1):
        v.append(shapes.LineSegment(points[i], points[i+1], **kwds))
    return sum(v)
