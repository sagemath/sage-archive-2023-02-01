
def ensure_subs(f):
    if not hasattr(f, 'subs'):
        from sage.calculus.all import SR
        return SR(f)
    return f


colors = {
    "red"   : (1,0,0),
    "orange": (1,.5,0),
    "yellow": (1,1,0),
    "green" : (0,1,0),
    "blue"  : (0,0,1),
    "purple": (.5,0,1),
    "white" : (1,1,1),
    "black" : (0,0,0),
    'brown': (0.65, 0.165, 0.165),
    "grey"  : (.5,.5,.5),
    "gray"  : (.5,.5,.5),
    "lightblue" : (0.4,0.4,1),
    "automatic": (0.4,0.4,1)
}

def rgbcolor(c):
    if isinstance(c, tuple):
        return (float(c[0]), float(c[1]), float(c[2]))
    if isinstance(c, str):
        try:
            return colors[c]
        except KeyError:
            pass

    raise ValueError, "unknown color '%s'"%c

