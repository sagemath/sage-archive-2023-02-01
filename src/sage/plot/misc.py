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
    if isinstance(c, Color):
        return c.rgb()
    if isinstance(c, tuple):
        return (float(c[0]), float(c[1]), float(c[2]))
    if isinstance(c, str):
        if len(c) == 7 and c[0] == '#':  # html hex color
            # we use Integer instead of 0x eval for security reasons
            return tuple([int(c[i:i+2], base=16)/float(256) for i in [1,3,5]])
        try:
            return colors[c]
        except KeyError:
            pass

    raise ValueError, "unknown color '%s'"%c


class Color:
    def __init__(self, r='#0000ff', g=None, b=None):
        if g is None and b is None:
            self.__rgb = rgbcolor(r)
        else:
            self.__rgb = (float(r),float(g),float(b))

    def __repr__(self):
        return "RGB color %s"%(self.__rgb,)

    def rgb(self):
        return self.__rgb

    def html_color(self):
        s = '#'
        for z in self.__rgb:
            h = '%x'%int(z*256)
            if len(h) > 2:
                h = 'ff'
            elif len(h) == 1:
                h = '0' + h
            s += h
        return s
