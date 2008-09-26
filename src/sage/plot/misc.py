from functools import wraps
from copy import copy


def ensure_subs(f):
    if not hasattr(f, 'subs'):
        from sage.calculus.all import SR
        return SR(f)
    return f


colors = {
    "red"   : (1.0,0.0,0.0),
    "orange": (1.0,.5,0.0),
    "yellow": (1.0,1.0,0.0),
    "green" : (0,1.0,0.0),
    "blue"  : (0.0,0.0,1.0),
    "purple": (.5,0.0,1.0),
    "white" : (1.0,1.0,1.0),
    "black" : (0.0,0.0,0.0),
    'brown': (0.65, 0.165, 0.165),
    "grey"  : (.5,.5,.5),
    "gray"  : (.5,.5,.5),
    "lightblue" : (0.4,0.4,1.0),
    "automatic": (0.4,0.4,1.0)
}

def rgbcolor(c):
    """
    Return the rgbcolor corresponding to c, which can be a Color, tuple, or string.

    INPUT:
        c -- a Color, 3-tuple, or string (HTML hex color).

    OUTPUT:
        rgb tuple of floats between 0 and 1.
    """
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
        """
        A color object.

        INPUT:
            r,g,b -- either a triple of floats between 0 and 1, OR
            r -- a color string or HTML color hex string

        EXAMPLES:
            sage: Color('purple')
            RGB color (0.5, 0.0, 1.0)
            sage: Color(0.5,0,1)
            RGB color (0.5, 0.0, 1.0)
            sage: Color('#8000ff')
            RGB color (0.5, 0.0, 0.99609375)
        """
        if g is None and b is None:
            self.__rgb = rgbcolor(r)
        else:
            self.__rgb = (float(r),float(g),float(b))

    def __repr__(self):
        """
        Return string representation of this RGB color.

        EXAMPLES:
            sage: Color('#8000ff').__repr__()
            'RGB color (0.5, 0.0, 0.99609375)'
        """
        return "RGB color %s"%(self.__rgb,)

    def rgb(self):
        """
        Return underlying RGB tuple.

        OUTPUT:
            3-tuple

        EXAMPLES:
            sage: Color('#8000ff').rgb()
            (0.5, 0.0, 0.99609375)
        """
        return self.__rgb

    def html_color(self):
        """
        Return color formated as an HTML hex color.

        OUTPUT:
            string of length 7.

        EXAMPLES:
            sage: Color('yellow').html_color()
            '#ffff00'
        """
        s = '#'
        for z in self.__rgb:
            h = '%x'%int(z*256)
            if len(h) > 2:
                h = 'ff'
            elif len(h) == 1:
                h = '0' + h
            s += h
        return s

class options(object):
    def __init__(self, **options):
        """
        A decorator for functions which allows for default options to be
        set and reset by the end user.  Additionally, if one needs to, one
        can get at the original keyword arguments passed into the
        decorator.

        TESTS:
            sage: from sage.plot.misc import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: o.options
            {'rgbcolor': (0, 0, 1)}
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: o.original_opts
            True
            sage: loads(dumps(o)).options
            {'rgbcolor': (0, 0, 1)}
        """
        self.options = options
        self.original_opts = options.pop('__original_opts', False)

    def __call__(self, func):
        """
        EXAMPLES:
            sage: from sage.plot.misc import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: def f(*args, **kwds): print args, list(sorted(kwds.items()))
            sage: f1 = o(f)
            sage: f1()
            () [('rgbcolor', (0, 0, 1))]
            sage: f1(rgbcolor=1)
            () [('rgbcolor', 1)]
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: f2 = o(f)
            sage: f2(alpha=1)
            () [('__original_opts', {'alpha': 1}), ('alpha', 1), ('rgbcolor', (0, 0, 1))]

        """
        @wraps(func)
        def wrapper(*args, **kwds):
            options = copy(wrapper.options)
            if self.original_opts:
                options['__original_opts'] = kwds
            options.update(kwds)
            return func(*args, **options)


        def reset():
            """
            Reset the options to the defaults.

            EXAMPLES:
                sage: from sage.plot.misc import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds): print args, list(sorted(kwds.items()))
                sage: f = o(f)
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.options
                {'rgbcolor': (1, 1, 1)}
                sage: f()
                () [('rgbcolor', (1, 1, 1))]
                sage: f.reset()
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f()
                () [('rgbcolor', (0, 0, 1))]
            """
            wrapper.options = copy(self.options)

        wrapper.options = copy(self.options)
        wrapper.reset = reset

        return wrapper

class rename_keyword(object):
    def __init__(self, **renames):
        """
        A decorator which renames keyword arguments.

        EXAMPLES:
            sage: from sage.plot.misc import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: r.renames
            {'color': 'rgbcolor'}
            sage: loads(dumps(r)).renames
            {'color': 'rgbcolor'}

        """
        self.renames = renames

    def __call__(self, func):
        """
        EXAMPLES:
            sage: from sage.plot.misc import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: def f(*args, **kwds): print args, kwds
            sage: f = r(f)
            sage: f()
            () {}
            sage: f(alpha=1)
            () {'alpha': 1}
            sage: f(rgbcolor=1)
            () {'rgbcolor': 1}
            sage: f(color=1)
            () {'rgbcolor': 1}
        """
        @wraps(func)
        def wrapper(*args, **kwds):
            for old_name, new_name in self.renames.items():
                if kwds.has_key(old_name) and not kwds.has_key(new_name):
                    kwds[new_name] = kwds[old_name]
                    del kwds[old_name]
            return func(*args, **kwds)
        return wrapper
