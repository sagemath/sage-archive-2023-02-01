from functools import wraps
from copy import copy
from math import modf

from sage.misc.misc import verbose

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

def to_mpl_color(c):
    """
    Convert a tuple or string to a matplotlib rgb color tuple.

    INPUT:
        c -- string or 3-tuple

    OUTPUT:
        3-tuple of floats between 0 and 1.

    EXAMPLES:
        sage: from sage.plot.plot import to_mpl_color
        sage: to_mpl_color('#fa0')
        (1.0, 0.66666666666666663, 0.0)
        sage: to_mpl_color('#ffffe1')
        (1.0, 1.0, 0.88235294117647056)
        sage: to_mpl_color('blue')
        (0.0, 0.0, 1.0)
        sage: to_mpl_color([1,1/2,1/3])
        (1.0, 0.5, 0.33333333333333331)
        sage: to_mpl_color([1,2,255])   # WARNING -- numbers are reduced mod 1!!
        (1.0, 0.0, 0.0)
    """
    if isinstance(c, Color):
        c = c.rgb()

    if isinstance(c, str):
        if len(c) > 0 and c[0] == '#':
            # it is some sort of html like color, e.g, #00ffff or #ab0
            h = c[1:]
            if len(h) == 3:
                h = '%s%s%s%s%s%s'%(h[0],h[0], h[1],h[1], h[2],h[2])
            elif len(h) != 6:
                raise ValueError, "color hex string (= '%s') must have length 3 or 6"%h
            return tuple([eval('0x%s'%h[i:i+2])/float(255) for i in [0,2,4]])
        else:
            try:
                return colors[c]
            except KeyError:
                raise ValueError, "unknown color '%s'"%c

    elif isinstance(c, (list, tuple)):
        c = list(c)
        if len(c) != 3:
            raise ValueError, "color tuple must have 3 entries, one for each RGB channel"
        for i in range(len(c)):
            s = float(c[i])
            if s != 1:
                s = modf(s)[0]
                if s < 0:
                    s += 1
            c[i] = s

    else:
        raise TypeError, "c must be a list, tuple, or string"

    return tuple(c)

# If you change what this accepts, remember to change the documentation
# about cmap where it is used and to test these classes.
def get_cmap(cmap):
    r"""
    Returns the colormap corresponding to cmap.

    INPUT:
        cmap -- a colormap description (type cmap_help() for more information)

    EXAMPLES:
        sage: from sage.plot.misc import get_cmap
        sage: get_cmap('jet')
        <matplotlib.colors.LinearSegmentedColormap instance at 0x...>
        sage: get_cmap([(0,0,0), (0.5,0.5,0.5), (1,1,1)])
        <matplotlib.colors.ListedColormap instance at 0x...>
        sage: get_cmap(['green', 'lightblue', 'blue'])
        <matplotlib.colors.ListedColormap instance at 0x...>
        sage: get_cmap('jolies')
        Traceback (most recent call last):
        ...
        RuntimeError: Color map jolies not known (type import matplotlib.cm; matplotlib.cm.datad.keys() for valid names)
        sage: get_cmap('mpl')
        Traceback (most recent call last):
        ...
        RuntimeError: Color map mpl not known (type import matplotlib.cm; matplotlib.cm.datad.keys() for valid names)
    """
    #cm is the matplotlib color map module
    from matplotlib import cm
    from matplotlib.colors import ListedColormap, Colormap
    if isinstance(cmap, Colormap):
        return cmap
    elif isinstance(cmap, str):
        if not cmap in cm.datad.keys():
            raise RuntimeError, "Color map %s not known (type import matplotlib.cm; matplotlib.cm.datad.keys() for valid names)"%cmap
        return cm.__dict__[cmap]
    elif isinstance(cmap, (list, tuple)):
        cmap = map(rgbcolor, cmap)
        return ListedColormap(cmap)

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

class suboptions(object):
    def __init__(self, name, **options):
        """
        A decorator for functions which collects all keywords
        starting with name_ and collects them into a dictionary
        which will be passed on to the wrapped function as a
        dictionary called name_options.

        The keyword arguments passed into the constructor are taken
        to be default for the name_options dictionary.

        EXAMPLES:
            sage: from sage.plot.misc import suboptions
            sage: s = suboptions('arrow', size=2)
            sage: s.name
            'arrow_'
            sage: s.options
            {'size': 2}
        """
        self.name = name + "_"
        self.options = options

    def __call__(self, func):
        """
        Returns a wrapper around func

        EXAMPLES:
            sage: from sage.plot.misc import suboptions
            sage: def f(*args, **kwds): print list(sorted(kwds.items()))
            sage: f = suboptions('arrow', size=2)(f)
            sage: f(size=2)
            [('arrow_options', {'size': 2}), ('size', 2)]
            sage: f(arrow_size=3)
            [('arrow_options', {'size': 3})]
            sage: f(arrow_options={'size':4})
            [('arrow_options', {'size': 4})]
            sage: f(arrow_options={'size':4}, arrow_size=5)
            [('arrow_options', {'size': 5})]

        """
        @wraps(func)
        def wrapper(*args, **kwds):
            suboptions = copy(self.options)
            suboptions.update(kwds.pop(self.name+"options", {}))

            #Collect all the relevant keywords in kwds
            #and put them in suboptions
            for key, value in kwds.items():
                if key.startswith(self.name):
                    suboptions[key[len(self.name):]] = value
                    del kwds[key]

            kwds[self.name + "options"] = suboptions

            return func(*args, **kwds)

        from sage.misc.sageinspect import sage_getsource
        wrapper._sage_src_ = lambda: sage_getsource(func)

        return wrapper

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


        def defaults():
            """
            Return the default options.

            EXAMPLES:
                sage: from sage.plot.misc import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds): print args, list(sorted(kwds.items()))
                sage: f = o(f)
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.defaults()
                {'rgbcolor': (0, 0, 1)}
            """
            return copy(self.options)

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
        wrapper.reset.__doc__ = """
        Reset the options to the defaults.

        Defaults:
        %s
        """%self.options

        wrapper.defaults = defaults
        wrapper.defaults.__doc__ = """
        Return the default options.

        Defaults:
        %s
        """%self.options


        from sage.misc.sageinspect import sage_getsource
        wrapper._sage_src_ = lambda: sage_getsource(func)

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

        from sage.misc.sageinspect import sage_getsource
        wrapper._sage_src_ = lambda: sage_getsource(func)

        return wrapper
