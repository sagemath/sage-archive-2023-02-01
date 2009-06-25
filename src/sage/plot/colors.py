#*****************************************************************************
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

from colorsys import hsv_to_rgb
from math import modf

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

    -  ``c`` - a string or 3-tuple

    OUTPUT: a 3-tuple of floats between 0 and 1.

    EXAMPLES::

        sage: from sage.plot.colors import to_mpl_color
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

    -  ``cmap`` - a colormap description (type ``import
        matplotlib.cm; matplotlib.cm.datad.keys()`` for valid names)

    EXAMPLES::

        sage: from sage.plot.colors import get_cmap
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


def hue(h, s=1, v=1):
    """
    hue(h,s=1,v=1) where ``h`` stands for hue, ``s`` stands for saturation,
    ``v`` stands for value. hue returns a tuple of rgb intensities (r, g,
    b).  All values are in the range 0 to 1.  Inputs ``h=0`` and ``h=1``
    yield red, while intermediate values yield orange, yellow, etc. with
    increasing input.

    The specifics of how hue, saturation, and value turn into an RGB
    color can be accessed from the source code of ``hsv_to_rgb`` by
    entering first ``from colorsys import hsv_to_rgb`` and then
    ``hsv_to_rgb??``.

    INPUT:

    -  ``h, s, v`` - real numbers between 0 and 1. Note that if any are
       not in this range they are automatically normalized to be in
       this range by reducing them modulo 1.

    OUTPUT: A valid RGB tuple.

    EXAMPLES::

        sage: hue(0.6)
        (0.0, 0.40000000000000036, 1.0)

    ``hue`` is an easy way of getting a broader range of colors for
    graphics::

        sage: plot(sin, -2, 2, rgbcolor=hue(0.6))
    """
    h = float(h); s = float(s); v = float(v)
    if h != 1:
        h = modf(h)[0]
        if h < 0:
            h += 1
    if s != 1:
        s = modf(s)[0]
        if s < 0:
            s += 1
    if v != 1:
        v = modf(v)
        if v < 0:
            v += 1
    c = hsv_to_rgb(h, s, v)
    return (float(c[0]),float(c[1]),float(c[2]))

def float_to_html(r, g, b):
    """
    This is a function to present tuples of RGB floats as HTML-happy
    hex for matplotlib. The input values should each be between 0 and 1.

    This may not seem necessary, but there are some  odd cases where
    matplotlib is just plain schizophrenic -- for an example, do

    EXAMPLES::

        sage: vertex_colors = {(1.0, 0.8571428571428571, 0.0): [4, 5, 6], (0.28571428571428559, 0.0, 1.0): [14, 15, 16], (1.0, 0.0, 0.0): [0, 1, 2, 3], (0.0, 0.57142857142857162, 1.0): [12, 13], (1.0, 0.0, 0.85714285714285676): [17, 18, 19], (0.0, 1.0, 0.57142857142857162): [10, 11], (0.28571428571428581, 1.0, 0.0): [7, 8, 9]}
        sage: graphs.DodecahedralGraph().show(vertex_colors=vertex_colors)

    Notice how the colors don't respect the partition at all.....

    ::

        sage: from sage.plot.colors import float_to_html
        sage: float_to_html(1.,1.,0.)
        '#ffff00'
        sage: float_to_html(.03,.06,.02)
        '#070f05'
    """ # TODO: figure out why this is necessary
    from sage.rings.integer import Integer
    from math import floor
    rr = Integer(int(floor(r*255))).str(base=16)
    gg = Integer(int(floor(g*255))).str(base=16)
    bb = Integer(int(floor(b*255))).str(base=16)
    rr = '0'*(2-len(rr)) + rr
    gg = '0'*(2-len(gg)) + gg
    bb = '0'*(2-len(bb)) + bb
    return '#' + rr\
               + gg\
               + bb

def rainbow(n, format='hex'):
    """
    Given an integer `n`, returns a list of colors, represented
    in HTML hex, that changes smoothly in hue from one end of the
    spectrum to the other. Written in order to easily represent vertex
    partitions on graphs, but useful elsewhere as well.

    AUTHORS:

    - Robert L. Miller

    - Karl-Dieter Crisman (directly use hsv_to_rgb for hues)

    INPUT:

    -  ``n`` - how many colors in the list of colors

    -  ``format`` - optional argument for what format to return the list of
       colors in (``hex`` or ``rgbtuple``)

    OUTPUT: A list of ``n`` colors, in either HTML hex or as RGB tuples.

    EXAMPLE::

        sage: from sage.plot.colors import rainbow
        sage: rainbow(7)
        ['#ff0000', '#ffda00', '#48ff00', '#00ff91', '#0091ff', '#4800ff', '#ff00da']
        sage: rainbow(7, 'rgbtuple')
        [(1, 0.0, 0.0), (1, 0.8571428571428571, 0.0), (0.28571428571428581, 1, 0.0), (0.0, 1, 0.57142857142857117), (0.0, 0.57142857142857162, 1), (0.2857142857142847, 0.0, 1), (1, 0.0, 0.85714285714285765)]
    """
    from sage.rings.integer import Integer
    n = Integer(n) # In case n is a Python int and i/n below would give 0!
    R = []
    for i in range(n):
        R.append(hsv_to_rgb(i/n,1,1))
    if format == 'rgbtuple':
        return R
    elif format == 'hex':
        for j in range(len(R)):
            R[j]=float_to_html(*R[j])
        return R

def rgbcolor(c):
    """
    Return the rgbcolor corresponding to c, which can be a Color, an
    RGB tuple, or a string.

    INPUT:

    -  ``c`` - a Color, 3-tuple, or string (HTML hex color).

    OUTPUT: rgb tuple of floats between 0 and 1.

    EXAMPLES::

        sage: from sage.plot.colors import rgbcolor
        sage: rgbcolor('purple')
        (0.5, 0.0, 1.0)
        sage: rgbcolor('#0033ea')
        (0.0, 0.19921875, 0.9140625)
        sage: rgbcolor((1,1/2,1/3))
        (1.0, 0.5, 0.33333333333333331)
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

        -  ``r,g,b`` - either a triple of floats between 0 and 1,
           OR ``r`` - a color string or HTML color hex string

        EXAMPLES::

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

        EXAMPLES::

            sage: Color('#8000ff').__repr__()
            'RGB color (0.5, 0.0, 0.99609375)'
        """
        return "RGB color %s"%(self.__rgb,)

    def rgb(self):
        """
        Return underlying RGB tuple.

        OUTPUT: 3-tuple

        EXAMPLES::

            sage: Color('#8000ff').rgb()
            (0.5, 0.0, 0.99609375)
        """
        return self.__rgb

    def html_color(self):
        """
        Return color formatted as an HTML hex color.

        OUTPUT: string of length 7.

        EXAMPLES::

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
