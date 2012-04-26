r"""
Graphics objects

This file contains the definition of the classes :class:`Graphics` and
:class:`GraphicsArray`.  Usually, you don't create these classes directly
(although you can do it), you would use :func:`plot` or
:func:`graphics_array` instead.

AUTHORS:

- Jeroen Demeyer (2012-04-19): split off this file from plot.py (:trac:`12857`)

"""

#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>
#       Copyright (C) 2006-2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Jason Grout
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import sage.misc.misc
from sage.misc.html import html
from sage.structure.sage_object import SageObject
from sage.misc.decorators import suboptions
from colors import rgbcolor

ALLOWED_EXTENSIONS = ['.eps', '.pdf', '.png', '.ps', '.sobj', '.svg']
DEFAULT_DPI = 100
DOCTEST_MODE_FILE = os.path.join(sage.misc.misc.SAGE_TMP, 'test.png')
SHOW_DEFAULT = True

def show_default(default=None):
    r"""
    Set the default for showing plots using any plot commands. If
    called with no arguments, returns the current default.

    If this is ``True`` (the default) then any plot object
    when displayed will be displayed as an actual plot instead of text,
    i.e., the show command is not needed.

    EXAMPLES: The default starts out as ``True``::

        sage: show_default()
        True

    We set it to ``False``.

    ::

        sage: show_default(False)

    We see that it is ``False``.

    ::

        sage: show_default()
        False

    Now plot commands will not display their plots by default.

    Turn back on default display.

    ::

        sage: show_default(True)
    """
    global SHOW_DEFAULT
    if default is None:
        return SHOW_DEFAULT
    SHOW_DEFAULT = bool(default)

# If do_verify is True, options are checked when drawing a
# GraphicsPrimitive.  See primitive.py
do_verify = True

def is_Graphics(x):
    """
    Return True if `x` is a Graphics object.

    EXAMPLES::

        sage: from sage.plot.graphics import is_Graphics
        sage: is_Graphics(1)
        False
        sage: is_Graphics(disk((0.0, 0.0), 1, (0, pi/2)))
        True
    """
    return isinstance(x, Graphics)

class Graphics(SageObject):
    """
    The Graphics object is an empty list of graphics objects It is
    useful to use this object when initializing a for loop where
    different graphics object will be added to the empty object.

    EXAMPLES::

        sage: G = Graphics(); print G
        Graphics object consisting of 0 graphics primitives
        sage: c = circle((1,1), 1)
        sage: G+=c; print G
        Graphics object consisting of 1 graphics primitive

    Here we make a graphic of embedded isosceles triangles, coloring
    each one with a different color as we go::

        sage: h=10; c=0.4; p=0.5;
        sage: G = Graphics()
        sage: for x in srange(1,h+1):
        ...        l = [[0,x*sqrt(3)],[-x/2,-x*sqrt(3)/2],[x/2,-x*sqrt(3)/2],[0,x*sqrt(3)]]
        ...        G+=line(l,color=hue(c + p*(x/h)))
        sage: G.show(figsize=[5,5])

    TESTS:

    From :trac:`4604`, ensure Graphics can handle 3d objects::

        sage: g = Graphics()
        sage: g += sphere((1, 1, 1), 2)
        sage: g.show()

    We check that graphics can be pickled (we can't use equality on
    graphics so we just check that the load/dump cycle gives a
    :class:`Graphics` instance)::

        sage: g = Graphics()
        sage: g2 = loads(dumps(g))
        sage: g2.show()

    ::

        sage: isinstance(g2, Graphics)
        True
    """

    def __init__(self):
        """
        Create a new empty Graphics objects with all the defaults.

        EXAMPLES::

            sage: G = Graphics()
        """
        self.__fontsize = 10
        self.__show_axes = True
        self.__show_legend = False
        self.__legend_opts = {}
        self.__axes_color = (0, 0, 0)
        self.__axes_label_color = (0, 0, 0)
        self.__tick_label_color = (0, 0, 0)
        self.__axes_width = 0.8
        self.__objects = []
        self._extra_kwds = {}
        self.__bbox_extra_artists = []

    def set_aspect_ratio(self, ratio):
        """
        Set the aspect ratio, which is the ratio of height and width
        of a unit square (i.e., height/width of a unit square), or
        'automatic' (expand to fill the figure).

        INPUT:


        -  ``ratio`` - a positive real number or 'automatic'


        EXAMPLES: We create a plot of the upper half of a circle, but it
        doesn't look round because the aspect ratio is off::

            sage: P = plot(sqrt(1-x^2),(x,-1,1)); P

        So we set the aspect ratio and now it is round::

            sage: P.set_aspect_ratio(1)
            sage: P.aspect_ratio()
            1.0
            sage: P

        Note that the aspect ratio is inherited upon addition (which takes
        the max of aspect ratios of objects whose aspect ratio has been
        set)::

            sage: P + plot(sqrt(4-x^2),(x,-2,2))

        In the following example, both plots produce a circle that looks
        twice as tall as wide::

            sage: Q = circle((0,0), 0.5); Q.set_aspect_ratio(2)
            sage: (P + Q).aspect_ratio(); P+Q
            2.0
            sage: (Q + P).aspect_ratio(); Q+P
            2.0
        """
        if ratio != 'auto' and ratio != 'automatic':
            ratio = float(ratio)
            if ratio <= 0:
                raise ValueError, "the aspect ratio must be positive or 'automatic'"
        else:
            ratio = 'automatic'
        self._extra_kwds['aspect_ratio'] = ratio

    def aspect_ratio(self):
        """
        Get the current aspect ratio, which is the ratio of height to
        width of a unit square, or 'automatic'.

        OUTPUT: a positive float (height/width of a unit square), or 'automatic'
        (expand to fill the figure).

        EXAMPLES:

        The default aspect ratio for a new blank Graphics object is 'automatic'::

            sage: P = Graphics()
            sage: P.aspect_ratio()
            'automatic'

        The aspect ratio can be explicitly set different than the object's default::

            sage: P = circle((1,1), 1)
            sage: P.aspect_ratio()
            1.0
            sage: P.set_aspect_ratio(2)
            sage: P.aspect_ratio()
            2.0
            sage: P.set_aspect_ratio('automatic')
            sage: P.aspect_ratio()
            'automatic'
        """
        return self._extra_kwds.get('aspect_ratio', 'automatic')

    def legend(self, show=None):
        r"""
        Set whether or not the legend is shown by default.

        INPUT:

        -  ``show`` - (default: None) a boolean

        If called with no input, return the current legend setting.

        EXAMPLES:

        By default no legend is displayed::

            sage: P = plot(sin)
            sage: P.legend()
            False

        But if we put a label then the legend is shown::

            sage: P = plot(sin, legend_label='sin')
            sage: P.legend()
            True

        We can turn it on or off::

            sage: P.legend(False)
            sage: P.legend()
            False
            sage: P.legend(True)
            sage: P # show with the legend
        """
        if show is None:
            return self.__show_legend
        else:
            self.__show_legend = bool(show)

    def set_legend_options(self, **kwds):
        r"""
        Set various legend options.

        INPUT:

        - ``title`` - (default: None) string, the legend title

        - ``ncol`` - (default: 1) positive integer, the number of columns

        - ``columnspacing`` - (default: None) the spacing between columns

        - ``borderaxespad`` - (default: None) float, length between the axes and the legend

        - ``back_color`` - (default: (0.9, 0.9, 0.9)) This parameter can be a string
          denoting a color or an RGB tuple. The string can be a color name
          as in ('red', 'green', 'yellow', ...) or a floating point number
          like '0.8' which gets expanded to (0.8, 0.8, 0.8). The
          tuple form is just a floating point RGB tuple with all values ranging
          from 0 to 1.

        - ``handlelength`` - (default: 0.05) float, the length of the legend handles

        - ``handletextpad`` - (default: 0.5) float, the pad between the legend handle and text

        - ``labelspacing`` - (default: 0.02) float, vertical space between legend entries

        - ``loc`` - (default: 'best') May be a string, an integer or a tuple. String or
              integer inputs must be one of the following:

          - 0, 'best'

          - 1, 'upper right'

          - 2, 'upper left'

          - 3, 'lower left'

          - 4, 'lower right'

          - 5, 'right'

          - 6, 'center left'

          - 7, 'center right'

          - 8, 'lower center'

          - 9, 'upper center'

          - 10, 'center'

          - Tuple arguments represent an absolute (x, y) position on the plot
            in axes coordinates (meaning from 0 to 1 in each direction).

        - ``markerscale`` - (default: 0.6) float, how much to scale the markers in the legend.

        - ``numpoints`` - (default: 2) integer, the number of points in the legend for line

        - ``borderpad`` - (default: 0.6) float, the fractional whitespace inside the legend border
          (between 0 and 1)

        - ``font_family`` - (default: 'sans-serif') string, one of 'serif', 'sans-serif',
          'cursive', 'fantasy', 'monospace'

        - ``font_style`` - (default: 'normal') string, one of 'normal', 'italic', 'oblique'

        - ``font_variant`` - (default: 'normal') string, one of 'normal', 'small-caps'

        - ``font_weight`` - (default: 'medium') string, one of 'black', 'extra bold', 'bold',
          'semibold', 'medium', 'normal', 'light'

        - ``font_size`` - (default: 'medium') string, one of 'xx-small', 'x-small', 'small',
          'medium', 'large', 'x-large', 'xx-large' or an absolute font size (e.g. 12)

        -  ``shadow`` - (default: False) boolean - draw a shadow behind the legend

        - ``fancybox`` - (default: False) a boolean.  If True, draws a frame with a round
          fancybox.

        These are all keyword arguments.

        OUTPUT: a dictionary of all current legend options

        EXAMPLES:

        By default, no options are set::

            sage: p = plot(tan, legend_label='tan')
            sage: p.set_legend_options()
            {}

        We build a legend with a shadow::

            sage: p.set_legend_options(shadow=True)
            sage: p.set_legend_options()['shadow']
            True

        To set the legend position to the center of the plot, all these
        methods are roughly equivalent::

            sage: p.set_legend_options(loc='center'); p

        ::

            sage: p.set_legend_options(loc=10); p

        ::

            sage: p.set_legend_options(loc=(0.5,0.5)); p # aligns the bottom of the box to the center
        """
        if len(kwds) == 0:
            return self.__legend_opts
        else:
            self.__legend_opts.update(kwds)


    def get_axes_range(self):
        """
        Returns a dictionary of the range of the axes for this graphics
        object.  This is fall back to the ranges in get_minmax_data() for
        any value which the user has not explicitly set.

        .. warning::

           Changing the dictionary returned by this function does not
           change the axes range for this object.  To do that, use the
           :meth:`set_axes_range` method.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: list(sorted(L.get_axes_range().items()))
            [('xmax', 3.0), ('xmin', 1.0), ('ymax', 5.0), ('ymin', -4.0)]
            sage: L.set_axes_range(xmin=-1)
            sage: list(sorted(L.get_axes_range().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 5.0), ('ymin', -4.0)]
        """
        axes_range = self.get_minmax_data()
        axes_range.update(self._get_axes_range_dict())
        return axes_range

    def set_axes_range(self, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Set the ranges of the `x` and `y` axes.

        INPUT:


        -  ``xmin, xmax, ymin, ymax`` - floats


        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.set_axes_range(-1, 20, 0, 2)
            sage: d = L.get_axes_range()
            sage: d['xmin'], d['xmax'], d['ymin'], d['ymax']
            (-1.0, 20.0, 0.0, 2.0)
        """
        l = locals()
        axes_range = self._get_axes_range_dict()
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            if l[name] is not None:
                axes_range[name] = float(l[name])

    axes_range = set_axes_range

    def _get_axes_range_dict(self):
        """
        Returns the underlying dictionary used to store the user's
        custom ranges for the axes on this object.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L._get_axes_range_dict()
            {}
            sage: L.set_axes_range(xmin=-1)
            sage: L._get_axes_range_dict()
            {'xmin': -1.0}
        """
        try:
            return self.__axes_range
        except AttributeError:
            self.__axes_range = {}
            return self.__axes_range

    def fontsize(self, s=None):
        """
        Set the font size of axes labels and tick marks.

        INPUT:


        -  ``s`` - integer, a font size in points.


        If called with no input, return the current fontsize.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.fontsize()
            10
            sage: L.fontsize(20)
            sage: L.fontsize()
            20

        All the numbers on the axes will be very large in this plot::

            sage: L
        """
        if s is None:
            try:
                return self.__fontsize
            except AttributeError:
                self.__fontsize = 10
                return self.__fontsize
        self.__fontsize = int(s)

    def axes(self, show=None):
        """
        Set whether or not the `x` and `y` axes are shown
        by default.

        INPUT:


        -  ``show`` - bool


        If called with no input, return the current axes setting.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])

        By default the axes are displayed.

        ::

            sage: L.axes()
            True

        But we turn them off, and verify that they are off

        ::

            sage: L.axes(False)
            sage: L.axes()
            False

        Displaying L now shows a triangle but no axes.

        ::

            sage: L
        """
        if show is None:
            try:
                return self.__show_axes
            except AttributeError:
                self.__show_axes = True
                return self.__show_axes
        self.__show_axes = bool(show)

    def axes_color(self, c=None):
        """
        Set the axes color.

        If called with no input, return the current axes_color setting.

        INPUT:


        -  ``c`` - an RGB color 3-tuple, where each tuple entry
           is a float between 0 and 1


        EXAMPLES: We create a line, which has like everything a default
        axes color of black.

        ::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.axes_color()
            (0, 0, 0)

        We change the axes color to red and verify the change.

        ::

            sage: L.axes_color((1,0,0))
            sage: L.axes_color()
            (1.0, 0.0, 0.0)

        When we display the plot, we'll see a blue triangle and bright red
        axes.

        ::

            sage: L
        """
        if c is None:
            try:
                return self.__axes_color

            except AttributeError:
                self.__axes_color = (0.0, 0.0, 0.0)
                return self.__axes_color
        self.__axes_color = rgbcolor(c)

    def axes_labels(self, l=None):
        """
        Set the axes labels.

        INPUT:


        -  ``l`` - (default: None) a list of two strings or
           None


        OUTPUT: a 2-tuple of strings

        If l is None, returns the current ``axes_labels``,
        which is itself by default None. The default labels are both
        empty.

        EXAMPLES: We create a plot and put x and y axes labels on it.

        ::

            sage: p = plot(sin(x), (x, 0, 10))
            sage: p.axes_labels(['$x$','$y$'])
            sage: p.axes_labels()
            ('$x$', '$y$')

        Now when you plot p, you see x and y axes labels::

            sage: p

        Notice that some may prefer axes labels which are not
        typeset::

            sage: plot(sin(x), (x, 0, 10), axes_labels=['x','y'])
        """
        if l is None:
            try:
                return self.__axes_labels
            except AttributeError:
                self.__axes_labels = None
                return self.__axes_labels
        if not isinstance(l, (list, tuple)):
            raise TypeError, "l must be a list or tuple"
        if len(l) != 2:
            raise ValueError, "l must have length 2"
        self.__axes_labels = (str(l[0]), str(l[1]))

    def axes_label_color(self, c=None):
        r"""
        Set the color of the axes labels.

        The axes labels are placed at the edge of the x and y axes, and are
        not on by default (use the ``axes_labels`` command to
        set them; see the example below). This function just changes their
        color.

        INPUT:


        -  ``c`` - an RGB 3-tuple of numbers between 0 and 1


        If called with no input, return the current axes_label_color
        setting.

        EXAMPLES: We create a plot, which by default has axes label color
        black.

        ::

            sage: p = plot(sin, (-1,1))
            sage: p.axes_label_color()
            (0, 0, 0)

        We change the labels to be red, and confirm this::

            sage: p.axes_label_color((1,0,0))
            sage: p.axes_label_color()
            (1.0, 0.0, 0.0)

        We set labels, since otherwise we won't see anything.

        ::

            sage: p.axes_labels(['$x$ axis', '$y$ axis'])

        In the plot below, notice that the labels are red::

            sage: p
        """
        if c is None:
            try:
                return self.__axes_label_color
            except AttributeError:
                self.__axes_label_color = (0, 0, 0)
                return self.__axes_label_color
        self.__axes_label_color = rgbcolor(c)


    def axes_width(self, w=None):
        r"""
        Set the axes width. Use this to draw a plot with really fat or
        really thin axes.

        INPUT:


        -  ``w`` - a float


        If called with no input, return the current
        ``axes_width`` setting.

        EXAMPLE: We create a plot, see the default axes width (with funny
        Python float rounding), then reset the width to 10 (very fat).

        ::

            sage: p = plot(cos, (-3,3))
            sage: p.axes_width()
            0.8
            sage: p.axes_width(10)
            sage: p.axes_width()
            10.0

        Finally we plot the result, which is a graph with very fat axes.

        ::

            sage: p
        """
        if w is None:
            try:
                return self.__axes_width
            except AttributeError:
                self.__axes_width = True
                return self.__axes_width
        self.__axes_width = float(w)

    def tick_label_color(self, c=None):
        """
        Set the color of the axes tick labels.

        INPUT:


        -  ``c`` - an RGB 3-tuple of numbers between 0 and 1


        If called with no input, return the current tick_label_color
        setting.

        EXAMPLES::

            sage: p = plot(cos, (-3,3))
            sage: p.tick_label_color()
            (0, 0, 0)
            sage: p.tick_label_color((1,0,0))
            sage: p.tick_label_color()
            (1.0, 0.0, 0.0)
            sage: p
        """
        if c is None:
            try:
                return self.__tick_label_color
            except AttributeError:
                self.__tick_label_color = (0, 0, 0)
                return self.__tick_label_color
        self.__tick_label_color = rgbcolor(c)

    def _repr_(self):
        r"""
        Show this graphics objects.

        If the ``show_default`` function has been called with
        True (the default), then you'll see this graphics object displayed.
        Otherwise you'll see a text representation of it.

        EXAMPLES: We create a plot and call ``_repr_`` on it,
        which causes it to be displayed as a plot::

            sage: P = plot(cos, (-1,1))
            sage: P._repr_()
            ''

        Just doing this also displays the plot::

            sage: P

        Note that printing P with the ``print`` statement does
        not display the plot::

            sage: print P
            Graphics object consisting of 1 graphics primitive

        Now we turn off showing plots by default::

            sage: show_default(False)

        Now we just get a string. To show P you would have to do
        ``show(P)``.

        ::

            sage: P._repr_()
            'Graphics object consisting of 1 graphics primitive'
            sage: P
            Graphics object consisting of 1 graphics primitive

        Finally, we turn ``show_default`` back on::

            sage: show_default(True)
        """
        if SHOW_DEFAULT:
            self.show()
            return ''
        else:
            return self.__str__()

    def __str__(self):
        r"""
        Return string representation of this plot.

        EXAMPLES::

            sage: S = circle((0,0), 2); S.__str__()
            'Graphics object consisting of 1 graphics primitive'
            sage: print S
            Graphics object consisting of 1 graphics primitive

        .. warning::

           ``__str__`` is not called when printing lists of graphics
           objects, which can be confusing, since they will all pop
           up. One workaround is to call ``show_default``:

        For example, below when we do ``print v`` two plots are
        displayed::

            sage: v = [circle((0,0), 2), circle((2,3), 1)]
            sage: print v
            [, ]

        However, if we call ``show_default`` then we see the
        text representations of the graphics::

            sage: show_default(False)
            sage: print v
            [Graphics object consisting of 1 graphics primitive, Graphics object consisting of 1 graphics primitive]
            sage: v
            [Graphics object consisting of 1 graphics primitive,
             Graphics object consisting of 1 graphics primitive]

        ::

            sage: show_default(True)
        """
        pr, i = '', 0
        for x in self:
            pr += '\n\t%s -- %s'%(i, x)
            i += 1
        s = "Graphics object consisting of %s graphics primitives"%(len(self))
        if len(self) == 1:
            s = s[:-1]
        return s

    def __getitem__(self, i):
        """
        Returns the ith graphics primitive object:

        EXAMPLE::

            sage: G = circle((1,1),2) + circle((2,2),5); print G
            Graphics object consisting of 2 graphics primitives
            sage: G[1]
            Circle defined by (2.0,2.0) with r=5.0
        """
        return self.__objects[i]

    def __len__(self):
        """
        If G is of type Graphics, then len(G) gives the number of distinct
        graphics primitives making up that object.

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
        """
        return len(self.__objects)

    def __delitem__(self, i):
        """
        If G is of type Graphics, then del(G[i]) removes the ith distinct
        graphic primitive making up that object.

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
            sage: del(G[2])
            sage: print G
            Graphics object consisting of 2 graphics primitives
            sage: len(G)
            2
        """
        del self.__objects[int(i)]

    def __setitem__(self, i, x):
        """
        You can replace a GraphicPrimitive (point, line, circle, etc...) in
        a Graphics object G with any other GraphicPrimitive

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives

        ::

            sage: p = polygon([[1,3],[2,-2],[1,1],[1,3]]); print p
            Graphics object consisting of 1 graphics primitive

        ::

            sage: G[1] = p[0]
            sage: G    # show the plot
        """
        from sage.plot.primitive import GraphicPrimitive
        if not isinstance(x, GraphicPrimitive):
            raise TypeError, "x must be a GraphicPrimitive"
        self.__objects[int(i)] = x

    def __radd__(self, other):
        """
        Compute and return other + this graphics object.

        This only works when other is a Python int equal to 0. In all other
        cases a TypeError is raised. The main reason for this function is
        to make summing a list of graphics objects easier.

        EXAMPLES::

            sage: S = circle((0,0), 2)
            sage: print int(0) + S
            Graphics object consisting of 1 graphics primitive
            sage: print S + int(0)
            Graphics object consisting of 1 graphics primitive

        The following would fail were it not for this function::

            sage: v = [circle((0,0), 2), circle((2,3), 1)]
            sage: print sum(v)
            Graphics object consisting of 2 graphics primitives
        """
        if isinstance(other, (int, long)) and other == 0:
            return self
        raise TypeError

    def __add__(self, other):
        """
        If you have any Graphics object G1, you can always add any other
        amount of Graphics objects G2,G3,... to form a new Graphics object:
        G4 = G1 + G2 + G3.

        The xmin, xmax, ymin, and ymax properties of the graphics objects
        are expanded to include all objects in both scenes. If the aspect
        ratio property of either or both objects are set, then the larger
        aspect ratio is chosen, with 'automatic' being overridden by a
        numeric aspect ratio.

        If one of the graphics object is set to show a legend, then the
        resulting object will also be set to show a legend.  None of the
        legend options are carried over.

        EXAMPLES::

            sage: g1 = plot(abs(sqrt(x^3-1)), (x,1,5), frame=True)
            sage: g2 = plot(-abs(sqrt(x^3-1)), (x,1,5), color='red')
            sage: g1 + g2  # displays the plot

        TESTS:

        Extra keywords to show are propagated::

            sage: (g1 + g2)._extra_kwds=={'aspect_ratio': 'automatic', 'frame': True}
            True
            sage: g1.set_aspect_ratio(2)
            sage: (g1+g2).aspect_ratio()
            2.0
            sage: g2.set_aspect_ratio(3)
            sage: (g1+g2).aspect_ratio()
            3.0
        """
        if isinstance(other, int) and other == 0:
            return self
        if not isinstance(other, Graphics):
            from sage.plot.plot3d.base import Graphics3d
            if isinstance(other, Graphics3d):
                return self.plot3d() + other
            raise TypeError, "other (=%s) must be a Graphics objects"%other
        g = Graphics()
        g.__objects = self.__objects + other.__objects
        g.__show_legend = self.__show_legend or other.__show_legend
        g._extra_kwds.update(self._extra_kwds)
        g._extra_kwds.update(other._extra_kwds)
        if self.aspect_ratio()=='automatic':
            g.set_aspect_ratio(other.aspect_ratio())
        elif other.aspect_ratio()=='automatic':
            g.set_aspect_ratio(self.aspect_ratio())
        else:
            g.set_aspect_ratio( max(self.aspect_ratio(), other.aspect_ratio()))
        return g

    def add_primitive(self, primitive):
        """
        Adds a primitive to this graphics object.

        EXAMPLES:

        We give a very explicit example::

            sage: G = Graphics()
            sage: from sage.plot.line import Line
            sage: from sage.plot.arrow import Arrow
            sage: L = Line([3,4,2,7,-2],[1,2,e,4,5.],{'alpha':1,'thickness':2,'rgbcolor':(0,1,1),'legend_label':''})
            sage: A = Arrow(2,-5,.1,.2,{'width':3,'head':0,'rgbcolor':(1,0,0),'linestyle':'dashed','zorder':8,'legend_label':''})
            sage: G.add_primitive(L)
            sage: G.add_primitive(A)
            sage: G
        """
        self.__objects.append(primitive)

    def plot(self, *args, **kwds):
        """
        Draw a 2D plot of this graphics object, which just returns this
        object since this is already a 2D graphics object.

        EXAMPLES::

            sage: S = circle((0,0), 2)
            sage: S.plot() is S
            True
        """
        return self

    def plot3d(self, z=0, **kwds):
        """
        Returns an embedding of this 2D plot into the xy-plane of 3D space,
        as a 3D plot object. An optional parameter z can be given to
        specify the z-coordinate.

        EXAMPLES::

            sage: sum([plot(z*sin(x), 0, 10).plot3d(z) for z in range(6)]) # long time
        """
        from sage.plot.plot3d.base import Graphics3dGroup
        g = Graphics3dGroup([g.plot3d(**kwds) for g in self.__objects])
        if z:
            g = g.translate(0,0,z)
        return g

    @classmethod
    def _extract_kwds_for_show(cls, kwds, ignore=[]):
        """
        Extract keywords relevant to show() from the provided dictionary.

        EXAMPLES::

            sage: kwds = {'f': lambda x: x, 'xmin': 0, 'figsize': [1,1], 'plot_points': (40, 40)}
            sage: G_kwds = Graphics._extract_kwds_for_show(kwds, ignore='xmin')
            sage: kwds # Note how this action modifies the passed dictionary
            {'xmin': 0, 'plot_points': (40, 40), 'f': <function <lambda> at ...>}
            sage: G_kwds
            {'figsize': [1, 1]}

        This method is intended to be used with _set_extra_kwds(). Here is an
        idiom to ensure the correct keywords will get passed on to show()::

            sage: options = {} # Usually this will come from an argument
            sage: g = Graphics()
            sage: g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
        """
        result = {}
        for option in cls.SHOW_OPTIONS:
            if option not in ignore:
                try:
                    result[option] = kwds.pop(option)
                except KeyError:
                    pass
        return result

    def _set_extra_kwds(self, kwds):
        """
        Set a dictionary of keywords that will get passed on to show().

        TESTS::

            sage: g = Graphics()
            sage: g._extra_kwds
            {}
            sage: g._set_extra_kwds({'figsize': [10,10]})
            sage: g._extra_kwds
            {'figsize': [10, 10]}
            sage: g.show() # Now the (blank) plot will be extra large
        """
        self._extra_kwds = kwds

    # This dictionary has the default values for the keywords to show(). When
    # show is invoked with keyword arguments, those arguments are merged with
    # this dictionary to create a set of keywords with the defaults filled in.
    # Then, those keywords are passed on to save().

    # NOTE: If you intend to use a new parameter in show(), you should update
    # this dictionary to contain the default value for that parameter.

    SHOW_OPTIONS = dict(xmin=None, xmax=None, ymin=None, ymax=None,
                        figsize=None, fig_tight=True,
                        filename=None,
                        dpi=DEFAULT_DPI, axes=None, axes_labels=None,frame=False,
                        fontsize=None,
                        aspect_ratio=None,
                        gridlines=None, gridlinesstyle=None,
                        vgridlinesstyle=None, hgridlinesstyle=None,transparent=False,
                        show_legend=None, legend_options={},
                        axes_pad=.02, ticks_integer=False,
                        ticks=None, tick_formatter=None)

    @suboptions('legend', numpoints=2, borderpad=0.6, markerscale=0.6, shadow=False,
                labelspacing=0.02, handlelength=0.05, handletextpad=0.5, borderaxespad=None,
                loc='best', font_size='medium', font_family='sans-serif', font_style='normal',
                font_weight='medium', font_variant='normal', back_color=(0.9, 0.9, 0.9),
                title=None, ncol=1, columnspacing=None, fancybox=False)
    def show(self, **kwds):
        """
        Show this graphics image with the default image viewer.

        OPTIONAL INPUT:

        - ``filename`` - (default: None) string

        - ``dpi`` - dots per inch

        - ``figsize`` - [width, height]

        - ``fig_tight`` - (default: True) whether to clip the drawing
          tightly around drawn objects.  If True, then the resulting
          image will usually not have dimensions corresponding to
          ``figsize``.  If False, the resulting image will have
          dimensions corresponding to ``figsize``.

        - ``aspect_ratio`` - the perceived height divided by the
          perceived width. For example, if the aspect ratio is set to ``1``, circles
          will look round and a unit square will appear to have sides
          of equal length, and if the aspect ratio is set ``2``, vertical units will be
          twice as long as horizontal units, so a unit square will be twice as
          high as it is wide.  If set to ``'automatic'``, the aspect ratio
          is determined by ``figsize`` and the picture fills the figure.

        - ``axes`` - (default: True)

        - ``axes_labels`` - (default: None) list (or tuple) of two
          strings; the first is used as the label for the horizontal
          axis, and the second for the vertical axis.

        - ``fontsize`` - (default: current setting -- 10) positive
          integer; used for axes labels; if you make this very large,
          you may have to increase figsize to see all labels.

        - ``frame`` - (default: False) draw a frame around the image

        - ``gridlines`` - (default: None) can be any of the following:

          - None, False: do not add grid lines.

          - True, "automatic", "major": add grid lines at major ticks of the axes.

          - "minor": add grid at major and minor ticks.

          - [xlist,ylist]: a tuple or list containing
            two elements, where xlist (or ylist) can be
            any of the following.


            - None, False: don't add horizontal (or vertical) lines.

            - True, "automatic", "major": add horizontal (or vertical) grid lines at
              the major ticks of the axes.

            - "minor": add horizontal (or vertical) grid lines at major and minor ticks of
              axes.

            - an iterable yielding numbers n or pairs (n,opts), where n
              is the coordinate of the line and opt is a dictionary of
              MATPLOTLIB options for rendering the line.


        - ``gridlinesstyle, hgridlinesstyle, vgridlinesstyle`` -
          (default: None) a dictionary of MATPLOTLIB options for the
          rendering of the grid lines, the horizontal grid lines or the
          vertical grid lines, respectively.

        - ``linkmode`` - (default: False) If True a string containing a link
            to the produced file is returned.

        - ``transparent`` - (default: False) If True, make the background transparent.

        - ``axes_pad`` - (default: 0.02) The percentage of the axis
          range that is added to each end of each axis.  This helps
          avoid problems like clipping lines because of line-width,
          etc.  To get axes that are exactly the specified limits, set
          ``axes_pad`` to zero.

        - ``ticks_integer`` - (default: False) guarantee that the ticks
          are integers (the ``ticks`` option, if specified, will
          override this)

        - ``ticks`` - A matplotlib locator for the major ticks, or
          a number. There are several options.  For more information about
          locators, type ``from matplotlib import ticker`` and then
          ``ticker?``.

          - If this is a locator object, then it is the locator for
            the horizontal axis.  A value of None means use the default
            locator.

          - If it is a list of two locators, then the first is for the
            horizontal axis and one for the vertical axis.  A value of
            None means use the default locator (so a value of
            [None, my_locator] uses my_locator for the vertical axis and
            the default for the horizontal axis).

          - If in either case above one of the entries is a number `m`
            (something which can be coerced to a float), it will be
            replaced by a MultipleLocator which places major ticks at
            integer multiples of `m`.  See examples.

          - If in either case above one of the entries is a list of
            numbers, it will be replaced by a FixedLocator which places
            ticks at the locations specified.  This includes the case of
            of the empty list, which will give no ticks.  See examples.

        - ``tick_formatter`` - A matplotlib formatter for the major
          ticks. There are several options.  For more information about
          formatters, type ``from matplotlib import ticker`` and then
          ``ticker?``.

          If the value of this keyword is a single item, then this will
          give the formatting for the horizontal axis *only* (except for
          the ``"latex"`` option).  If it is a list or tuple, the first
          is for the horizontal axis, the second for the vertical axis.
          The options are below:

          - If one of the entries is a formatter object, then it used.
            A value of None means to use the default locator (so using
            ``tick_formatter=[None, my_formatter]`` uses my_formatter
            for the vertical axis and the default for the horizontal axis).

          - If one of the entries is a symbolic constant such as `\pi`,
            `e`, or `sqrt(2)`, ticks will be formatted nicely at rational
            multiples of this constant.

          .. warning:: This should only be used with the ``ticks`` option
             using nice rational multiples of that constant!

          - If one of the entries is the string ``"latex"``, then the
            formatting will be nice typesetting of the ticks.  This is
            intended to be used when the tick locator for at least one of
            the axes is a list including some symbolic elements.  See examples.

        - ``show_legend`` - (default: None) If True, show the legend

        - ``legend_*`` - all the options valid for :meth:`set_legend_options` prefixed with ``legend_``

        EXAMPLES::

            sage: c = circle((1,1), 1, color='red')
            sage: c.show(xmin=-1, xmax=3, ymin=-1, ymax=3)

        You could also just make the picture larger by changing ``figsize``::

            sage: c.show(figsize=8, xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can turn off the drawing of the axes::

            sage: show(plot(sin,-4,4), axes=False)

        You can also label the axes.  Putting something in dollar
        signs formats it as a mathematical expression::

            sage: show(plot(sin,-4,4), axes_labels=('$x$','$y$'))

        You can turn on the drawing of a frame around the plots::

            sage: show(plot(sin,-4,4), frame=True)

        You can make the background transparent::

            sage: plot(sin(x), (x, -4, 4), transparent=True)

        Add grid lines at the major ticks of the axes.

        ::

            sage: c = circle((0,0), 1)
            sage: c.show(gridlines=True)
            sage: c.show(gridlines="automatic")
            sage: c.show(gridlines="major")

        Add grid lines at the major and minor ticks of the axes.

        ::

            sage: u,v = var('u v')
            sage: f = exp(-(u^2+v^2))
            sage: p = plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2))
            sage: p.show(gridlines="minor")

        Add only horizontal or vertical grid lines.

        ::

            sage: p = plot(sin,-10,20)
            sage: p.show(gridlines=[None, "automatic"])
            sage: p.show(gridlines=["minor", False])

        Add grid lines at specific positions (using lists/tuples).

        ::

            sage: x, y = var('x, y')
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3)-4*(x^2+y^2-2*x)^2, \
            ...             (x,-2,2), (y,-2,2), plot_points=1000)
            sage: p.show(gridlines=[[1,0],[-1,0,1]])

        Add grid lines at specific positions (using iterators).

        ::

            sage: def maple_leaf(t):
            ...     return (100/(100+(t-pi/2)^8))*(2-sin(7*t)-cos(30*t)/2)
            sage: p = polar_plot(maple_leaf, -pi/4, 3*pi/2, color="red",plot_points=1000) # long time
            sage: p.show(gridlines=( [-3,-2.75,..,3], xrange(-1,5,2) )) # long time

        Add grid lines at specific positions (using functions).

        ::

            sage: y = x^5 + 4*x^4 - 10*x^3 - 40*x^2 + 9*x + 36
            sage: p = plot(y, -4.1, 1.1)
            sage: xlines = lambda a,b: [z for z,m in y.roots()]
            sage: p.show(gridlines=[xlines, [0]], frame=True, axes=False)

        Change the style of all the grid lines.

        ::

            sage: b = bar_chart([-3,5,-6,11], color='red')
            sage: b.show(gridlines=([-1,-0.5,..,4],True),
            ...     gridlinesstyle=dict(color="blue", linestyle=":"))

        Change the style of the horizontal or vertical grid lines
        separately.

        ::

            sage: p = polar_plot(2 + 2*cos(x), 0, 2*pi, color=hue(0.3))
            sage: p.show(gridlines=True,
            ...     hgridlinesstyle=dict(color="orange", linewidth=1.0),
            ...     vgridlinesstyle=dict(color="blue", linestyle=":"))

        Change the style of each grid line individually.

        ::

            sage: x, y = var('x, y')
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3)-4*(x^2+y^2-2*x)^2,
            ...             (x,-2,2), (y,-2,2), plot_points=1000)
            sage: p.show(gridlines=(
            ...    [
            ...     (1,{"color":"red","linestyle":":"}),
            ...     (0,{"color":"blue","linestyle":"--"})
            ...    ],
            ...    [
            ...     (-1,{"color":"red","linestyle":":"}),
            ...     (0,{"color":"blue","linestyle":"--"}),
            ...     (1,{"color":"red","linestyle":":"}),
            ...    ]
            ...    ),
            ...    gridlinesstyle=dict(marker='x',color="black"))

        Grid lines can be added to contour plots.

        ::

            sage: f = sin(x^2 + y^2)*cos(x)*sin(y)
            sage: c = contour_plot(f, (x, -4, 4), (y, -4, 4), plot_points=100)
            sage: c.show(gridlines=True, gridlinesstyle={'linestyle':':','linewidth':1, 'color':'red'})

        Grid lines can be added to matrix plots.

        ::

            sage: M = MatrixSpace(QQ,10).random_element()
            sage: matrix_plot(M).show(gridlines=True)

        By default, Sage increases the horizontal and vertical axes
        limits by a certain percentage in all directions.  This is
        controlled by the ``axes_pad`` parameter.  Increasing the range
        of the axes helps avoid problems with lines and dots being
        clipped because the linewidth extends beyond the axes.  To get
        axes limits that are exactly what is specified, set
        ``axes_pad`` to zero.  Compare the following two examples

        ::

            sage: plot(sin(x), (x, -pi, pi),thickness=2)+point((pi, -1), pointsize=15)
            sage: plot(sin(x), (x, -pi, pi),thickness=2,axes_pad=0)+point((pi, -1), pointsize=15)

        Via matplotlib, Sage allows setting of custom ticks.  See above
        for more details.

        ::

            sage: plot(sin(pi*x), (x, -8, 8)) # Labels not so helpful
            sage: plot(sin(pi*x), (x, -8, 8), ticks=2) # Multiples of 2
            sage: plot(sin(pi*x), (x, -8, 8), ticks=[[-7,-3,0,3,7],[-1/2,0,1/2]]) # Your choices
            sage: plot(sin(pi*x), (x, -8, 8), ticks=[[],[]]) # No ticks at all!

        This can be very helpful in showing certain features of plots. ::

            sage: plot(1.5/(1+e^(-x)), (x, -10, 10)) # doesn't quite show value of inflection point

        ::

            sage: plot(1.5/(1+e^(-x)), (x, -10, 10), ticks=[None, 1.5/4]) # It's right at f(x)=0.75!

        But be careful to leave enough room for at least two major ticks, so that
        the user can tell what the scale is.

        ::

            sage: plot(x^2,(x,1,8),ticks=6)
            Traceback (most recent call last):
            ...
            ValueError: Expand the range of the independent variable to allow two multiples of your tick locator (option `ticks`).

        We can also do custom formatting if you need it.  See above for full
        details.

        ::

            sage: plot(2*x+1,(x,0,5),ticks=[[0,1,e,pi,sqrt(20)],2],tick_formatter="latex")

        This is particularly useful when setting custom ticks in multiples
        of `\pi`.

        ::

            sage: plot(sin(x),(x,0,2*pi),ticks=pi/3,tick_formatter=pi)

        But keep in mind that you will get exactly the formatting you asked
        for if you specify both formatters.  The first syntax is recommended
        for best style in that case. ::

            sage: plot(arcsin(x),(x,-1,1),ticks=[None,pi/6],tick_formatter=["latex",pi]) # Nice-looking!

        ::

            sage: plot(arcsin(x),(x,-1,1),ticks=[None,pi/6],tick_formatter=[None,pi]) # Not so nice-looking

        """

        # This option should not be passed on to save().
        linkmode = kwds.pop('linkmode', False)

        if sage.plot.plot.DOCTEST_MODE:
            kwds.pop('filename', None)
            self.save(DOCTEST_MODE_FILE, **kwds)
        elif sage.plot.plot.EMBEDDED_MODE:
            kwds.setdefault('filename', sage.misc.misc.graphics_filename())
            self.save(**kwds)
            if linkmode == True:
                return "<img src='cell://%s'>" % kwds['filename']
            else:
                html("<img src='cell://%s'>" % kwds['filename'])
        else:
            kwds.setdefault('filename', sage.misc.misc.tmp_filename() + '.png')
            self.save(**kwds)
            os.system('%s %s 2>/dev/null 1>/dev/null &'
                      % (sage.misc.viewer.browser(), kwds['filename']))

    def xmin(self, xmin=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.xmin()
            -1.0
            sage: g.xmin(-3)
            sage: g.xmin()
            -3.0
        """
        if xmin is None:
            return self.get_axes_range()['xmin']
        else:
            self.set_axes_range(xmin=xmin)

    def xmax(self, xmax=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.xmax()
            3.0
            sage: g.xmax(10)
            sage: g.xmax()
            10.0
        """
        if xmax is None:
            return self.get_axes_range()['xmax']
        else:
            self.set_axes_range(xmax=xmax)

    def ymin(self, ymin=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.ymin()
            1.0
            sage: g.ymin(-3)
            sage: g.ymin()
            -3.0
        """
        if ymin is None:
            return self.get_axes_range()['ymin']
        else:
            self.set_axes_range(ymin=ymin)

    def ymax(self, ymax=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.ymax()
            2.0
            sage: g.ymax(10)
            sage: g.ymax()
            10.0
        """
        if ymax is None:
            return self.get_axes_range()['ymax']
        else:
            self.set_axes_range(ymax=ymax)


    def get_minmax_data(self):
        """
        Return a dictionary whose keys give the xmin, xmax, ymin, and ymax
        data for this graphic.

        .. warning::

           The returned dictionary is mutable, but changing it does
           not change the xmin/xmax/ymin/ymax data.  The minmax data is a function
           of the primitives which make up this Graphics object.  To change the
           range of the axes, call methods :meth:`xmin`, :meth:`xmax`,
           :meth:`ymin`, :meth:`ymax`, or :meth:`set_axes_range`.

        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: list(sorted(g.get_minmax_data().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 2.0), ('ymin', 1.0)]

        Note that changing ymax doesn't change the output of get_minmax_data::

            sage: g.ymax(10)
            sage: list(sorted(g.get_minmax_data().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 2.0), ('ymin', 1.0)]
        """
        objects = self.__objects
        if objects:
            minmax_data = [o.get_minmax_data() for o in objects]
            xmin = min(d['xmin'] for d in minmax_data)
            xmax = max(d['xmax'] for d in minmax_data)
            ymin = min(d['ymin'] for d in minmax_data)
            ymax = max(d['ymax'] for d in minmax_data)
            # check for NaN's: weird thing -- only way I know to check if a float
            # is a NaN is to check if it is not equal to itself.
            if xmin!=xmin:
                xmin=0; sage.misc.misc.verbose("xmin was NaN (setting to 0)", level=0)
            if xmax!=xmax:
                xmax=0; sage.misc.misc.verbose("xmax was NaN (setting to 0)", level=0)
            if ymin!=ymin:
                ymin=0; sage.misc.misc.verbose("ymin was NaN (setting to 0)", level=0)
            if ymax!=ymax:
                ymax=0; sage.misc.misc.verbose("ymax was NaN (setting to 0)", level=0)
        else:
            xmin = xmax = ymin = ymax = 0

        if xmin == xmax:
            xmin -= 1
            xmax += 1
        if ymin == ymax:
            ymin -= 1
            ymax += 1
        return {'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax}

    def matplotlib(self, filename=None,
                   xmin=None, xmax=None, ymin=None, ymax=None,
                   figsize=None, figure=None, sub=None,
                   axes=None, axes_labels=None, fontsize=None,
                   frame=False, verify=True,
                   aspect_ratio = None,
                   gridlines=None, gridlinesstyle=None,
                   vgridlinesstyle=None, hgridlinesstyle=None,
                   show_legend=None, legend_options={},
                   axes_pad=0.02, ticks_integer=None,
                   tick_formatter=None, ticks=None):
        r"""
        Return a matplotlib figure object representing the graphic

        EXAMPLES::

            sage: c = circle((1,1),1)
            sage: print c.matplotlib()
            Figure(640x480)

        To obtain the first matplotlib axes object inside of the
        figure, you can do something like the following.

        ::

            sage: p=plot(sin(x), (x, -2*pi, 2*pi))
            sage: figure=p.matplotlib()
            sage: axes=figure.axes[0]

        For input parameters, see the documentation for the
        :meth:`show` method (this function accepts all except the
        transparent argument).

        TESTS:

        We verify that :trac:`10291` is fixed::

          sage: p = plot(sin(x), (x, -2*pi, 2*pi))
          sage: figure = p.matplotlib()
          sage: axes_range = p.get_axes_range()
          sage: figure = p.matplotlib()
          sage: axes_range2 = p.get_axes_range()
          sage: axes_range == axes_range2
          True
        """
        if not isinstance(ticks, (list, tuple)):
            ticks = (ticks, None)

        from sage.symbolic.ring import SR
        if not isinstance(tick_formatter, (list, tuple)):  # make sure both formatters typeset or both don't
            if tick_formatter == "latex" or tick_formatter in SR:
                tick_formatter = (tick_formatter, "latex")
            else:
                tick_formatter = (tick_formatter, None)

        self.set_axes_range(xmin, xmax, ymin, ymax)
        d = self.get_axes_range()
        xmin = d['xmin']
        xmax = d['xmax']
        ymin = d['ymin']
        ymax = d['ymax']

        x_pad=(xmax-xmin)*float(axes_pad)
        y_pad=(ymax-ymin)*float(axes_pad)

        xmin-=x_pad
        xmax+=x_pad
        ymin-=y_pad
        ymax+=y_pad

        global do_verify
        do_verify = verify

        if axes is None:
            axes = self.__show_axes

        from matplotlib.figure import Figure
        from matplotlib import rcParams
        self.fontsize(fontsize)
        self.axes_labels(l=axes_labels)

        if figsize is not None and not isinstance(figsize, (list, tuple)):
            default_width, default_height=rcParams['figure.figsize']
            figsize=(figsize, default_height*figsize/default_width)

        if figure is None:
            figure=Figure(figsize=figsize)

        #the incoming subplot instance
        subplot = sub
        if not subplot:
            subplot = figure.add_subplot(111)
        if aspect_ratio is None:
            aspect_ratio=self.aspect_ratio()
        if aspect_ratio == 'automatic':
            subplot.set_aspect('auto', adjustable='box')
        else:
            subplot.set_aspect(aspect_ratio, adjustable='box')
        #add all the primitives to the subplot
        for g in self.__objects:
            g._render_on_subplot(subplot)
            if hasattr(g, '_bbox_extra_artists'):
                self.__bbox_extra_artists.extend(g._bbox_extra_artists)

        #add the legend if requested
        if show_legend is None:
            show_legend = self.__show_legend

        if show_legend:
            from matplotlib.font_manager import FontProperties
            lopts = dict()
            lopts.update(legend_options)
            lopts.update(self.__legend_opts)
            prop = FontProperties(family=lopts.pop('font_family'), weight=lopts.pop('font_weight'), \
                    size=lopts.pop('font_size'), style=lopts.pop('font_style'), variant=lopts.pop('font_variant'))
            color = lopts.pop('back_color')
            leg = subplot.legend(prop=prop, **lopts)
            if leg is None:
                sage.misc.misc.warn("legend requested but no items are labeled")
            else:
                # color
                lframe = leg.get_frame()
                lframe.set_facecolor(color)


        subplot.set_xlim([xmin, xmax])
        subplot.set_ylim([ymin,ymax])

        locator_options=dict(nbins=9,steps=[1,2,5,10],integer=ticks_integer)


        if axes is None:
            axes = self.__show_axes

        for spine in subplot.spines.values():
            spine.set_color(self.__axes_color)
            spine.set_linewidth(self.__axes_width)


        if frame:
            # For now, set the formatter to the old one, since that is
            # sort of what we are used to.  We should eventually look at
            # the default one to see if we like it better.

            from matplotlib.ticker import OldScalarFormatter, MaxNLocator, MultipleLocator, FixedLocator, NullLocator, Locator
            x_locator, y_locator = ticks
            if x_locator is None:
                x_locator = MaxNLocator(**locator_options)
            elif isinstance(x_locator,Locator):
                pass
            elif x_locator == []:
                x_locator = NullLocator()
            elif isinstance(x_locator,list):
                x_locator = FixedLocator(x_locator)
            else: # x_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(xmax/x_locator)-ceil(xmin/x_locator)>1:
                    x_locator=MultipleLocator(float(x_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the independent variable to allow two multiples of your tick locator (option `ticks`).')
            if y_locator is None:
                y_locator = MaxNLocator(**locator_options)
            elif isinstance(y_locator,Locator):
                pass
            elif y_locator == []:
                y_locator = NullLocator()
            elif isinstance(y_locator,list):
                y_locator = FixedLocator(y_locator)
            else: # y_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(ymax/y_locator)-ceil(ymin/y_locator)>1:
                    y_locator=MultipleLocator(float(y_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the dependent variable to allow two multiples of your tick locator (option `ticks`).')

            x_formatter, y_formatter = tick_formatter
            from matplotlib.ticker import FuncFormatter
            from sage.misc.latex import latex
            if x_formatter is None:
                x_formatter = OldScalarFormatter()
            elif x_formatter in SR:
                from misc import _multiple_of_constant
                x_const = x_formatter
                x_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,x_const))
            elif x_formatter == "latex":
                x_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))
            if y_formatter is None:
                y_formatter = OldScalarFormatter()
            elif y_formatter in SR:
                from misc import _multiple_of_constant
                y_const = y_formatter
                y_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,y_const))
            elif y_formatter == "latex":
                y_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))

            subplot.xaxis.set_major_locator(x_locator)
            subplot.yaxis.set_major_locator(y_locator)
            subplot.xaxis.set_major_formatter(x_formatter)
            subplot.yaxis.set_major_formatter(y_formatter)

            subplot.set_frame_on(True)
            if axes:
                if ymin<=0 and ymax>=0:
                    subplot.axhline(color=self.__axes_color,
                                    linewidth=self.__axes_width)
                if xmin<=0 and xmax>=0:
                    subplot.axvline(color=self.__axes_color,
                                    linewidth=self.__axes_width)

        elif axes:
            ymiddle=False
            xmiddle=False
            if xmin>0:
                subplot.spines['right'].set_visible(False)
                subplot.spines['left'].set_position(('outward',10))
                subplot.yaxis.set_ticks_position('left')
                subplot.yaxis.set_label_position('left')
                yaxis='left'
            elif xmax<0:
                subplot.spines['left'].set_visible(False)
                subplot.spines['right'].set_position(('outward',10))
                subplot.yaxis.set_ticks_position('right')
                subplot.yaxis.set_label_position('right')
                yaxis='right'
            else:
                subplot.spines['left'].set_position('zero')
                subplot.yaxis.set_ticks_position('left')
                subplot.yaxis.set_label_position('left')
                subplot.spines['right'].set_visible(False)
                ymiddle=True
                yaxis='left'

            if ymin>0:
                subplot.spines['top'].set_visible(False)
                subplot.spines['bottom'].set_position(('outward',10))
                subplot.xaxis.set_ticks_position('bottom')
                subplot.xaxis.set_label_position('bottom')
                xaxis='bottom'
            elif ymax<0:
                subplot.spines['bottom'].set_visible(False)
                subplot.spines['top'].set_position(('outward',10))
                subplot.xaxis.set_ticks_position('top')
                subplot.xaxis.set_label_position('top')
                xaxis='top'
            else:
                subplot.spines['bottom'].set_position('zero')
                subplot.xaxis.set_ticks_position('bottom')
                subplot.xaxis.set_label_position('bottom')
                subplot.spines['top'].set_visible(False)
                xmiddle=True
                xaxis='bottom'

            # For now, set the formatter to the old one, since that is
            # sort of what we are used to.  We should eventually look at
            # the default one to see if we like it better.

            from matplotlib.ticker import OldScalarFormatter, MaxNLocator, MultipleLocator, FixedLocator, NullLocator, Locator
            x_locator, y_locator = ticks
            if x_locator is None:
                x_locator = MaxNLocator(**locator_options)
            elif isinstance(x_locator,Locator):
                pass
            elif x_locator == []:
                x_locator = NullLocator()
            elif isinstance(x_locator,list):
                x_locator = FixedLocator(x_locator)
            else: # x_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(xmax/x_locator)-ceil(xmin/x_locator)>1:
                    x_locator=MultipleLocator(float(x_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the independent variable to allow two multiples of your tick locator (option `ticks`).')
            if y_locator is None:
                y_locator = MaxNLocator(**locator_options)
            elif isinstance(y_locator,Locator):
                pass
            elif y_locator == []:
                y_locator = NullLocator()
            elif isinstance(y_locator,list):
                y_locator = FixedLocator(y_locator)
            else: # y_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(ymax/y_locator)-ceil(ymin/y_locator)>1:
                    y_locator=MultipleLocator(float(y_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the dependent variable to allow two multiples of your tick locator (option `ticks`).')

            x_formatter, y_formatter = tick_formatter
            from matplotlib.ticker import FuncFormatter
            from sage.misc.latex import latex
            from sage.symbolic.ring import SR
            if x_formatter is None:
                x_formatter = OldScalarFormatter()
            elif x_formatter in SR:
                from misc import _multiple_of_constant
                x_const = x_formatter
                x_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,x_const))
            elif x_formatter == "latex":
                x_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))
            if y_formatter is None:
                y_formatter = OldScalarFormatter()
            elif y_formatter in SR:
                from misc import _multiple_of_constant
                y_const = y_formatter
                y_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,y_const))
            elif y_formatter == "latex":
                y_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))

            subplot.xaxis.set_major_locator(x_locator)
            subplot.yaxis.set_major_locator(y_locator)
            subplot.xaxis.set_major_formatter(x_formatter)
            subplot.yaxis.set_major_formatter(y_formatter)

            # Make ticklines go on both sides of the axes
            #             if xmiddle:
            #                 for t in subplot.xaxis.get_majorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(8)
            #                 for t in subplot.xaxis.get_minorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(4)

            #             if ymiddle:
            #                 for t in subplot.yaxis.get_majorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(8)
            #                 for t in subplot.yaxis.get_minorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(4)

            # Make the zero tick labels disappear if the axes cross
            # inside the picture
            if xmiddle and ymiddle:
                from sage.plot.plot import SelectiveFormatter
                subplot.yaxis.set_major_formatter(SelectiveFormatter(subplot.yaxis.get_major_formatter(),skip_values=[0]))
                subplot.xaxis.set_major_formatter(SelectiveFormatter(subplot.xaxis.get_major_formatter(),skip_values=[0]))

        else:
            for spine in subplot.spines.values():
                spine.set_visible(False)
            from matplotlib.ticker import NullFormatter, NullLocator
            subplot.xaxis.set_major_formatter(NullFormatter())
            subplot.yaxis.set_major_formatter(NullFormatter())
            subplot.xaxis.set_major_locator(NullLocator())
            subplot.yaxis.set_major_locator(NullLocator())

        if frame or axes:
            # Make minor tickmarks, unless we specify fixed ticks or no ticks
            from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullLocator
            if isinstance(x_locator, (NullLocator, FixedLocator)):
                subplot.xaxis.set_minor_locator(NullLocator())
            else:
                subplot.xaxis.set_minor_locator(AutoMinorLocator())
            if isinstance(y_locator, (NullLocator, FixedLocator)):
                subplot.yaxis.set_minor_locator(NullLocator())
            else:
                subplot.yaxis.set_minor_locator(AutoMinorLocator())

            ticklabels=subplot.xaxis.get_majorticklabels() + \
                subplot.xaxis.get_minorticklabels() + \
                subplot.yaxis.get_majorticklabels() + \
                subplot.yaxis.get_minorticklabels()
            for ticklabel in ticklabels:
                ticklabel.set_fontsize(self.__fontsize)
                ticklabel.set_color(self.__tick_label_color)

            ticklines=subplot.xaxis.get_majorticklines() + \
                subplot.xaxis.get_minorticklines() + \
                subplot.yaxis.get_majorticklines() + \
                subplot.yaxis.get_minorticklines()
            for tickline in ticklines:
                tickline.set_color(self.__axes_color)


        if gridlines is not None:
            if isinstance(gridlines, (list, tuple)):
                vgridlines,hgridlines=gridlines
            else:
                hgridlines=gridlines
                vgridlines=gridlines

            if gridlinesstyle is None:
                # Set up the default grid style
                gridlinesstyle=dict(color='black',linestyle=':',linewidth=0.5)

            vgridstyle=gridlinesstyle.copy()
            if vgridlinesstyle is not None:
                vgridstyle.update(vgridlinesstyle)

            hgridstyle=gridlinesstyle.copy()
            if hgridlinesstyle is not None:
                hgridstyle.update(hgridlinesstyle)

            if hgridlines=='minor':
                hgridstyle['which']='both'
            if vgridlines=='minor':
                vgridstyle['which']='both'

            if hasattr(hgridlines, '__iter__'):
                hlines=iter(hgridlines)
                hgridstyle.pop("minor",None)
                for hline in hlines:
                    if isinstance(hline, (list, tuple)):
                        hl, style=hline
                        st=hgridstyle.copy()
                        st.update(style)
                    else:
                        hl=hline
                        st=hgridstyle
                    subplot.axhline(hl,**st)
            else:
                if hgridlines not in (None, False):
                    subplot.yaxis.grid(True, **hgridstyle)

            if hasattr(vgridlines, '__iter__'):
                vlines=iter(vgridlines)
                vgridstyle.pop("minor",None)
                for vline in vlines:
                    if isinstance(vline, (list, tuple)):
                        vl, style=vline
                        st=vgridstyle.copy()
                        st.update(style)
                    else:
                        vl=vline
                        st=vgridstyle
                    subplot.axvline(vl,**st)
            else:
                if vgridlines not in (None, False):
                    subplot.xaxis.grid(True, **vgridstyle)



        if self.__axes_labels is not None:
            label_options={}
            label_options['color']=self.__axes_label_color
            label_options['size']=self.__fontsize
            subplot.set_xlabel(self.__axes_labels[0], **label_options)
            subplot.set_ylabel(self.__axes_labels[1], **label_options)


            if axes is True and frame is False:
                # We set the label positions according to where we are
                # drawing the axes.
                if xaxis=='bottom':
                    yaxis_labely=subplot.get_ylim()[1]
                    yaxis_labeloffset=8
                    yaxis_vert='bottom'
                    xaxis_labely=0
                    xaxis_vert='baseline'
                else:
                    yaxis_labely=subplot.get_ylim()[0]
                    yaxis_labeloffset=-8
                    yaxis_vert='top'
                    xaxis_labely=1
                    xaxis_vert='top'

                if yaxis=='left':
                    xaxis_labelx=subplot.get_xlim()[1]
                    xaxis_labeloffset=8
                    xaxis_horiz='left'
                    yaxis_labelx=0
                else:
                    xaxis_labelx=subplot.get_xlim()[0]
                    xaxis_labeloffset=-8
                    xaxis_horiz='right'
                    yaxis_labelx=1

                from matplotlib.transforms import offset_copy
                xlabel=subplot.xaxis.get_label()
                xlabel.set_horizontalalignment(xaxis_horiz)
                xlabel.set_verticalalignment(xaxis_vert)
                trans=subplot.spines[xaxis].get_transform()
                labeltrans=offset_copy(trans, figure, x=xaxis_labeloffset, y=0, units='points')
                subplot.xaxis.set_label_coords(x=xaxis_labelx,y=xaxis_labely,transform=labeltrans)

                ylabel=subplot.yaxis.get_label()
                ylabel.set_horizontalalignment('center')
                ylabel.set_verticalalignment(yaxis_vert)
                ylabel.set_rotation('horizontal')
                trans=subplot.spines[yaxis].get_transform()
                labeltrans=offset_copy(trans, figure, x=0, y=yaxis_labeloffset, units='points')
                subplot.yaxis.set_label_coords(x=yaxis_labelx,y=yaxis_labely,transform=labeltrans)

        # This option makes the xlim and ylim limits not take effect
        # todo: figure out which limits were specified, and let the
        # free limits autoscale
        #subplot.autoscale_view(tight=True)
        return figure

    # ALLOWED_EXTENSIONS is the list of recognized formats.
    # filename argument is written explicitly so that it can be used as a
    # positional one, which is a very likely usage for this function.
    @suboptions('legend', numpoints=2, borderpad=0.6, markerscale=0.6, shadow=False,
                labelspacing=0.02, handlelength=0.05, handletextpad=0.5, borderaxespad=None,
                loc='best', font_size='medium', font_family='sans-serif', font_style='normal',
                font_weight='medium', font_variant='normal', back_color=(0.9, 0.9, 0.9),
                title=None, ncol=1, columnspacing=None, fancybox=False)
    def save(self, filename=None, **kwds):
        r"""
        Save the graphics to an image file.

        INPUT:

        - ``filename`` -- a string (default: autogenerated), the filename and
          the image format given by the extension, which can be one of the
          following:

            * ``.eps``,

            * ``.pdf``,

            * ``.png``,

            * ``.ps``,

            * ``.sobj`` (for a Sage object you can load later),

            * ``.svg``,

            * empty extension will be treated as ``.sobj``.

        All other keyword arguments will be passed to the plotter.

        OUTPUT:

        - none.

        EXAMPLES::

            sage: c = circle((1,1), 1, color='red')
            sage: filename = os.path.join(SAGE_TMP, 'test.png')
            sage: c.save(filename, xmin=-1, xmax=3, ymin=-1, ymax=3)

        To make a figure bigger or smaller, use ``figsize``::

            sage: c.save(filename, figsize=5, xmin=-1, xmax=3, ymin=-1, ymax=3)

        By default, the figure grows to include all of the graphics and text,
        so the final image may not be exactly the figure size you specified.
        If you want a figure to be exactly a certain size, specify the keyword
        ``fig_tight=False``::

            sage: c.save(filename, figsize=[8,4], fig_tight=False,
            ...       xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can also pass extra options to the plot command instead of this
        method, e.g. ::

            sage: plot(x^2 - 5, (x, 0, 5), ymin=0).save(
            ...       sage.misc.misc.tmp_filename() + '.png')

        will save the same plot as the one shown by this command::

            sage: plot(x^2 - 5, (x, 0, 5), ymin=0)

        (This test verifies that :trac:`8632` is fixed.)

        TESTS:

        Legend labels should save correctly::

            sage: P = plot(x,(x,0,1),legend_label='$xyz$')
            sage: P.set_legend_options(back_color=(1,0,0))
            sage: P.set_legend_options(loc=7)
            sage: filename=os.path.join(SAGE_TMP, 'test.png')
            sage: P.save(filename)

        This plot should save with the frame shown, showing :trac:`7524`
        is fixed (same issue as :trac:`7981` and :trac:`8632`)::

            sage: var('x,y')
            (x, y)
            sage: a = plot_vector_field((x,-y),(x,-1,1),(y,-1,1))
            sage: filename=os.path.join(SAGE_TMP, 'test2.png')
            sage: a.save(filename)
        """
        options = dict()
        options.update(self.SHOW_OPTIONS)
        options.update(self._extra_kwds)
        options.update(kwds)
        dpi = options.pop('dpi')
        transparent = options.pop('transparent')
        fig_tight = options.pop('fig_tight')

        if filename is None:
            filename = options.pop('filename')
        if filename is None:
            filename = sage.misc.misc.graphics_filename()
        ext = os.path.splitext(filename)[1].lower()

        if ext not in ALLOWED_EXTENSIONS:
            raise ValueError("allowed file extensions for images are '"
                             + "', '".join(ALLOWED_EXTENSIONS) + "'!")
        elif ext in ['', '.sobj']:
            SageObject.save(self, filename)
        else:
            figure = self.matplotlib(**options)
            # You can output in PNG, PS, EPS, PDF, or SVG format, depending on the file extension.
            # matplotlib looks at the file extension to see what the renderer should be.
            # The default is FigureCanvasAgg for PNG's because this is by far the most
            # common type of files rendered, like in the notebook, for example.
            # if the file extension is not '.png', then matplotlib will handle it.
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            figure.set_canvas(FigureCanvasAgg(figure))
            # this messes up the aspect ratio!
            #figure.canvas.mpl_connect('draw_event', pad_for_tick_labels)

            # tight_layout adjusts the *subplot* parameters so ticks aren't cut off, etc.
            figure.tight_layout()

            if fig_tight is True:
                figure.savefig(filename, dpi=dpi, bbox_inches='tight',
                    bbox_extra_artists=self.__bbox_extra_artists,
                    transparent=transparent)
            else:
                figure.savefig(filename, dpi=dpi,
                           transparent=transparent)


class GraphicsArray(SageObject):
    """
    GraphicsArray takes a (`m` x `n`) list of lists of
    graphics objects and plots them all on one canvas.
    """
    def __init__(self, array):
        """
        Constructor for ``GraphicsArray`` class.  Normally used only
        via :func:`graphics_array` function.

        INPUT: a list or list of lists/tuples, all of which are graphics objects

        EXAMPLES::

            sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in range(10)]
            sage: G = graphics_array(L)
            sage: G.ncols()
            10
            sage: M = [[plot(x^2)],[plot(x^3)]]
            sage: H = graphics_array(M)
            sage: str(H[1])
            'Graphics object consisting of 1 graphics primitive'

        TESTS::

            sage: L = [[plot(sin),plot(cos)],[plot(tan)]]
            sage: graphics_array(L)
            Traceback (most recent call last):
            ...
            TypeError: array (=[[, ], []]) must be a list of lists of Graphics objects
            sage: G = plot(x,(x,0,1))
            sage: graphics_array(G)
            Traceback (most recent call last):
            ...
            TypeError: array (=Graphics object consisting of 1 graphics primitive) must be a list of lists of Graphics objects
            sage: G = [[plot(x,(x,0,1)),x]]
            sage: graphics_array(G)
            Traceback (most recent call last):
            ...
            TypeError: every element of array must be a Graphics object
        """
        if not isinstance(array, (list, tuple)):
            raise TypeError,"array (=%s) must be a list of lists of Graphics objects"%(array)
        array = list(array)
        self._glist = []
        self._rows = len(array)
        if self._rows > 0:
            if not isinstance(array[0], (list, tuple)):
                array = [array]
                self._rows = 1
            self._cols = len(array[0])
        else:
            self._cols = 0
        self._dims = self._rows*self._cols
        for row in array: #basically flatten the list
            if not isinstance(row, (list, tuple)) or len(row) != self._cols:
                raise TypeError,"array (=%s) must be a list of lists of Graphics objects"%(array)
            for g in row:
                if not isinstance(g, Graphics):
                    raise TypeError, "every element of array must be a Graphics object"
                self._glist.append(g)
        self._figsize = None

    def _repr_(self):
        """
        Representation of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: G # plot shown is default (indirect doctest)

        We can make commands not display their plots by default. ::

            sage: show_default(False)
            sage: graphics_array(L) # indirect doctest
            Graphics Array of size 1 x 6
            sage: show_default(True)
        """
        if SHOW_DEFAULT:
            self.show()
            return ''
        else:
            return self.__str__()

    def __str__(self):
        """
        String representation of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: str(G)
            'Graphics Array of size 2 x 3'

        We can make commands not display their plots by default. ::

            sage: show_default(False)
            sage: graphics_array(L)
            Graphics Array of size 1 x 6
            sage: show_default(True)
        """
        return "Graphics Array of size %s x %s"%(self._rows, self._cols)

    def nrows(self):
        """
        Number of rows of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: G.nrows()
            2
            sage: graphics_array(L).nrows()
            1
        """
        return self._rows

    def ncols(self):
        """
        Number of columns of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: G.ncols()
            3
            sage: graphics_array(L).ncols()
            6
        """
        return self._cols

    def __getitem__(self, i):
        """
        Return the ``i``th element of the list of graphics
        in the (flattened) array.

        EXAMPLES:

        We can access and view individual plots::

            sage: M = [[plot(x^2)],[plot(x^3)]]
            sage: H = graphics_array(M)
            sage: H[1]

        They can also be represented::

            sage: str(H[1])
            'Graphics object consisting of 1 graphics primitive'

        Another example::

            sage: L = [plot(sin(k*x),(x,-pi,pi))+circle((k,k),1,color='red') for k in range(10)]
            sage: G = graphics_array(L,5,2)
            sage: str(G[3])
            'Graphics object consisting of 2 graphics primitives'
            sage: G[3]
        """
        i = int(i)
        return self._glist[i]

    def __setitem__(self, i, g):
        """
        Set the ``i``th element of the list of graphics
        in the (flattened) array.

        EXAMPLES::

            sage: M = [[plot(x^2)],[plot(x^3)]]
            sage: H = graphics_array(M)
            sage: str(H[1])
            'Graphics object consisting of 1 graphics primitive'

        We can check this is one primitive::

            sage: H[1] # the plot of x^3

        Now we change it::

            sage: H[1] = circle((1,1),2)+points([(1,2),(3,2),(5,5)],color='purple')
            sage: str(H[1])
            'Graphics object consisting of 2 graphics primitives'

        And we visually check that it's different::

            sage: H[1] # a circle and some purple points
        """
        i = int(i)
        self._glist[i] = g

    def __set_figsize__(self, ls):
        """
        Set the figsize of all plots in the array.

        This is normally only used via the ``figsize`` keyword in
        :meth:`save` or :meth:`show`.

        EXAMPLES::

            sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in [1..3]]
            sage: G = graphics_array(L)
            sage: G.show(figsize=[5,3])  # smallish and compact

        ::

            sage: G.show(figsize=[10,20])  # bigger and tall and thin; long time (2s on sage.math, 2012)

        ::

            sage: G.show(figsize=8)  # figure as a whole is a square
        """
        # if just one number is passed in for figsize, as documented
        if not isinstance(ls,list):
            ls = [ls,ls]
        # now the list is a list
        m = int(ls[0])
        n = int(ls[1])
        self._figsize = [m,n]

    def __len__(self):
        """
        Total number of elements of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: G.ncols()
            3
            sage: graphics_array(L).ncols()
            6
        """
        return len(self._glist)

# This does not work, and can never have worked!
# To make this work, one would also change the
# dimensions of the array, but it's not clear there
# is a canonical way to do this.
#
#    def append(self, g):
#        """
#        Appends a graphic to the array.
#        """
#        self._glist.append(g)

    def append(self, g):
        """
        Appends a graphic to the array.  Currently
        not implemented.

        TESTS::

            sage: from sage.plot.graphics import GraphicsArray
            sage: G = GraphicsArray([plot(sin),plot(cos)])
            sage: G.append(plot(tan))
            Traceback (most recent call last):
            ...
            NotImplementedError: Appending to a graphics array is not yet implemented
        """
        raise NotImplementedError('Appending to a graphics array is not yet implemented')


    def _render(self, filename, dpi=None, figsize=None, axes=None, **args):
        r"""
        ``_render`` loops over all graphics objects in the array
        and adds them to the subplot.  This is only used internally
        when the plot is actually saved or shown.

        EXAMPLES::

            sage: graphics_array([[plot(sin), plot(cos)], [plot(tan), plot(sec)]])
        """
        #glist is a list of Graphics objects:
        glist = self._glist
        rows = self._rows
        cols = self._cols
        dims = self._dims
        #make a blank matplotlib Figure:
        from matplotlib.figure import Figure
        figure = Figure(figsize)
        global do_verify
        do_verify = True
        for i,g in zip(range(1, dims+1), glist):
            subplot = figure.add_subplot(rows, cols, i)
            g.matplotlib(filename, figure=figure, sub=subplot,
                         verify=do_verify, axes = axes, **args)
        g.save(filename, dpi=dpi, figure=figure, sub=subplot,
               verify=do_verify, axes = axes, **args)

    def save(self, filename=None, dpi=DEFAULT_DPI, figsize=None,
             axes = None, **args):
        """
        Save the ``graphics_array`` to (for now) a png called
        'filename'.

        OPTIONAL INPUT:

        -  ``filename`` - (default: None) string

        -  ``dpi`` - dots per inch

        -  ``figsize`` - width or [width, height]

        -  ``axes`` - (default: True)

        EXAMPLES::

            sage: F = sage.misc.misc.tmp_filename()+'.png'
            sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in [1..3]]
            sage: G = graphics_array(L)
            sage: G.save(F,500,axes=False)  # long time (6s on sage.math, 2012)
        """
        if (figsize is not None): self.__set_figsize__(figsize)
        self._render(filename, dpi=dpi, figsize=self._figsize, axes = axes, **args)

    def show(self, filename=None, dpi=DEFAULT_DPI, figsize=None,
             axes = None, **args):
        r"""
        Show this graphics array using the default viewer.

        OPTIONAL INPUT:

        -  ``filename`` - (default: None) string

        -  ``dpi`` - dots per inch

        -  ``figsize`` - width or [width, height]

        -  ``axes`` - (default: True)

        -  ``fontsize`` - positive integer

        -  ``frame`` - (default: False) draw a frame around the
           image

        EXAMPLES: This draws a graphics array with four trig plots and no
        axes in any of the plots.

        ::

            sage: G = graphics_array([[plot(sin), plot(cos)], [plot(tan), plot(sec)]])
            sage: G.show(axes=False)
        """
        if (figsize is not None): self.__set_figsize__(figsize)
        if sage.plot.plot.DOCTEST_MODE:
            self.save(DOCTEST_MODE_FILE,
                      dpi=dpi, figsize=self._figsize, axes = axes, **args)
            return
        if sage.plot.plot.EMBEDDED_MODE:
            self.save(filename, dpi=dpi, figsize=self._figsize, axes = axes, **args)
            return
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self._render(filename, dpi=dpi, figsize=self._figsize, axes = axes, **args)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(
                         sage.misc.viewer.browser(), filename))


