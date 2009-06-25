"""
Text in plots
"""
#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#
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
from sage.plot.primitive import GraphicPrimitive
from sage.plot.misc import options, rename_keyword
from sage.plot.colors import to_mpl_color

class Text(GraphicPrimitive):
    """
    Base class for Text graphics primitive.

    TESTS:

    We test creating some text::

        sage: text("I like Fibonacci",(3,5))
    """
    def __init__(self, string, point, options):
        """
        Initializes base class Text.

        EXAMPLES::

            sage: T = text("I like Fibonacci", (3,5))
            sage: t = T[0]
            sage: t.string
            'I like Fibonacci'
            sage: t.x
            3.0
            sage: t.options()['fontsize']
            10
        """
        self.string = string
        self.x = float(point[0])
        self.y = float(point[1])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data. Notice
        that, for text, the box is just the location itself.

        EXAMPLES::

            sage: T = text("Where am I?",(1,1))
            sage: t=T[0]
            sage: t.get_minmax_data()['ymin']
            1.0
            sage: t.get_minmax_data()['ymax']
            1.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data([self.x], [self.y], dict=True)

    def _repr_(self):
        """
        String representation of Text primitive.

        EXAMPLES::

            sage: T = text("I like cool constants", (pi,e))
            sage: t=T[0];t
            Text 'I like cool constants' at the point (3.14159265359,2.71828182846)
        """
        return "Text '%s' at the point (%s,%s)"%(self.string, self.x, self.y)

    def _allowed_options(self):
        """
        Return the allowed options for the Text class.

        EXAMPLES::

            sage: T = text("ABC",(1,1),zorder=3)
            sage: T[0]._allowed_options()['fontsize']
            'How big the text is.'
            sage: T[0]._allowed_options()['zorder']
            'The layer level in which to draw'
        """
        return {'fontsize': 'How big the text is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'axis_coords':'Uses axis coordinates -- (0,0) lower left and (1,1) upper right',
                'vertical_alignment': 'how to align vertically: top, center, bottom',
                'horizontal_alignment':'how to align horizontally: left, center, right',
                'zorder':'The layer level in which to draw'}

    def _plot3d_options(self, options=None):
        """
        Translate 2D plot options into 3D plot options.

        EXAMPLES::

            sage: T = text("ABC",(1,1))
            sage: t = T[0]
            sage: t.options()['rgbcolor']
            (0.0, 0.0, 1.0)
            sage: s=t.plot3d()
            sage: s.jmol_repr(s.testing_render_params())[0][1]
            'color atom  [0,0,255]'
        """
        if options == None:
            options = dict(self.options())
        options_3d = {}
        # TODO: figure out how to implement rather than ignore
        if 'fontsize' in options:
            del options['fontsize']
        if 'vertical_alignment' in options:
            del options['vertical_alignment']
        if 'horizontal_alignment' in options:
            del options['horizontal_alignment']
        if 'axis_coords' in options:
            del options['axis_coords']
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        Plots 2D text in 3D.

        EXAMPLES::

            sage: T = text("ABC",(1,1))
            sage: t = T[0]
            sage: s=t.plot3d()
            sage: s.jmol_repr(s.testing_render_params())[0][2]
            'label "ABC"'
            sage: s._trans
            (1.0, 1.0, 0)
        """
        from sage.plot.plot3d.shapes2 import text3d
        options = self._plot3d_options()
        options.update(kwds)
        return text3d(self.string, (self.x, self.y, 0), **options)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: t1 = text("Hello",(1,1), vertical_alignment="top", fontsize=30, rgbcolor='black')
            sage: t2 = text("World", (1,1), horizontal_alignment="left",fontsize=20, zorder=-1)
            sage: t1 + t2   # render the sum
        """
        options = self.options()
        opts = {}
        opts['color'] = options['rgbcolor']
        opts['fontsize'] = int(options['fontsize'])
        opts['verticalalignment'] = options['vertical_alignment']
        opts['horizontalalignment'] = options['horizontal_alignment']
        if 'zorder' in options:
            opts['zorder'] = options['zorder']
        if options['axis_coords']:
            opts['transform'] = subplot.transAxes
        subplot.text(self.x, self.y, self.string, **opts)


@options(fontsize=10, rgbcolor=(0,0,1), horizontal_alignment='center',
         vertical_alignment='center', axis_coords=False)
def text(string, xy, **options):
    r"""
    Returns a 2D text graphics object at the point `(x,y)`.

    Type ``text.options`` for a dictionary of options for 2D text.

    2D OPTIONS:

    - ``fontsize`` - How big the text is

    - ``rgbcolor`` - The color as an RGB tuple

    - ``hue`` - The color given as a hue

    - ``vertical_alignment`` - How to align vertically: top, center, bottom

    - ``horizontal_alignment`` - How to align horizontally: left, center, right

    - ``axis_coords`` - (default: False) if True, use axis coordinates, so that
      (0,0) is the lower left and (1,1) upper right, regardless of the x and y
      range of plotted values.

    EXAMPLES::

        sage: text("Sage is really neat!!",(2,12))

    The same text in larger font and colored red::

        sage: text("Sage is really neat!!",(2,12),fontsize=20,rgbcolor=(1,0,0))

    Same text but guaranteed to be in the lower left no matter what::

        sage: text("Sage is really neat!!",(0,0), axis_coords=True, horizontal_alignment='left')

    You can also align text differently::

        sage: t1 = text("Hello",(1,1), vertical_alignment="top")
        sage: t2 = text("World", (1,0.5), horizontal_alignment="left")
        sage: t1 + t2   # render the sum

    You can save text as part of PDF output::

        sage: text("sage", (0,0), rgbcolor=(0,0,0)).save(SAGE_TMP + 'a.pdf')

    Text must be 2D (use the text3d command for 3D text)::

        sage: t = text("hi",(1,2,3))
        Traceback (most recent call last):
        ...
        ValueError: use text3d instead for text in 3d
        sage: t = text3d("hi",(1,2,3))
    """
    try:
        x, y = xy
    except ValueError:
        if isinstance(xy, (list, tuple)) and len(xy) == 3:
            raise ValueError, "use text3d instead for text in 3d"
        raise
    from sage.plot.plot import Graphics
    options['rgbcolor'] = to_mpl_color(options['rgbcolor'])
    point = (float(x), float(y))
    g = Graphics()
    g.add_primitive(Text(string, point, options))
    return g
