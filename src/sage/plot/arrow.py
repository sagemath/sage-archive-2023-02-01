#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Emily Kirkman
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
from sage.plot.misc import options, rename_keyword, to_mpl_color

class CurveArrow(GraphicPrimitive):
    def __init__(self, path, options):
        """
        Returns an arrow graphics primitive along the provided path (bezier curve).

        EXAMPLES:
            sage: from sage.plot.arrow import CurveArrow
            sage: b = CurveArrow(path=[[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]],options={})
            sage: b
            CurveArrow from (0, 0) to (0, 0)
        """
        import numpy as np
        self.path = path
        codes = [1] + (len(self.path[0])-1)*[len(self.path[0])]
        vertices = self.path[0]
        for curve in self.path[1:]:
            vertices += curve
            codes += (len(curve))*[len(curve)+1]
        self.codes = codes
        self.vertices = np.array(vertices, np.float)
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES:
            sage: from sage.plot.arrow import CurveArrow
            sage: b = CurveArrow(path=[[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]],options={})
            sage: d = b.get_minmax_data()
            sage: d['xmin']
            0.0
            sage: d['xmax']
            1.0
        """
        return {'xmin': self.vertices[:,0].min(),
                'xmax': self.vertices[:,0].max(),
                'ymin': self.vertices[:,1].min(),
                'ymax': self.vertices[:,1].max()}

    def _allowed_options(self):
        """
        Return the dictionary of allowed options for the curve arrow graphics primitive.

        EXAMPLES:
             sage: from sage.plot.arrow import CurveArrow
             sage: list(sorted(CurveArrow(path=[[(0,0),(2,3)]],options={})._allowed_options().iteritems()))
             [('arrowsize', 'The size of the arrowhead'),
             ('arrowstyle', 'todo'),
             ('head', '2-d only: Which end of the path to draw the head (one of 0 (start), 1 (end) or 2 (both)'),
             ('hue', 'The color given as a hue.'),
             ('linestyle', "2d only: The style of the line, which is one of 'dashed', 'dotted', 'solid', 'dashdot'."),
             ('rgbcolor', 'The color as an rgb tuple.'),
             ('width', 'The width of the shaft of the arrow, in points.'),
             ('zorder', '2-d only: The layer level in which to draw')]
        """
        return {'width':'The width of the shaft of the arrow, in points.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'arrowstyle': 'todo',
                'linestyle': 'todo',
                'arrowsize':'The size of the arrowhead',
                'zorder':'2-d only: The layer level in which to draw',
                'head':'2-d only: Which end of the path to draw the head (one of 0 (start), 1 (end) or 2 (both)',
                'linestyle':"2d only: The style of the line, which is one of 'dashed', 'dotted', 'solid', 'dashdot'."}

    def _repr_(self):
        """
        Text representation of an arrow graphics primitive.

        EXAMPLES:
            sage: from sage.plot.arrow import CurveArrow
            sage: CurveArrow(path=[[(0,0),(1,4),(2,3)]],options={})._repr_()
            'CurveArrow from (0, 0) to (2, 3)'
        """
        return "CurveArrow from %s to %s"%(self.path[0][0],self.path[-1][-1])

    def _render_on_subplot(self, subplot):
        """
        Render this arrow in a subplot.  This is the key function that
        defines how this arrow graphics primitive is rendered in
        matplotlib's library.

        EXAMPLES:
        This function implicitly ends up rendering this arrow on a matplotlib subplot:
            sage: arrow(path=[[(0,1), (2,-1), (4,5)]])
        """
        options = self.options()
        width = float(options['width'])
        head = options.pop('head')
        if head == 0: style = '<|-'
        elif head == 1: style = '-|>'
        elif head == 2: style = '<|-|>'
        else: raise KeyError('head parameter must be one of 0 (start), 1 (end) or 2 (both).')
        arrowsize = float(options.get('arrowsize',5))
        head_width=arrowsize
        head_length=arrowsize*2.0
        color = to_mpl_color(options['rgbcolor'])
        from matplotlib.patches import FancyArrowPatch
        from matplotlib.path import Path
        bpath = Path(self.vertices, self.codes)
        p = FancyArrowPatch(path=bpath,
                            lw=width, arrowstyle='%s,head_width=%s,head_length=%s'%(style,head_width, head_length),
                            fc=color, ec=color, linestyle=options['linestyle'])
        p.set_zorder(options['zorder'])
        subplot.add_patch(p)
        return p


class Arrow(GraphicPrimitive):
    """
    Primitive class that initializes the (line) arrow graphics type

    EXAMPLES:
    We create an arrow graphics object, then take the 0th entry
    in it to get the actual Arrow graphics primitive:
        sage: P = arrow((0,1), (2,3))[0]
        sage: type(P)
        <class 'sage.plot.arrow.Arrow'>
        sage: P
        Arrow from (0.0,1.0) to (2.0,3.0)
    """
    def __init__(self, xtail, ytail, xhead, yhead, options):
        """
        Create an arrow graphics primitive.

        EXAMPLES:
            sage: from sage.plot.arrow import Arrow
            sage: Arrow(0,0,2,3,{})
            Arrow from (0.0,0.0) to (2.0,3.0)
        """
        self.xtail = float(xtail)
        self.xhead = float(xhead)
        self.ytail = float(ytail)
        self.yhead = float(yhead)
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a bounding box for this arrow.

        EXAMPLES:
            sage: d = arrow((1,1), (5,5)).get_minmax_data()
            sage: d['xmin']
            1.0
            sage: d['xmax']
            5.0
        """
        return {'xmin': min(self.xtail, self.xhead),
                'xmax': max(self.xtail, self.xhead),
                'ymin': min(self.ytail, self.yhead),
                'ymax': max(self.ytail, self.yhead)}


    def _allowed_options(self):
        """
        Return the dictionary of allowed options for the line arrow graphics primitive.

        EXAMPLES:
             sage: from sage.plot.arrow import Arrow
             sage: list(sorted(Arrow(0,0,2,3,{})._allowed_options().iteritems()))
             [('arrowshorten', 'The length in points to shorten the arrow.'),
             ('arrowsize', 'The size of the arrowhead'),
             ('head', '2-d only: Which end of the path to draw the head (one of 0 (start), 1 (end) or 2 (both)'),
             ('hue', 'The color given as a hue.'),
             ('linestyle', "2d only: The style of the line, which is one of 'dashed', 'dotted', 'solid', 'dashdot'."),
             ('rgbcolor', 'The color as an rgb tuple.'),
             ('width', 'The width of the shaft of the arrow, in points.'),
             ('zorder', '2-d only: The layer level in which to draw')]
        """
        return {'width':'The width of the shaft of the arrow, in points.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'arrowshorten':'The length in points to shorten the arrow.',
                'arrowsize':'The size of the arrowhead',
                'zorder':'2-d only: The layer level in which to draw',
                'head':'2-d only: Which end of the path to draw the head (one of 0 (start), 1 (end) or 2 (both)',
                'linestyle':"2d only: The style of the line, which is one of 'dashed', 'dotted', 'solid', 'dashdot'."}

    def _plot3d_options(self, options=None):
        if options == None:
            options = self.options()
        options = dict(self.options())
        options_3d = {}
        if 'width' in options:
            options_3d['thickness'] = options['width']
            del options['width']
        # ignore zorder and head in 3d plotting
        if 'zorder' in options:
            del options['zorder']
        if 'head' in options:
            del options['head']
        if 'linestyle' in options:
            del options['linestyle']
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        EXAMPLE:
            sage: arrow((0,0),(1,1)).plot3d()
        """
        from sage.plot.plot3d.shapes2 import line3d
        options = self._plot3d_options()
        options.update(kwds)
        return line3d([(self.xtail, self.ytail, 0), (self.xhead, self.yhead, 0)], arrow_head=True, **options)

    def _repr_(self):
        """
        Text representation of an arrow graphics primitive.

        EXAMPLES:
            sage: from sage.plot.arrow import Arrow
            sage: Arrow(0,0,2,3,{})._repr_()
            'Arrow from (0.0,0.0) to (2.0,3.0)'
        """
        return "Arrow from (%s,%s) to (%s,%s)"%(self.xtail, self.ytail, self.xhead, self.yhead)

    def _render_on_subplot(self, subplot):
        """
        Render this arrow in a subplot.  This is the key function that
        defines how this arrow graphics primitive is rendered in
        matplotlib's library.

        EXAMPLES:
        This function implicitly ends up rendering this arrow on a matplotlib subplot:
            sage: arrow((0,1), (2,-1))
        """
        options = self.options()
        head = options.pop('head')
        if head == 0: style = '<|-'
        elif head == 1: style = '-|>'
        elif head == 2: style = '<|-|>'
        else: raise KeyError('head parameter must be one of 0 (start), 1 (end) or 2 (both).')
        width = float(options['width'])
        arrowshorten_end = float(options.get('arrowshorten',0))/2.0+width*2
        arrowsize = float(options.get('arrowsize',5))
        head_width=arrowsize
        head_length=arrowsize*2.0
        color = to_mpl_color(options['rgbcolor'])
        from matplotlib.patches import FancyArrowPatch
        p = FancyArrowPatch((self.xtail, self.ytail), (self.xhead, self.yhead),
                            lw=width, arrowstyle='%s,head_width=%s,head_length=%s'%(style,head_width, head_length),
                            shrinkA=arrowshorten_end, shrinkB=arrowshorten_end,
                            fc=color, ec=color, linestyle=options['linestyle'])
        p.set_zorder(options['zorder'])
        subplot.add_patch(p)
        return p

@rename_keyword(color='rgbcolor')
@options(width=2, rgbcolor=(0,0,1),zorder=2, head = 1, linestyle='solid')
def arrow(tailpoint=None, headpoint=None, path=None, **options):
    """
    If tailpoint and headpoint are provided, returns an arrow from (xmin, ymin)
    to (xmax, ymax).  If tailpoint or headpoint is None and path is not None,
    returns an arrow along the path.  (See further info on paths in bezier_path).

    INPUT
        tailpoint -- the starting point of the arrow
        headpoint -- where the arrow is pointing to
        path -- the list of points and control points (see bezier_path for detail) that
                the arrow will follow from source to destination
        head -- 0, 1 or 2, whether to draw the head at the start (0), end (1) or both (2)
                of the path (using 0 will swap headpoint and tailpoint).  This is ignored
                in 3-d plotting.
        width -- (default 2) the width of the arrow shaft, in points
        color -- (default (0,0,1)) the color of the arrow (as an rgb tuple or a string)
        hue -- the color of the arrow (as a number)
        arrowsize -- the size of the arrowhead
        arrowshorten -- the length in points to shorten the arrow (ignored if using path
                parameter)
        zorder -- the layer level to draw the arrow-- note that this is ignored in 3-d
                plotting.

    EXAMPLES:

    A straight, blue arrow
       sage: arrow((1, 1), (3, 3))

    Make a red arrow:
       sage: arrow((-1, -1), (2, 3), color=(1,0,0))
       sage: arrow((-1, -1), (2, 3), color='red')

    You can change the width of an arrow:
        sage: arrow((1, 1), (3, 3), width=5, arrowsize=15)

    A pretty circle of arrows:
        sage: sum([arrow((0,0), (cos(x),sin(x)), hue=x/(2*pi)) for x in [0..2*pi,step=0.1]]).show(aspect_ratio=1)

    If we want to draw the arrow between objects, for example, the
    boundaries of two lines, we can use the arrowshorten option
    to make the arrow shorter by a certain number of points.
        sage: line([(0,0),(1,0)],thickness=10)+line([(0,1),(1,1)], thickness=10)+arrow((0.5,0),(0.5,1), arrowshorten=10,rgbcolor=(1,0,0))


    """
    from sage.plot.plot import Graphics
    g = Graphics()
    if headpoint is not None and tailpoint is not None:
        xtail, ytail = tailpoint
        xhead, yhead = headpoint
        g.add_primitive(Arrow(xtail, ytail, xhead, yhead, options=options))
    elif path is not None:
        g.add_primitive(CurveArrow(path, options=options))
    else:
        raise TypeError('Arrow requires either both headpoint and tailpoint or a path parameter.')
    return g
