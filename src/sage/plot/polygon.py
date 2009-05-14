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
from sage.plot.primitive import GraphicPrimitive_xydata
from sage.plot.misc import options, rename_keyword, to_mpl_color

class Polygon(GraphicPrimitive_xydata):
    """
    Primitive class for the Polygon graphics type.  For information
    on actual plotting, please see polygon?, polygon2d?, or polygon3d?

    INPUT:

    - xdata - list of `x`-coordinates of points defining Polygon

    - ydata - list of `y`-coordinates of points defining Polygon

    - options - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via ``circle``::

        sage: from sage.plot.polygon import Polygon
        sage: P = Polygon([1,2,3],[2,3,2],{'alpha':.5})
        sage: P
        Polygon defined by 3 points
        sage: P.options()['alpha']
        0.500000000000000
        sage: P.ydata
        [2, 3, 2]

    TESTS:

    We test creating polygons::

        sage: polygon([(0,0), (1,1), (0,1)])
        sage: polygon([(0,0,1), (1,1,1), (2,0,1)])
    """
    def __init__(self, xdata, ydata, options):
        """
        Initializes base class Polygon.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,1), (-1,3)], thickness=2)
            sage: P[0].xdata
            [0.0, 1.0, -1.0]
            sage: P[0].options()['thickness']
            2
        """
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive_xydata.__init__(self, options)

    def _repr_(self):
        """
        String representation of Polygon primitive.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,1), (-1,3)])
            sage: p=P[0]; p
            Polygon defined by 3 points
        """
        return "Polygon defined by %s points"%len(self)

    def __getitem__(self, i):
        """
        Returns `i`th vertex of Polygon primitive, starting count
        from 0th vertex.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,1), (-1,3)])
            sage: p=P[0]
            sage: p[0]
            (0.0, 0.0)
        """
        return self.xdata[i], self.ydata[i]

    def __setitem__(self, i, point):
        """
        Changes `i`th vertex of Polygon primitive, starting count
        from 0th vertex.  Note that this only changes a vertex,
        but does not create new vertices.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p=P[0]
            sage: [p[i] for i in range(4)]
            [(0.0, 0.0), (1.0, 2.0), (0.0, 1.0), (-1.0, 2.0)]
            sage: p[2]=(0,.5)
            sage: p[2]
            (0.0, 0.5)
        """
        i = int(i)
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        """
        Returns number of vertices of Polygon primitive.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p=P[0]
            sage: len(p)
            4
        """
        return len(self.xdata)

    def _allowed_options(self):
        """
        Return the allowed options for the Polygon class.

        EXAMPLES::

            sage: P = polygon([(1,1), (1,2), (2,2), (2,1)], alpha=.5)
            sage: P[0]._allowed_options()['alpha']
            'How transparent the figure is.'
        """
        return {'alpha':'How transparent the figure is.',
                'thickness': 'How thick the border line is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'zorder':'The layer level in which to draw'}

    def _plot3d_options(self, options=None):
        """
        Translate 2d plot options into 3d plot options.

        EXAMPLES::

            sage: P = polygon([(1,1), (1,2), (2,2), (2,1)], alpha=.5)
            sage: p=P[0]; p
            Polygon defined by 4 points
            sage: q=p.plot3d()
            sage: q.texture.opacity
            0.500000000000000
        """
        if options == None:
            options = dict(self.options())
        if 'thickness' in options:
            del options['thickness']
        if 'zorder' in options:
            del options['zorder']
        return GraphicPrimitive_xydata._plot3d_options(self, options)

    def plot3d(self, z=0, **kwds):
        """
        Plots a 2D polygon in 3D, with default height zero.

        INPUT:


        -  ``z`` - optional 3D height above `xy`-plane, or a list of
           heights corresponding to the list of 2D polygon points.

        EXAMPLES:

        A pentagon::

            sage: polygon([(cos(t), sin(t)) for t in srange(0, 2*pi, 2*pi/5)]).plot3d()

        Showing behavior of the optional parameter z::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p = P[0]; p
            Polygon defined by 4 points
            sage: q = p.plot3d()
            sage: q.obj_repr(q.testing_render_params())[2]
            ['v 0 0 0', 'v 1 2 0', 'v 0 1 0', 'v -1 2 0']
            sage: r = p.plot3d(z=3)
            sage: r.obj_repr(r.testing_render_params())[2]
            ['v 0 0 3', 'v 1 2 3', 'v 0 1 3', 'v -1 2 3']
            sage: s = p.plot3d(z=[0,1,2,3])
            sage: s.obj_repr(s.testing_render_params())[2]
            ['v 0 0 0', 'v 1 2 1', 'v 0 1 2', 'v -1 2 3']

        TESTS:

        Heights passed as a list should have same length as
        number of points::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p = P[0]
            sage: q = p.plot3d(z=[2,-2])
            Traceback (most recent call last):
            ...
            ValueError: Incorrect number of heights given
        """
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        options = self._plot3d_options()
        options.update(kwds)
        zdata=[]
        if type(z) is list:
            zdata=z
        else:
            zdata=[z]*len(self.xdata)
        if len(zdata)==len(self.xdata):
            return IndexFaceSet([[(x, y, z) for x, y, z in zip(self.xdata, self.ydata, zdata)]], **options)
        else:
            raise ValueError, 'Incorrect number of heights given'

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
        """
        import matplotlib.patches as patches
        options = self.options()
        p = patches.Polygon([(self.xdata[i],self.ydata[i]) for i in xrange(len(self.xdata))])
        p.set_linewidth(float(options['thickness']))
        a = float(options['alpha'])
        p.set_alpha(a)
        c = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(c)
        p.set_facecolor(c)
        subplot.add_patch(p)

def polygon(points, **options):
    """
    Returns either a 2-dimensional or 3-dimensional polygon depending
    on value of points.

    For information regarding additional arguments, see either polygon2d?
    or polygon3d?.

    EXAMPLES::

        sage: polygon([(0,0), (1,1), (0,1)])
        sage: polygon([(0,0,1), (1,1,1), (2,0,1)])
    """
    try:
        return polygon2d(points, **options)
    except ValueError:
        from sage.plot.plot3d.shapes2 import polygon3d
        return polygon3d(points, **options)

@options(alpha=1, rgbcolor=(0,0,1), thickness=0)
def polygon2d(points, **options):
    r"""
    Returns a polygon defined by ``points``.

    Type ``polygon.options`` for a dictionary of the default
    options for polygons.  You can change this to change
    the defaults for all future polygons.  Use ``polygon.reset()``
    to reset to the default options.

    EXAMPLES:

    We create a purple-ish polygon::

        sage: polygon2d([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))

    Some modern art -- a random polygon::

        sage: v = [(randrange(-5,5), randrange(-5,5)) for _ in range(10)]
        sage: polygon2d(v)

    A purple hexagon::

        sage: L = [[cos(pi*i/3),sin(pi*i/3)] for i in range(6)]
        sage: polygon2d(L, rgbcolor=(1,0,1))

    A green deltoid::

        sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8,3/4,1/2))

    A blue hypotrochoid::

        sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8,1/4,1/2))

    Another one::

        sage: n = 4; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8,1/4,3/4))

    A purple epicycloid::

        sage: m = 9; b = 1
        sage: L = [[m*cos(pi*i/100)+b*cos((m/b)*pi*i/100),m*sin(pi*i/100)-b*sin((m/b)*pi*i/100)] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(7/8,1/4,3/4))

    A brown astroid::

        sage: L = [[cos(pi*i/100)^3,sin(pi*i/100)^3] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(3/4,1/4,1/4))

    And, my favorite, a greenish blob::

        sage: L = [[cos(pi*i/100)*(1+cos(pi*i/50)), sin(pi*i/100)*(1+sin(pi*i/50))] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8, 3/4, 1/2))

    This one is for my wife::

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,100)]
        sage: polygon2d(L, rgbcolor=(1,1/4,1/2))

    AUTHORS:

    - David Joyner (2006-04-14): the long list of examples above.

    """
    from sage.plot.plot import xydata_from_point_list, Graphics
    xdata, ydata = xydata_from_point_list(points)
    g = Graphics()
    g.add_primitive(Polygon(xdata, ydata, options))
    return g
