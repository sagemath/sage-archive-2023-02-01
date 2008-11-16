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
from sage.plot.misc import options, rename_keyword, to_mpl_color
from math import pi

class Disk(GraphicPrimitive):
    """
    Primitive class that initializes the
    disk graphics type
    """
    def __init__(self, point, r, angle, options):
        self.x = float(point[0])
        self.y = float(point[1])
        self.r = float(r)
        self.rad1 = float(angle[0])
        self.rad2 = float(angle[1])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES:
            sage: d = disk((5,5), 1, (pi/2, pi), rgbcolor=(0,0,0))
            sage: d = d.get_minmax_data()
            sage: d['xmin']
            4.0
            sage: d['ymin']
            4.0
            sage: d['xmax']
            6.0

        """
        from sage.plot.plot import minmax_data
        return minmax_data([self.x - self.r, self.x + self.r],
                           [self.y - self.r, self.y + self.y],
                           dict=True)

    def _allowed_options(self):
        return {'alpha':'How transparent the line is.',
                'fill': 'Whether or not to fill the polygon.',
                'thickness':'How thick the border of the polygon is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

    def _repr_(self):
        return "Disk defined by (%s,%s) with r=%s spanning (%s, %s) radians"%(self.x,
        self.y, self.r, self.rad1, self.rad2)

    def _render_on_subplot(self, subplot):
        import matplotlib.patches as patches
        options = self.options()
        deg1 = self.rad1*(360.0/(2.0*pi)) #convert radians to degrees
        deg2 = self.rad2*(360.0/(2.0*pi))
        p = patches.Wedge((float(self.x), float(self.y)), float(self.r), float(deg1),
                            float(deg2))
        p.set_linewidth(float(options['thickness']))
        p.set_fill(options['fill'])
        p.set_alpha(options['alpha'])
        c = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(c)
        p.set_facecolor(c)
        subplot.add_patch(p)

@options(alpha=1, fill=True, rgbcolor=(0,0,1), thickness=0)
def disk(point, radius, angle, **options):
    r"""
    A disk at a point = $(x,y)$ with radius = $r$
    spanning (in radians) angle=$(rad1, rad2)$

    Type \code{disk.options} to see all options.

    EXAMPLES:
    Make some dangerous disks:

        sage: bl = disk((0.0,0.0), 1, (pi, 3*pi/2), rgbcolor=(1,1,0))
        sage: tr = disk((0.0,0.0), 1, (0, pi/2), rgbcolor=(1,1,0))
        sage: tl = disk((0.0,0.0), 1, (pi/2, pi), rgbcolor=(0,0,0))
        sage: br = disk((0.0,0.0), 1, (3*pi/2, 2*pi), rgbcolor=(0,0,0))
        sage: P  = tl+tr+bl+br
        sage: P.show(figsize=(4,4),xmin=-2,xmax=2,ymin=-2,ymax=2)

    """
    from sage.plot.plot import Graphics
    g = Graphics()
    g.add_primitive(Disk(point, radius, angle, options))
    return g
