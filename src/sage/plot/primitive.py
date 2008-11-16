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
from sage.structure.sage_object import SageObject
from sage.misc.misc import verbose

class GraphicPrimitive(SageObject):
    """
    Base class for graphics primitives, e.g., things that knows how to draw
    themselves in 2d.

    EXAMPLES:
    We create an object that derives from GraphicPrimitive:
        sage: P = line([(-1,-2), (3,5)])
        sage: P[0]
        Line defined by 2 points
        sage: type(P[0])
        <class 'sage.plot.line.Line'>
    """
    def __init__(self, options):
        """
        Create a base class GraphicsPrimitive.  All this does is
        set the options.

        EXAMPLES:
        We indirectly test this function.
            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({})
            Graphics primitive
        """

        self.__options = options

    def _allowed_options(self):
        """
        Return the allowed options for a graphics primitive.

        OUTPUT:
            -- a reference to a dictionary.

        EXAMPLES:
            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({})._allowed_options()
            {}
        """
        return {}

    def plot3d(self, **kwds):
        raise NotImplementedError, "3d plotting not implemented for %s" % type(self)

    def _plot3d_options(self, options=None):
        """
        Translate 2d plot options into 3d plot options.
        """
        if options == None:
            options = self.options()
        options_3d = {}
        if 'rgbcolor' in options:
            options_3d['rgbcolor'] = options['rgbcolor']
            del options['rgbcolor']
        if 'alpha' in options:
            options_3d['opacity'] = options['alpha']
            del options['alpha']
        if len(options) != 0:
            raise NotImplementedError, "Unknown plot3d equivalent for %s" % ", ".join(options.keys())
        return options_3d

    def options(self):
        """
        Return the dictionary of options for this graphics primitive.

        By default this function verifies that the options are all
        valid; if any aren't a verbose message is printed with level 0.

        EXAMPLES:
            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({}).options()
            {}
        """
        from sage.plot.plot import do_verify, hue
        O = dict(self.__options)
        if do_verify:
            A = self._allowed_options()
            t = False
            K = A.keys() + ['xmin', 'xmax', 'ymin', 'ymax', 'axes']
            for k in O.keys():
                if not k in K:
                    do_verify = False
                    verbose("WARNING: Ignoring option '%s'=%s"%(k,O[k]), level=0)
                    t = True
            if t:
                s = "\nThe allowed options for %s are:\n"%self
                K.sort()
                for k in K:
                    if A.has_key(k):
                        s += "    %-15s%-60s\n"%(k,A[k])
                verbose(s, level=0)


        if 'hue' in O:
            t = O['hue']
            if not isinstance(t, (tuple,list)):
                t = [t,1,1]
            O['rgbcolor'] = hue(*t)
            del O['hue']
        return O

    def _repr_(self):
        """
        String representation of this graphics primitive.

        EXAMPLES:
            sage: from sage.plot.primitive import GraphicPrimitive
            sage: GraphicPrimitive({})._repr_()
            'Graphics primitive'
        """
        return "Graphics primitive"



class GraphicPrimitive_xydata(GraphicPrimitive):
    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES:
            sage: d = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))[0].get_minmax_data()
            sage: d['ymin']
            0.0
            sage: d['xmin']
            1.0

            sage: d = point((3, 3), rgbcolor=hue(0.75))[0].get_minmax_data()
            sage: d['xmin']
            3.0
            sage: d['ymin']
            3.0

            sage: l = line([(100, 100), (120, 120)])[0]
            sage: d = l.get_minmax_data()
            sage: d['xmin']
            100.0
            sage: d['xmax']
            120.0

        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xdata, self.ydata, dict=True)

