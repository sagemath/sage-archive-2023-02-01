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
from sage.misc.misc import verbose

class MatrixPlot(GraphicPrimitive):
    """
    Primitive class that initializes the
    matrix_plot graphics type
    """
    def __init__(self, xy_data_array, xrange, yrange, options):
        self.xrange = xrange
        self.yrange = yrange
        self.xy_data_array = xy_data_array
        self.xy_array_row = len(xy_data_array)
        self.xy_array_col = len(xy_data_array[0])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES:
            sage: m = matrix_plot(matrix([[1,3,5,1],[2,4,5,6],[1,3,5,7]]))[0]
            sage: list(sorted(m.get_minmax_data().items()))
            [('xmax', 4), ('xmin', 0), ('ymax', 3), ('ymin', 0)]

        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xrange, self.yrange, dict=True)

    def _allowed_options(self):
        return {'cmap':"""The colormap, one of (autumn, bone, cool, copper,
                gray, hot, hsv, jet, pink, prism, spring, summer, winter)"""}

    def _repr_(self):
        return "MatrixPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        options = self.options()
        given_cmap = options['cmap']
        #cm is the matplotlib color map module
        from matplotlib import cm
        from matplotlib.colors import LinearSegmentedColormap as C
        key_error = False
        try:
            cmap = cm.__dict__[given_cmap]
        except KeyError:
            key_error = True

        if key_error or not isinstance(cmap, C):
            possibilities = ', '.join([str(x) for x in cm.__dict__.keys() if
                                       isinstance(cm.__dict__[x], C)])
            verbose("The possible color maps include: %s"%possibilities, level=0)
            raise RuntimeError, "Color map %s not known"%given_cmap

        subplot.imshow(self.xy_data_array, cmap=cmap, interpolation='nearest', extent=(0,self.xrange[1],0,self.yrange[1]))


@options(cmap='gray')
def matrix_plot(mat, **options):
    r"""
    A plot of a given matrix or 2D array.

    Each ($i$th, $j$th) matrix element is given a different
    color value depending on its relative size compared
    to the other elements in the matrix.

    The tick marks drawn on the frame axes denote the
    ($i$th, $j$th) element of the matrix.

    EXAMPLES:

    A matrix over ZZ colored with different grey levels:

        sage: matrix_plot(matrix([[1,3,5,1],[2,4,5,6],[1,3,5,7]]))

    Here we make a random matrix over RR and use cmap='hsv'
    to color the matrix elements different RGB colors:

        sage: matrix_plot(random_matrix(RDF, 50), cmap='hsv')

    Another random plot, but over GF(389):
        sage: matrix_plot(random_matrix(GF(389), 10), cmap='Oranges')

    TESTS:

        sage: matrix_plot(random_matrix(RDF, 50), cmap='jolies')
        Traceback (most recent call last):
        ...
        RuntimeError: Color map jolies not known

        sage: matrix_plot(random_matrix(RDF, 50), cmap='mpl')
        Traceback (most recent call last):
        ...
        RuntimeError: Color map mpl not known
    """
    from sage.plot.plot import Graphics
    from sage.matrix.all import is_Matrix
    from matplotlib.numerix import array
    if not is_Matrix(mat) or (isinstance(mat, (list, tuple)) and isinstance(mat[0], (list, tuple))):
        raise TypeError, "mat must be of type Matrix or a two dimensional array"

    if is_Matrix(mat):
        xrange = (0, mat.ncols())
        yrange = (0, mat.nrows())
    else:
        xrange = (0, len(mat[0]))
        yrange = (0, len(mat))
    xy_data_array = [array(r, dtype=float) for r in mat]

    g = Graphics()
    g.add_primitive(MatrixPlot(xy_data_array, xrange, yrange, options))
    return g
