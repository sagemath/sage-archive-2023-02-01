"""
Matrix Plots
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
from sage.plot.colors import to_mpl_color, get_cmap

class MatrixPlot(GraphicPrimitive):
    """
    Primitive class for the matrix plot graphics type.  See
    ``matrix_plot?`` for help actually doing matrix plots.

    INPUT:

    - ``xy_data_array`` - list of lists giving matrix values corresponding to
      the grid

    - ``xrange`` - tuple of 2 floats indicating range for horizontal direction
      (number of columns in the matrix)

    - ``yrange`` - tuple of 2 floats indicating range for vertical direction
      (number of rows in the matrix)

    - ``options`` - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via ``matrix_plot``::

        sage: from sage.plot.matrix_plot import MatrixPlot
        sage: M = MatrixPlot([[1,3],[2,4]],(1,2),(2,3),options={'cmap':'winter'})
        sage: M
        MatrixPlot defined by a 2 x 2 data grid
        sage: M.yrange
        (2, 3)
        sage: M.xy_data_array
        [[1, 3], [2, 4]]
        sage: M.options()
        {'cmap': 'winter'}

    TESTS:

    We test creating a matrix plot::

        sage: M = matrix_plot([[mod(i,5)^j for i in range(5)] for j in range(1,6)])
    """
    def __init__(self, xy_data_array, xrange, yrange, options):
        """
        Initializes base class MatrixPlot.

        EXAMPLES::

            sage: M = matrix_plot([[mod(i,5)^j for i in range(5)] for j in range(1,6)], cmap='jet')
            sage: M[0].xrange
            (0, 5)
            sage: M[0].options()['cmap']
            'jet'
            sage: M[0].xy_array_row
            5
        """
        self.xrange = xrange
        self.yrange = yrange
        self.xy_data_array = xy_data_array
        self.xy_array_row = len(xy_data_array)
        self.xy_array_col = len(xy_data_array[0])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: m = matrix_plot(matrix([[1,3,5,1],[2,4,5,6],[1,3,5,7]]))[0]
            sage: list(sorted(m.get_minmax_data().items()))
            [('xmax', 4), ('xmin', 0), ('ymax', 3), ('ymin', 0)]

        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xrange, self.yrange, dict=True)

    def _allowed_options(self):
        """
        Return the allowed options for the MatrixPlot class.

        EXAMPLES::

            sage: M = matrix_plot([[sin(i*j) for i in range(5)] for j in range(5)])
            sage: isinstance(M[0]._allowed_options(),dict)
            True
        """
        return {'cmap':"""the name of a predefined colormap,
                        a list of colors, or an instance of a
                        matplotlib Colormap. Type: import matplotlib.cm; matplotlib.cm.datad.keys()
                        for available colormap names.""",
                'zorder':"The layer level in which to draw"}

    def _repr_(self):
        """
        String representation of MatrixPlot primitive.

        EXAMPLES::

            sage: M = matrix_plot([[sin(i*j) for i in range(5)] for j in range(5)])
            sage: m = M[0]; m
            MatrixPlot defined by a 5 x 5 data grid
        """
        return "MatrixPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: matrix_plot(random_matrix(RDF, 50), cmap='jet')
        """
        options = self.options()
        cmap = get_cmap(options['cmap'])

        subplot.imshow(self.xy_data_array, cmap=cmap, interpolation='nearest', extent=(0,self.xrange[1],0,self.yrange[1]))


@options(cmap='gray')
def matrix_plot(mat, **options):
    r"""
    A plot of a given matrix or 2D array.

    Each (`i`-th, `j`-th) matrix element is given a different
    color value depending on its relative size compared
    to the other elements in the matrix.

    The tick marks drawn on the frame axes denote the
    (`i`-th, `j`-th) element of the matrix.

    INPUT:

    - ``mat`` - a 2D matrix or array

    The following input must all be passed in as named parameters, if
    default not used:

    - ``cmap`` - a colormap (default: 'gray'), the name of
      a predefined colormap, a list of colors,
      or an instance of a matplotlib Colormap.
      Type: ``import matplotlib.cm; matplotlib.cm.datad.keys()``
      for available colormap names.

    EXAMPLES:

    A matrix over `\ZZ` colored with different grey levels::

        sage: matrix_plot(matrix([[1,3,5,1],[2,4,5,6],[1,3,5,7]]))

    Here we make a random matrix over `\RR` and use ``cmap='hsv'``
    to color the matrix elements different RGB colors::

        sage: matrix_plot(random_matrix(RDF, 50), cmap='hsv')

    Another random plot, but over `\GF{389}`::

        sage: m = random_matrix(GF(389), 10)
        sage: matrix_plot(m, cmap='Oranges')

    It also works if you lift it to the polynomial ring::

        sage: matrix_plot(m.change_ring(GF(389)['x']), cmap='Oranges')

    Here we plot a random sparse matrix::

        sage: sparse = matrix(dict([((randint(0, 10), randint(0, 10)), 1) for i in xrange(100)]))
        sage: matrix_plot(sparse)

    Plotting lists of lists also works::

        sage: matrix_plot([[1,3,5,1],[2,4,5,6],[1,3,5,7]])

    As does plotting of numpy arrays::

        sage: import numpy
        sage: matrix_plot(numpy.random.rand(10, 10))

    TESTS::

        sage: P.<t> = RR[]
        sage: matrix_plot(random_matrix(P, 3, 3))
        Traceback (most recent call last):
        ...
        TypeError: cannot coerce nonconstant polynomial to float

    ::

        sage: matrix_plot([1,2,3])
        Traceback (most recent call last):
        ...
        TypeError: mat must be a Matrix or a two dimensional array

    ::

        sage: matrix_plot([[sin(x), cos(x)], [1, 0]])
        Traceback (most recent call last):
        ...
        ValueError: can not convert array entries to floating point numbers
    """
    import numpy as np
    from sage.plot.plot import Graphics
    from sage.matrix.all import is_Matrix
    from sage.rings.all import RDF

    if is_Matrix(mat):
        mat = mat.change_ring(RDF).numpy()

    try:
        xy_data_array = np.asarray(mat, dtype = float)
    except TypeError:
        raise TypeError, "mat must be a Matrix or a two dimensional array"
    except ValueError:
        raise ValueError, "can not convert array entries to floating point numbers"

    if len(xy_data_array.shape) < 2:
        raise TypeError, "mat must be a Matrix or a two dimensional array"

    xrange = (0, xy_data_array.shape[1])
    yrange = (0, xy_data_array.shape[0])

    g = Graphics()
    g.add_primitive(MatrixPlot(xy_data_array, xrange, yrange, options))
    return g
