"""
Constructor for making matrices
"""


#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import sage.rings.all as rings
import sage.matrix.matrix_space as matrix_space
from sage.structure.sequence import Sequence

def matrix(arg0=None, arg1=None, arg2=None, arg3=None, sparse=False):
    """
    Create a matrix.

    INPUT:
    Supported formats
        1. matrix([sparse]):
               the 0x0 matrix over ZZ
        2. matrix(list_of_rows, [sparse]):
               matrix with each row constructed from the list_of_rows
        3. matrix(nrows, entries, [sparse]):
               matrix with given number of rows and flat list of entries
        4. matrix(nrows, ncols, entries, [sparse]):
               matrix with given number of rows and columns with flat list of entries
        5. matrix(ring, list_of_row, [sparse]):
               matrix over ring with rows the elements of the list_of_rows
        6. matrix(ring, nrows, entries, [sparse]):
               matrix over ring with given number of rows and entries from the flat list
        7. matrix(ring, nrows, ncols, entries, [sparse]):
               matrix over the given ring with given number of rows and columns and entries.

    In each case the sparse option is optional, must be explicitly named (i.e., sparse=True),
    and may be either True or false.

    OUTPUT:
        a matrix

    EXAMPLES:
    1. The degenerate matrix input format:
        sage: m = matrix(); m
        []
        sage: parent(m)
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m = matrix(sparse=True); m
        []
        sage: parent(m)
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring

    2. The "matrix(list_of_rows, [sparse])" input format.
       Notice that SAGE is careful to find a sensible common
       ring for all the entries (using the Sequence object):

        sage: m = matrix([[1/3,2+x],[3,4]]); m
        [  1/3 x + 2]
        [    3     4]
        sage: parent(m)
        Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

        sage: m = matrix([[1,2],[3,4/7]]); m
        [  1   2]
        [  3 4/7]
        sage: parent(m)
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        sage: m = matrix([[1,2,3], [RDF(2), CDF(1,2), 3]]); m
        [        1.0         2.0         3.0]
        [        2.0 1.0 + 2.0*I         3.0]
        sage: parent(m)
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field

        sage: V = GF(7)^2
        sage: m = matrix([V.0, V.0+2*V.1], sparse=True); m
        [1 0]
        [1 2]
        sage: parent(m)
        Full MatrixSpace of 2 by 2 sparse matrices over Finite Field of size 7

    3. matrix(nrows, entries, [sparse]):
        sage: matrix(2,[1,2,3, 4,5,6])
        [1 2 3]
        [4 5 6]
        sage: matrix(3, [1,2,  3/4, 5/6, 7*x, 8*x])
        [  1   2]
        [3/4 5/6]
        [7*x 8*x]

    The number of rows must divide the number of entries.
        sage: matrix(5,[1,2,3, 4,5,6])
        Traceback (most recent call last):
        ...
        TypeError: entries has the wrong length

    4. matrix(nrows, ncols, entries, [sparse]):
        sage: matrix(2,3, [1,2,3, 4,5,6])
        [1 2 3]
        [4 5 6]
        sage: matrix(3,2, [1,2,  3/4, 5/6, 7*x, 8*x])
        [  1   2]
        [3/4 5/6]
        [7*x 8*x]

    The length of the entry list must be the same as the number of rows times columns.
        sage: matrix(3,3, [1,2,  3/4, 5/6, 7*x, 8*x])
        Traceback (most recent call last):
        ...
        TypeError: entries has the wrong length

    5. matrix(ring, list_of_row, [sparse]):
        sage: m = matrix(QQ, [[1,2], [3,4]]); m
        [1 2]
        [3 4]
        sage: parent(m)
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: m = matrix(GF(7), [[1/3,2/3,1/2], [3/4,4/5,7]]); m
        [5 3 4]
        [6 5 0]
        sage: parent(m)
        Full MatrixSpace of 2 by 3 dense matrices over Finite Field of size 7

    6. matrix(ring, nrows, entries, [sparse]):
        sage: M = ZZ^4
        sage: m = matrix(QQ, [M([1,2,3,4]), M([-1,0,3,4])]); m
        [ 1  2  3  4]
        [-1  0  3  4]
        sage: parent(m)
        Full MatrixSpace of 2 by 4 dense matrices over Rational Field

    7. matrix(ring, nrows, ncols, entries, [sparse])
        sage: matrix(QQ, 2, 2, [1,2,3,4])
        [1 2]
        [3 4]

        sage: matrix(GF(5), 2, 3, range(6))
        [0 1 2]
        [3 4 0]

        sage: matrix(ZZ, 10, 10, range(100)).parent()
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring

        sage: matrix(ZZ, 10, 10, range(100), sparse=True).parent()
        Full MatrixSpace of 10 by 10 sparse matrices over Integer Ring

        sage: R = PolynomialRing(QQ, 9, 'x')
        sage: A = matrix(R, 3, 3, R.gens()); A
        [x0 x1 x2]
        [x3 x4 x5]
        [x6 x7 x8]
        sage: det(A)   # just for fun

    """
    if arg0 is None:
        # 1. matrix([sparse]):
        #       the 0x0 matrix over ZZ
        # the degenerate matrix
        ring = rings.ZZ
        nrows = 0
        ncols = 0
        entries = []

    elif isinstance(arg0, (list, tuple)):
        # 2. matrix(list_of_rows, [sparse]):
        #       matrix with each row constructed from the list_of_rows
        # list of rows
        if not (arg1 is None and arg2 is None and arg3 is None):
            raise TypeError, "invalid input"
        nrows = len(arg0)
        if nrows == 0:
            ncols = 0
        else:
            ncols = len(arg0[0])
        w = sum([list(v) for v in arg0], [])
        entries, ring = prepare(w)

    elif not rings.is_Ring(arg0) and isinstance(arg1, (list, tuple)) and arg2 is None and arg3 is None:
        # 3. matrix(nrows, entries, [sparse]):
        #       matrix with given number of rows and flat list of entries
        nrows = int(arg0)
        entries, ring = prepare(arg1)
        ncols = len(entries) // nrows

    elif not rings.is_Ring(arg0) and isinstance(arg2, (list, tuple)) and arg3 is None:
        # 4. matrix(nrows, ncols, entries, [sparse]):
        #       matrix with given number of rows and columns with flat list of entries
        nrows = int(arg0)
        ncols = int(arg1)
        entries, ring = prepare(arg2)

    elif rings.is_Ring(arg0) and isinstance(arg1, (list, tuple)) and arg2 is None and arg3 is None:
        # 5. matrix(ring, list_of_row, [sparse]):
        #       matrix over ring with rows the elements of the list_of_rows
        ring = arg0
        nrows = len(arg1)
        if nrows == 0:
            ncols = 0
        else:
            try:
                ncols = len(arg1[0])
            except TypeError:
                raise TypeError, "If making a matrix with the matrix(ring, list_of_row, [sparse]) constructor, the second input must be a list of rows."
        entries = sum([list(v) for v in arg1], [])

    elif rings.is_Ring(arg0) and isinstance(arg2, (list, tuple)) and arg3 is None:
        # 6. matrix(ring, nrows, entries, [sparse]):
        #       matrix over ring with given number of rows and entries from the flat list
        ring = arg0
        nrows = int(arg1)
        entries = arg2
        ncols = len(entries) // nrows

    elif rings.is_Ring(arg0):
        # 7. matrix(ring, nrows, ncols, entries, [sparse]):
        #       matrix over the given ring with given number of rows and columns and entries.
        ring = arg0
        nrows = int(arg1)
        ncols = int(arg2)
        entries = arg3

    else:

        raise TypeError, "unknown matrix constructor format.  Type matrix? for help"

    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse=sparse)(entries)


def prepare(w):
    entries = Sequence(w)
    ring = entries.universe()
    if ring is int or ring is long:
        ring = rings.ZZ
    elif ring is float:
        ring = rings.RDF
    elif ring is complex:
        ring = rings.CDF
    elif not rings.is_Ring(ring):
        raise TypeError, "unable to find a common ring for all elements"
    return entries, ring

Matrix = matrix
