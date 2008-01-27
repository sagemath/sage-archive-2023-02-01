"""
Matrix Constructor.
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
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.integer_ring import ZZ
from sage.misc.misc_c import running_total
from matrix import is_Matrix

import sage.categories.pushout

def matrix(arg0=None, arg1=None, arg2=None, arg3=None, sparse=None):
    """
    Create a matrix.

    INPUT:
    Supported formats
        1.  matrix([sparse=True]):
                the 0x0 matrix over ZZ
        2.  matrix(list_of_rows, [sparse=True]):
                matrix with each row constructed from the list_of_rows
        3.  matrix(nrows, entries, [sparse=True]):
                matrix with given number of rows and flat list of entries
        4.  matrix(nrows, ncols, entries, [sparse=True]):
                matrix with given number of rows and columns with flat list of entries
        5.  matrix(ring, list_of_row, [sparse=True]):
                matrix over ring with rows the elements of the list_of_rows
        6.  matrix(ring, nrows, entries, [sparse=True]):
                matrix over ring with given number of rows and entries from the flat list
        7.  matrix(ring, nrows, ncols, entries, [sparse=True]):
                matrix over the given ring with given number of rows and columns and entries.
        8.  matrix(numpy_array)
        9.  matrix(object)
        10. matrix(ring, object)
        11. matrix(object, ring)

    The sparse option is optional, must be explicitly named (i.e.,
    sparse=True), and may be either True or False.

    The entries can instead be a dictionary of key:value pairs of the
    form (i,j):x, where i,j are integers instead of a list.  If sparse
    is not specified and the entries are a dictionary, it default to
    True.

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

        sage: x = polygen(QQ)
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
        sage: det(A)
        -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8

    CREATING SPARSE MATRICES FROM DICTS:
        sage: a = matrix({(1,2):10, (2,3):5/1})
        sage: print a, a.parent()
        [ 0  0  0  0]
        [ 0  0 10  0]
        [ 0  0  0  5] Full MatrixSpace of 3 by 4 sparse matrices over Rational Field
        sage: a = matrix({(1,2):10})
        sage: print a, a.is_sparse()
        [ 0  0  0]
        [ 0  0 10] True
        sage: a = matrix(3,{(1,2):10})
        sage: print a, a.is_sparse()
        [ 0  0  0]
        [ 0  0 10]
        [ 0  0  0] True
        sage: a = matrix(3,5,{(1,2):10})
        sage: print a, a.is_sparse()
        [ 0  0  0  0  0]
        [ 0  0 10  0  0]
        [ 0  0  0  0  0] True
        sage: a = matrix(QQ, 3, {(1,2):10})
        sage: print a, a.is_sparse()
        [ 0  0  0]
        [ 0  0 10]
        [ 0  0  0] True
        sage: a = matrix(QQ, 3, {(1,2):10}, sparse=True)
        sage: print a, a.is_sparse()
        [ 0  0  0]
        [ 0  0 10]
        [ 0  0  0] True
        sage: a = matrix(QQ, 3, 5, {(1,2):10}, sparse=True)
        sage: print a, a.is_sparse()
        [ 0  0  0  0  0]
        [ 0  0 10  0  0]
        [ 0  0  0  0  0] True

    8. Creating a matrix from a numpy array
        Any numpy array with a datatype of float or complex may be passed to matrix
        If the data type is float the result will be a matrix over the real double field.
        complex numpy arrays will give matrices over the complex double field.
        The data of the numpy array must be contiguous, so slices of other matrices will raise an exception.
        sage: import numpy
        sage: n=numpy.array([[1,2],[3,4]],float)
        sage: m=matrix(n)
        sage: n=numpy.array([[numpy.complex(0,1),numpy.complex(0,2)],[3,4]],complex)
        sage: m=matrix(n)
        sage: a = numpy.array([[1,2],[3,4]],'int32')
        sage: matrix(a)
        [1 2]
        [3 4]

    9.  Creating a matrix from Sage objects.
        sage: v = vector(ZZ, [1, 10, 100])
        sage: matrix(v)
        [  1  10 100]

        sage: matrix(GF(7), v)
        [1 3 2]

        sage: matrix(ZZ['x'], v)
        [  1  10 100]

    TESTS:
        sage: import numpy
        sage: a = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float32')
        sage: matrix(a)
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        sage: a = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float64')
        sage: matrix(a)
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        sage: a = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'complex64')
        sage: matrix(a)
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        sage: a = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'complex128')
        sage: matrix(a)
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]


        sage: a = matrix([[1,2],[3,4]])
        sage: b = matrix(a.numpy()); b
        [1 2]
        [3 4]
        sage: a == b
        True
        sage: c = matrix(a.numpy('float32')); c
        [1.0 2.0]
        [3.0 4.0]

    """
    if hasattr(arg0, '_matrix_'):
        if arg1 is None:
            try:
                arg1 = rings.ZZ
                return arg0._matrix_(arg1)
            except TypeError:
                return arg0._matrix_()
        return arg0._matrix_(arg1)

    if hasattr(arg1, '_matrix_'):
        return arg1._matrix_(arg0)

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

    elif isinstance(arg0, dict):
        # 2. matrix(dict_of_rows, sparse=True):
        if sparse is None: sparse = True
        entries, ring = prepare_dict(arg0)
        nrows = nrows_from_dict(entries)
        ncols = ncols_from_dict(entries)

    elif not rings.is_Ring(arg0) and isinstance(arg1, (list, tuple)) and arg2 is None and arg3 is None:
        # 3. matrix(nrows, entries, [sparse]):
        #       matrix with given number of rows and flat list of entries
        nrows = int(arg0)
        entries, ring = prepare(arg1)
        ncols = len(entries) // nrows

    elif not rings.is_Ring(arg0) and isinstance(arg1, dict) and arg2 is None and arg3 is None:
        # 3. matrix(nrows, entries, sparse=True):
        #       matrix with given number of rows and sparse dict of entries
        if sparse is None: sparse = True
        nrows = int(arg0)
        entries, ring = prepare_dict(arg1)
        ncols = ncols_from_dict(entries)

    elif not rings.is_Ring(arg0) and isinstance(arg2, (list, tuple)) and arg3 is None:
        # 4. matrix(nrows, ncols, entries, [sparse]):
        #       matrix with given number of rows and columns with flat list of entries
        nrows = int(arg0)
        ncols = int(arg1)
        entries, ring = prepare(arg2)

    elif not rings.is_Ring(arg0) and isinstance(arg2, dict) and arg3 is None:
        # 4. matrix(nrows, ncols, entries, sparse=True):
        #       matrix with given number of rows and columns with flat list of entries
        if sparse is None: sparse = True
        nrows = int(arg0)
        ncols = int(arg1)
        entries, ring = prepare_dict(arg2)

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

    elif rings.is_Ring(arg0) and isinstance(arg1, dict) and arg2 is None and arg3 is None:
        # 5. matrix(ring, dict, sparse=True):
        #       matrix over ring with rows the elements the dict
        if sparse is None: sparse = True
        ring = arg0
        entries = arg1
        nrows = nrows_from_dict(entries)

    elif rings.is_Ring(arg0) and isinstance(arg2, (list, tuple)) and arg3 is None:
        # 6. matrix(ring, nrows, entries, [sparse]):
        #       matrix over ring with given number of rows and entries from the flat list
        ring = arg0
        nrows = int(arg1)
        entries = arg2
        if nrows == 0:
            ncols = 0
        else:
            ncols = len(entries) // nrows

    elif rings.is_Ring(arg0) and isinstance(arg2, dict) and arg3 is None:
        # 6. matrix(ring, nrows, entries, sparse=True)
        #       matrix over ring with given number of rows and entries from the dict
        if sparse is None: sparse = True
        ring = arg0
        nrows = int(arg1)
        entries = arg2
        if nrows == 0:
            ncols = 0
        else:
            ncols = ncols_from_dict(entries)

    elif rings.is_Ring(arg0):
        # 7. matrix(ring, nrows, ncols, entries, [sparse]):
        #       matrix over the given ring with given number of rows and columns and entries.
        ring = arg0
        nrows = int(arg1)
        if arg2 is None:
            ncols = nrows
        else:
            ncols = int(arg2)
        entries = arg3
        if isinstance(entries, dict):
            if sparse is None: sparse = True

    else:
        import numpy
        if isinstance(arg0,numpy.ndarray):
            str_dtype = str(arg0.dtype)

            if not( arg0.flags.c_contiguous is True or arg0.flags.f_contiguous is True):
                raise TypeError('numpy matrix must be either c_contiguous or f_contiguous')
            if str_dtype.count('float32')==1:
                m=matrix(RDF,arg0.shape[0],arg0.shape[1],0)
                m._replace_self_with_numpy32(arg0)

            elif str_dtype.count('float64')==1:
                m=matrix(RDF,arg0.shape[0],arg0.shape[1],0)
                m._replace_self_with_numpy(arg0)

            elif str_dtype.count('complex64')==1:
                m=matrix(CDF,arg0.shape[0],arg0.shape[1],0)
                m._replace_self_with_numpy32(arg0)

            elif str_dtype.count('complex128')==1:
                m=matrix(CDF,arg0.shape[0],arg0.shape[1],0)
                m._replace_self_with_numpy(arg0)

            elif str_dtype.count('int') == 1:
                m = matrix(ZZ, map(list, list(arg0)))

            elif str_dtype.count('object') == 1:
                #Get the raw nested list from the numpy array
                #and feed it back into matrix
                try:
                    return matrix( map(list, list(arg0)) )
                except TypeError:
                    raise TypeError("cannot convert numpy matrix to SAGE matrix")
            else:
                raise TypeError("cannot convert numpy matrix to SAGE matrix")

            if arg0.flags.c_contiguous:
                return m
            else:
                return m.transpose()

        else:
            raise TypeError, "unknown matrix constructor format.  Type matrix? for help"

    if sparse is None:
        sparse = False

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

def prepare_dict(w):
    Z = w.items()
    X = [x for _, x in Z]
    entries, ring = prepare(X)
    return dict([(Z[i][0],entries[i]) for i in xrange(len(entries))]), ring

def nrows_from_dict(d):
    return max([0] + [ij[0] for ij in d.keys()]) + 1

def ncols_from_dict(d):
    return max([0] + [ij[1] for ij in d.keys()]) + 1

Matrix = matrix


def random_matrix(R, nrows, ncols=None, sparse=False, density=1, *args, **kwds):
    """
    Return a random matrix with entries in the ring R.

    INPUT:
        R -- a ring
        nrows -- integer; number of rows
        ncols -- (default: None); number of columns; if None defaults to nrows
        sparse -- (default; False); whether or not matrix is sparse.
        density -- integer (default: 1)
        *args, **kwds -- passed on to randomize function

    EXAMPLES:
        sage: A = random_matrix(ZZ,50,x=2^16)    # entries are up to 2^16 i size
        sage: A
        50 x 50 dense matrix over Integer Ring
    """
    if ncols is None:
        ncols = nrows
    A = matrix_space.MatrixSpace(R, nrows, ncols, sparse=sparse).zero_matrix()
    A.randomize(density=density, *args, **kwds)
    return A

def diagonal_matrix(arg0=None, arg1=None, arg2=None, sparse=None):
    """
    INPUT:
    Supported formats
        1. matrix(diagonal_entries, [sparse=True]):
               matrix with each row constructed from the list_of_rows
        2. matrix(nrows, diagonal_entries, [sparse=True]):
               matrix with each row constructed from the list_of_rows
        3. matrix(ring, diagonal_entries, [sparse=True]):
               matrix with each row constructed from the list_of_rows
        4. matrix(ring, nrows, diagonal_entries, [sparse=True]):
               matrix with given number of rows and flat list of entries
    The sparse option is optional, must be explicitly named (i.e.,
    sparse=True), and may be either True or False.

    EXAMPLES:
    Input format 1.
        sage: diagonal_matrix([1,2,3])
        [1 0 0]
        [0 2 0]
        [0 0 3]

    Input format 2.
        sage: diagonal_matrix(GF(3), [1,2,3])
        [1 0 0]
        [0 2 0]
        [0 0 0]

    Input format 3:
        sage: diagonal_matrix(3, [1,2])
        [1 0 0]
        [0 2 0]
        [0 0 0]

    Input format 4:
        sage: diagonal_matrix(GF(3), 3, [8,2])
        [2 0 0]
        [0 2 0]
        [0 0 0]
    """
    ring = None
    if isinstance(arg0, (list, tuple)):
        # Format 1
        v = arg0
        nrows = len(v)
    elif isinstance(arg0, (int, long, rings.Integer)):
        nrows = arg0
        v = arg1
    elif rings.is_Ring(arg0):
        ring = arg0
        if isinstance(arg1, (list, tuple)):
            v = arg1
            nrows = len(v)
        else:
            nrows = arg1
            v = arg2

    if isinstance(v, list):
        w = {}
        for i in range(len(v)):
            w[(i,i)] = v[i]
    else:
        w = v

    if ring is None:
        return matrix(nrows, nrows, w, sparse=sparse)
    else:
        return matrix(ring, nrows, nrows, w, sparse=sparse)


def identity_matrix(ring, n=0):
    """
    Return the n x n identity matrix over the given ring.

    EXAMPLES:
        sage: identity_matrix(QQ, 2)
        [1 0]
        [0 1]
        sage: identity_matrix(2)
        [1 0]
        [0 1]
    """
    if isinstance(ring, (int, long, rings.Integer)):
        n = ring
        ring = rings.ZZ
    return matrix_space.MatrixSpace(ring, n, n).identity_matrix()


def block_matrix(sub_matrices, nrows=None, ncols=None, subdivide=True):
    """
    Returns a larger matrix made by concatinating the sub_matrices
    (rows first, then columns). For example, the matrix

        [ A B ]
        [ C D ]

    is made up of submatrices A, B, C, and D.

    INPUT:
        sub_matrices -- matrices (must be of the correct size, or constants)
        nrows         -- (optional) the number of block rows
        ncols         -- (optional) the number of block cols
        subdivide    -- boolean, whether or not to add
                        subdivision information to the matrix

    EXAMPLES:
        sage: A = matrix(QQ, 2, 2, [3,9,6,10])
        sage: block_matrix([A, -A, ~A, 100*A])
        [    3     9|   -3    -9]
        [    6    10|   -6   -10]
        [-----------+-----------]
        [-5/12   3/8|  300   900]
        [  1/4  -1/8|  600  1000]

    One can use constant entries:
        sage: block_matrix([1, A, 0, 1])
        [ 1  0| 3  9]
        [ 0  1| 6 10]
        [-----+-----]
        [ 0  0| 1  0]
        [ 0  0| 0  1]

    One can specify the number of rows or columns (optional for square number of matrices):
        sage: block_matrix([A, -A, ~A, 100*A], ncols=4)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

        sage: block_matrix([A, -A, ~A, 100*A], nrows=1)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

    It handle baserings nicely too:
        sage: R.<x> = ZZ['x']
        sage: block_matrix([1/2, A, 0, x-1])
        [  1/2     0|    3     9]
        [    0   1/2|    6    10]
        [-----------+-----------]
        [    0     0|x - 1     0]
        [    0     0|    0 x - 1]
        sage: block_matrix([1/2, A, 0, x-1]).parent()
        Full MatrixSpace of 4 by 4 dense matrices over Univariate Polynomial Ring in x over Rational Field

    Subdivisions are optional:
        sage: B = matrix(QQ, 2, 3, range(6))
        sage: block_matrix([~A, B, B, ~A], subdivide=False)
        [-5/12   3/8     0     1     2]
        [  1/4  -1/8     3     4     5]
        [    0     1     2 -5/12   3/8]
        [    3     4     5   1/4  -1/8]
    """
    # determine the block dimensions
    n = ZZ(len(sub_matrices))
    if nrows is None:
        if ncols is None:
            if n.is_square():
                nrows = ncols = n.sqrt()
            else:
                raise ValueError, "Must specify rows or cols for non-square block matrix."
        else:
            nrows = int(n/ncols)
    elif ncols is None:
        ncols = int(n/nrows)
    if nrows * ncols != n:
        raise ValueError, "Given number of rows (%s), columns (%s) incompatable with number of submatrices (%s)" % (nrows, ncols, n)

    # determine the sub-block dimensions
    row_heights = [None] * nrows
    col_widths = [None] * ncols
    for i in range(nrows):
        for j in range(0, ncols):
            M = sub_matrices[i*ncols+j]
            if is_Matrix(M):
                if row_heights[i] is None:
                    row_heights[i] = M.nrows()
                if col_widths[j] is None:
                    col_widths[j] = M.ncols()

    if None in row_heights or None in col_widths:
        for i in range(nrows):
            for j in range(0, ncols):
                x = sub_matrices[i*ncols+j]
                if not is_Matrix(x) and x: # must be square matrix
                    if row_heights[i] is None:
                        row_heights[i] = col_widths[j]
                    if col_widths[j] is None:
                        col_widths[j] = row_heights[i]

        if None in row_heights or None in col_widths:
            raise ValueError, "Insufficient information to determine dimensions."

    # determine the base ring
    base = ZZ
    for M in sub_matrices:
        R = M.base_ring() if is_Matrix(M) else M.parent()
        if R is not ZZ:
            base = sage.categories.pushout.pushout(base, R)

    # finally concatinate
    for i in range(nrows):
        for j in range(ncols):
            # coerce
            M = sub_matrices[i*ncols+j]
            if is_Matrix(M):
                if M.base_ring() is not base:
                    M = M.change_ring(base)
            else:
                M = matrix(base, row_heights[i], col_widths[j], M)
            # append
            if j == 0:
                row = M
            else:
                row = row.augment(M)
        if i == 0:
            big = row
        else:
            big = big.stack(row)

    if subdivide:
        big.subdivide(running_total(row_heights[:-1]),
                      running_total(col_widths[:-1]))

    return big

def block_diagonal_matrix(*sub_matrices, **kwds):
    """
    Create a block matrix whose diagonal block entries are given by sub_matrices,
    with zero elsewhere.

    See also \code{block_matrix}.

    EXAMPLES:
        sage: A = matrix(ZZ, 2, [1,2,3,4])
        sage: block_diagonal_matrix(A, A)
        [1 2|0 0]
        [3 4|0 0]
        [---+---]
        [0 0|1 2]
        [0 0|3 4]

    The sub-matrices need not be square:
        sage: B = matrix(QQ, 2, 3, range(6))
        sage: block_diagonal_matrix(~A, B)
        [  -2    1|   0    0    0]
        [ 3/2 -1/2|   0    0    0]
        [---------+--------------]
        [   0    0|   0    1    2]
        [   0    0|   3    4    5]
    """
    if isinstance(sub_matrices, (list, tuple)) and len(sub_matrices) == 1:
        sub_matrices = sub_matrices[0]
    n = len(sub_matrices)
    entries = [ZZ.zero_element()] * n**2
    for i in range(n):
        entries[n*i+i] = sub_matrices[i]
    return block_matrix(entries, **kwds)
