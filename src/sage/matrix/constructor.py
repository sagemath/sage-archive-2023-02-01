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
import types
import sage.rings.all as rings
import sage.matrix.matrix_space as matrix_space
from sage.structure.element import is_Vector
from sage.structure.sequence import Sequence
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.all import ZZ, QQ
from sage.misc.misc_c import running_total
from matrix import is_Matrix
from copy import copy

import sage.categories.pushout

def matrix(*args, **kwds):
    """
    Create a matrix.

    INPUT: The matrix command takes the entries of a matrix, optionally
    preceded by a ring and the dimensions of the matrix, and returns a
    matrix.

    The entries of a matrix can be specified as a flat list of
    elements, a list of lists (i.e., a list of rows), a list of Sage
    vectors, a callable object, or a dictionary having positions as
    keys and matrix entries as values (see the examples). If you pass
    in a callable object, then you must specify the number of rows and
    columns. You can create a matrix of zeros by passing an empty list
    or the integer zero for the entries.  To construct a multiple of
    the identity (`cI`), you can specify square dimensions and pass in
    `c`. Calling matrix() with a Sage object may return something that
    makes sense. Calling matrix() with a NumPy array will convert the
    array to a matrix.

    The ring, number of rows, and number of columns of the matrix can
    be specified by setting the ring, nrows, or ncols parameters or by
    passing them as the first arguments to the function in the order
    ring, nrows, ncols. The ring defaults to ZZ if it is not specified
    or cannot be determined from the entries. If the numbers of rows
    and columns are not specified and cannot be determined, then an
    empty 0x0 matrix is returned.


    -  ``ring`` - the base ring for the entries of the
       matrix.

    -  ``nrows`` - the number of rows in the matrix.

    -  ``ncols`` - the number of columns in the matrix.

    -  ``sparse`` - create a sparse matrix. This defaults
       to True when the entries are given as a dictionary, otherwise
       defaults to False.


    OUTPUT:

    a matrix

    EXAMPLES::

        sage: m=matrix(2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

    ::

        sage: m=matrix(2,3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m=matrix(QQ,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix(QQ, 3, 3, lambda i, j: i+j); m
        [0 1 2]
        [1 2 3]
        [2 3 4]
        sage: m = matrix(3, lambda i,j: i-j); m
        [ 0 -1 -2]
        [ 1  0 -1]
        [ 2  1  0]

    ::

        sage: matrix(QQ,2,3,lambda x, y: x+y)
        [0 1 2]
        [1 2 3]
        sage: matrix(QQ,3,2,lambda x, y: x+y)
        [0 1]
        [1 2]
        [2 3]

    ::

        sage: v1=vector((1,2,3))
        sage: v2=vector((4,5,6))
        sage: m=matrix([v1,v2]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m=matrix(QQ,2,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m=matrix(QQ,2,3,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m=matrix({(0,1): 2, (1,1):2/5}); m; m.parent()
        [  0   2]
        [  0 2/5]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field

    ::

        sage: m=matrix(QQ,2,3,{(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field

    ::

        sage: import numpy
        sage: n=numpy.array([[1,2],[3,4]],float)
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0]
        [3.0 4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Real Double Field

    ::

        sage: v = vector(ZZ, [1, 10, 100])
        sage: m=matrix(v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices over Integer Ring
        sage: m=matrix(GF(7), v); m; m.parent()
        [1 3 2]
        Full MatrixSpace of 1 by 3 dense matrices over Finite Field of size 7

    ::

        sage: g = graphs.PetersenGraph()
        sage: m = matrix(g); m; m.parent()
        [0 1 0 0 1 1 0 0 0 0]
        [1 0 1 0 0 0 1 0 0 0]
        [0 1 0 1 0 0 0 1 0 0]
        [0 0 1 0 1 0 0 0 1 0]
        [1 0 0 1 0 0 0 0 0 1]
        [1 0 0 0 0 0 0 1 1 0]
        [0 1 0 0 0 0 0 0 1 1]
        [0 0 1 0 0 1 0 0 0 1]
        [0 0 0 1 0 1 1 0 0 0]
        [0 0 0 0 1 0 1 1 0 0]
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring

    ::

        sage: matrix(ZZ, 10, 10, range(100), sparse=True).parent()
        Full MatrixSpace of 10 by 10 sparse matrices over Integer Ring

    ::

        sage: R = PolynomialRing(QQ, 9, 'x')
        sage: A = matrix(R, 3, 3, R.gens()); A
        [x0 x1 x2]
        [x3 x4 x5]
        [x6 x7 x8]
        sage: det(A)
        -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8

    TESTS::

        sage: m=matrix(); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(QQ); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m=matrix(QQ,2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: m=matrix(QQ,2,3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix([]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(QQ,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m=matrix(2,2,1); m; m.parent()
        [1 0]
        [0 1]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: m=matrix(QQ,2,2,1); m; m.parent()
        [1 0]
        [0 1]
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: m=matrix(2,3,0); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m=matrix(QQ,2,3,0); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix([[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m=matrix(QQ,2,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix(QQ,3,[[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: Number of rows does not match up with specified number.
        sage: m=matrix(QQ,2,3,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix(QQ,2,4,[[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: Number of columns does not match up with specified number.
        sage: m=matrix([(1,2,3),(4,5,6)]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m=matrix([1,2,3,4,5,6]); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Integer Ring
        sage: m=matrix((1,2,3,4,5,6)); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Integer Ring
        sage: m=matrix(QQ,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Rational Field
        sage: m=matrix(QQ,3,2,[1,2,3,4,5,6]); m; m.parent()
        [1 2]
        [3 4]
        [5 6]
        Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        sage: m=matrix(QQ,2,4,[1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: entries has the wrong length
        sage: m=matrix(QQ,5,[1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        TypeError: entries has the wrong length
        sage: m=matrix({(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m=matrix(QQ,{(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m=matrix(QQ,3,{(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        [0 0 0]
        Full MatrixSpace of 3 by 3 sparse matrices over Rational Field
        sage: m=matrix(QQ,3,4,{(1,1): 2}); m; m.parent()
        [0 0 0 0]
        [0 2 0 0]
        [0 0 0 0]
        Full MatrixSpace of 3 by 4 sparse matrices over Rational Field
        sage: m=matrix(QQ,2,{(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m=matrix(QQ,1,{(1,1): 2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: m=matrix({}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m=matrix(QQ,{}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Rational Field
        sage: m=matrix(QQ,2,{}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m=matrix(QQ,2,3,{}); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field
        sage: m=matrix(2,{}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m=matrix(2,3,{}); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Integer Ring
        sage: m=matrix(0); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,2); m; m.parent()
        []
        Full MatrixSpace of 0 by 2 dense matrices over Integer Ring
        sage: m=matrix(2,0); m; m.parent()
        []
        Full MatrixSpace of 2 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,[1]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: entries has the wrong length
        sage: m=matrix(1,0,[]); m; m.parent()
        []
        Full MatrixSpace of 1 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,1,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 1 dense matrices over Integer Ring
        sage: m=matrix(0,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,{}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m=matrix(0,{(1,1):2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: m=matrix(2,0,{(1,1):2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: import numpy
        sage: n=numpy.array([[numpy.complex(0,1),numpy.complex(0,2)],[3,4]],complex)
        sage: m=matrix(n); m; m.parent()
        [1.0*I 2.0*I]
        [  3.0   4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field
        sage: n=numpy.array([[1,2],[3,4]],'int32')
        sage: m=matrix(n); m; m.parent()
        [1 2]
        [3 4]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: n = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float32')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n=numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float64')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n=numpy.array([[1,2,3],[4,5,6],[7,8,9]],'complex64')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: n=numpy.array([[1,2,3],[4,5,6],[7,8,9]],'complex128')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: a = matrix([[1,2],[3,4]])
        sage: b = matrix(a.numpy()); b
        [1 2]
        [3 4]
        sage: a == b
        True
        sage: c = matrix(a.numpy('float32')); c
        [1.0 2.0]
        [3.0 4.0]
        sage: matrix(numpy.array([[5]]))
        [5]
        sage: v = vector(ZZ, [1, 10, 100])
        sage: m=matrix(ZZ['x'], v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices over Univariate Polynomial Ring in x over Integer Ring
        sage: matrix(ZZ, 10, 10, range(100)).parent()
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring
        sage: m = matrix(GF(7), [[1/3,2/3,1/2], [3/4,4/5,7]]); m; m.parent()
        [5 3 4]
        [6 5 0]
        Full MatrixSpace of 2 by 3 dense matrices over Finite Field of size 7
        sage: m = matrix([[1,2,3], [RDF(2), CDF(1,2), 3]]); m; m.parent()
        [        1.0         2.0         3.0]
        [        2.0 1.0 + 2.0*I         3.0]
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field
        sage: m=matrix(3,3,1/2); m; m.parent()
        [1/2   0   0]
        [  0 1/2   0]
        [  0   0 1/2]
        Full MatrixSpace of 3 by 3 dense matrices over Rational Field
        sage: matrix([[1],[2,3]])
        Traceback (most recent call last):
        ...
        ValueError: List of rows is not valid (rows are wrong types or lengths)
        sage: matrix([[1],2])
        Traceback (most recent call last):
        ...
        ValueError: List of rows is not valid (rows are wrong types or lengths)
        sage: matrix(vector(RR,[1,2,3])).parent()
        Full MatrixSpace of 1 by 3 dense matrices over Real Field with 53 bits of precision
        sage: matrix(ZZ, [[0] for i in range(10^6)]).is_zero() # see #10158
        True

    AUTHORS:

    - ??: Initial implementation

    - Jason Grout (2008-03): almost a complete rewrite, with bits and
      pieces from the original implementation
    """
    args = list(args)
    sparse = kwds.get('sparse',False)
    # if the first object already knows how to make itself into a
    # matrix, try that, defaulting to a matrix over the integers.
    if len(args) == 1 and hasattr(args[0], '_matrix_'):
        try:
                return args[0]._matrix_(sparse=sparse)
        except TypeError:
                return args[0]._matrix_()
    elif len(args) == 2:
        if hasattr(args[0], '_matrix_'):
            try:
                return args[0]._matrix_(args[1], sparse=sparse)
            except TypeError:
                return args[0]._matrix_(args[1])
        elif hasattr(args[1], '_matrix_'):
            try:
                return args[1]._matrix_(args[0], sparse=sparse)
            except TypeError:
                return args[1]._matrix_(args[0])

    if len(args) == 0:
        # if nothing was passed return the empty matrix over the
        # integer ring.
        return matrix_space.MatrixSpace(rings.ZZ, 0, 0, sparse=sparse)([])

    if len(args) >= 1 and rings.is_Ring(args[0]):
        # A ring is specified
        if kwds.get('ring', args[0]) != args[0]:
            raise ValueError, "Specified rings are not the same"
        else:
            ring = args[0]
            args.pop(0)
    else:
        ring = kwds.get('ring', None)

    if len(args) >= 1:
        # check to see if the number of rows is specified
        try:
            import numpy
            if isinstance(args[0], numpy.ndarray):
                raise TypeError
            nrows = int(args[0])
            args.pop(0)
            if kwds.get('nrows', nrows) != nrows:
                raise ValueError, "Number of rows specified in two places and they are not the same"
        except TypeError:
            nrows = kwds.get('nrows', None)
    else:
        nrows = kwds.get('nrows', None)

    if len(args) >= 1:
        # check to see if additionally, the number of columns is specified
        try:
            import numpy
            if isinstance(args[0], numpy.ndarray):
                raise TypeError
            ncols = int(args[0])
            args.pop(0)
            if kwds.get('ncols', ncols) != ncols:
                raise ValueError, "Number of columns specified in two places and they are not the same"
        except TypeError:
            ncols = kwds.get('ncols', None)
    else:
        ncols = kwds.get('ncols', None)


    # Now we've taken care of initial ring, nrows, and ncols arguments.
    # We've also taken care of the Sage object case.

    # Now the rest of the arguments are a list of rows, a flat list of
    # entries, a callable, a dict, a numpy array, or a single value.
    if len(args) == 0:
        # If no entries are specified, pass back a zero matrix
        entries = 0
        entry_ring = rings.ZZ
    elif len(args) == 1:
        if isinstance(args[0], (types.FunctionType, types.LambdaType, types.MethodType)):
            if ncols is None and nrows is None:
                raise ValueError, "When passing in a callable, the dimensions of the matrix must be specified"
            if ncols is None:
                ncols = nrows
            elif nrows is None:
                nrows = ncols

            f = args[0]
            args[0] = [[f(i,j) for j in range(ncols)] for i in range(nrows)]

        if isinstance(args[0], (list, tuple)):
            if len(args[0]) == 0:
                # no entries are specified, pass back the zero matrix
                entries = 0
                entry_ring = rings.ZZ
            elif isinstance(args[0][0], (list, tuple)) or is_Vector(args[0][0]):
                # Ensure we have a list of lists, each inner list having the same number of elements
                first_len = len(args[0][0])
                if not all( (isinstance(v, (list, tuple)) or is_Vector(v)) and len(v) == first_len for v in args[0]):
                    raise ValueError, "List of rows is not valid (rows are wrong types or lengths)"
                # We have a list of rows or vectors
                if nrows is None:
                    nrows = len(args[0])
                elif nrows != len(args[0]):
                    raise ValueError, "Number of rows does not match up with specified number."
                if ncols is None:
                    ncols = len(args[0][0])
                elif ncols != len(args[0][0]):
                    raise ValueError, "Number of columns does not match up with specified number."

                entries = []
                for v in args[0]:
                    entries.extend(v)

            else:
                # We have a flat list; figure out nrows and ncols
                if nrows is None:
                    nrows = 1

                if nrows > 0:
                    if ncols is None:
                        ncols = len(args[0]) // nrows
                    elif ncols != len(args[0]) // nrows:
                        raise ValueError, "entries has the wrong length"
                elif len(args[0]) > 0:
                    raise ValueError, "entries has the wrong length"

                entries = args[0]

            if nrows > 0 and ncols > 0 and ring is None:
                entries, ring = prepare(entries)

        elif isinstance(args[0], dict):
            # We have a dictionary
            # default to sparse
            sparse = kwds.get('sparse', True)
            if len(args[0]) == 0:
                # no entries are specified, pass back the zero matrix
                entries = 0
            else:
                entries, entry_ring = prepare_dict(args[0])
                if nrows is None:
                    nrows = nrows_from_dict(entries)
                    ncols = ncols_from_dict(entries)
                # note that ncols can still be None if nrows is set --
                # it will be assigned nrows down below.

            # See the construction after the numpy case below.
        else:
            import numpy
            if isinstance(args[0], numpy.ndarray):
                num_array = args[0]
                str_dtype = str(num_array.dtype)

                if not( num_array.flags.c_contiguous is True or num_array.flags.f_contiguous is True):
                    raise TypeError('numpy matrix must be either c_contiguous or f_contiguous')
                if str_dtype.count('float32')==1:
                    m=matrix(RDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy32(num_array)

                elif str_dtype.count('float64')==1:
                    m=matrix(RDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy(num_array)

                elif str_dtype.count('complex64')==1:
                    m=matrix(CDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy32(num_array)

                elif str_dtype.count('complex128')==1:
                    m=matrix(CDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy(num_array)

                elif str_dtype.count('int') == 1:
                    m = matrix(ZZ, [list(row) for row in list(num_array)])

                elif str_dtype.count('object') == 1:
                    #Get the raw nested list from the numpy array
                    #and feed it back into matrix
                    try:
                        return matrix( [list(row) for row in list(num_array)])
                    except TypeError:
                        raise TypeError("cannot convert NumPy matrix to Sage matrix")
                else:
                    raise TypeError("cannot convert NumPy matrix to Sage matrix")

                return m
            elif nrows is not None and ncols is not None:
                # assume that we should just pass the thing into the
                # MatrixSpace constructor and hope for the best
                # This is not documented, but it is doctested
                if ring is None:
                    entry_ring = args[0].parent()
                entries = args[0]
            else:
                raise ValueError, "Invalid matrix constructor.  Type matrix? for help"
    else:
        raise ValueError, "Invalid matrix constructor.  Type matrix? for help"

    if nrows is None:
        nrows = 0
    if ncols is None:
        ncols = nrows


    if ring is None:
        try:
            ring = entry_ring
        except NameError:
            ring = rings.ZZ

    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse=sparse)(entries)



def prepare(w):
    """
    Given a list w of numbers, find a common ring that they all
    canonically map to, and return the list of images of the elements
    of w in that ring along with the ring.

    This is for internal use by the matrix function.

    INPUT:

    - ``w`` - list

    OUTPUT:

    list, ring

    EXAMPLES::

        sage: sage.matrix.constructor.prepare([-2, Mod(1,7)])
        ([5, 1], Ring of integers modulo 7)

    Notice that the elements must all canonically coerce to a common
    ring (since Sequence is called)::

        sage: sage.matrix.constructor.prepare([2/1, Mod(1,7)])
        Traceback (most recent call last):
        ...
        TypeError: unable to find a common ring for all elements
    """
    if 0 == len(w):
        return Sequence([], rings.ZZ), rings.ZZ
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
    """
    Given a dictionary w of numbers, find a common ring that they all
    canonically map to, and return the dictionary of images of the
    elements of w in that ring along with the ring.

    This is for internal use by the matrix function.

    INPUT:

    - ``w`` - dict

    OUTPUT:

    dict, ring

    EXAMPLES::

        sage: sage.matrix.constructor.prepare_dict({(0,1):2, (4,10):Mod(1,7)})
        ({(0, 1): 2, (4, 10): 1}, Ring of integers modulo 7)
    """
    Z = w.items()
    X = [x for _, x in Z]
    entries, ring = prepare(X)
    return dict([(Z[i][0],entries[i]) for i in xrange(len(entries))]), ring

def nrows_from_dict(d):
    """
    Given a dictionary that defines a sparse matrix, return the number
    of rows that matrix should have.

    This is for internal use by the matrix function.

    INPUT:

    - ``d`` - dict

    OUTPUT:

    integer

    EXAMPLES::

        sage: sage.matrix.constructor.nrows_from_dict({})
        0

    Here the answer is 301 not 300, since there is a 0-th row.

    ::
        sage: sage.matrix.constructor.nrows_from_dict({(300,4):10})
        301
    """
    if 0 == len(d):
        return 0
    return max([0] + [ij[0] for ij in d.keys()]) + 1

def ncols_from_dict(d):
    """
    Given a dictionary that defines a sparse matrix, return the number
    of columns that matrix should have.

    This is for internal use by the matrix function.

    INPUT:

    - ``d`` - dict

    OUTPUT:

    integer

    EXAMPLES::

        sage: sage.matrix.constructor.ncols_from_dict({})
        0

    Here the answer is 301 not 300, since there is a 0-th row.

    ::

        sage: sage.matrix.constructor.ncols_from_dict({(4,300):10})
        301
    """
    if 0 == len(d):
        return 0
    return max([0] + [ij[1] for ij in d.keys()]) + 1

Matrix = matrix


def random_matrix(ring, nrows, ncols=None, algorithm='randomize', *args, **kwds):
    r"""
    Return a random matrix with entries in a specified ring, and possibly with additional properties.

    INPUT:

    -  ``ring`` - base ring for entries of the matrix

    -  ``nrows`` - Integer; number of rows

    -  ``ncols`` - (default: ``None``); number of columns; if ``None``
       defaults to ``nrows``

    -  ``algorithm`` - (default: ``randomize``); determines what properties
       the matrix will have.  See examples below for possible additional
       arguments.

       -  ``randomize`` - randomize the elements of the matrix, possibly
          controlling the density of non-zero entries.

       -  ``echelon_form`` - creates a matrix in echelon form

       -  ``echelonizable`` - creates a matrix that has a predictable
          echelon form

       - ``subspaces`` - creates a matrix whose four subspaces, when
         explored, have reasonably sized, integral valued, entries.

       - ``unimodular`` - creates a matrix of determinant 1.

       - ``diagonalizable`` - creates a diagonalizable matrix whose
         eigenvectors, if computed by hand, will have only integer
         entries.

    -  ``*args, **kwds`` - arguments and keywords to describe additional properties.
       See more detailed documentation below.

    .. note::

        When constructing matrices with random entries and no additional properties
        (i.e. when ``algorithm='randomize'``), most of the randomness
        is controlled by the ``random_element`` method for elements of the
        base ring of the matrix, so the documentation of that method may be
        relevant or useful.  Also, the default is to not create zero entries,
        unless the ``density`` keyword is set to something strictly less
        than one.

    EXAMPLES:

    Random integer matrices.  With no arguments, the majority of the entries
    are -1 and 1, never zero, and rarely "large." ::

        sage: random_matrix(ZZ, 5, 5)
        [ -8   2   0   0   1]
        [ -1   2   1 -95  -1]
        [ -2 -12   0   0   1]
        [ -1   1  -1  -2  -1]
        [  4  -4  -6   5   0]

    The ``distribution`` keyword  set to ``uniform`` will limit values
    between -2 and 2, and never zero. ::

        sage: random_matrix(ZZ, 5, 5, distribution='uniform')
        [ 1  0 -2  1  1]
        [ 1  0  0  0  2]
        [-1 -2  0  2 -2]
        [-1 -1  1  1  2]
        [ 0 -2 -1  0  0]

    The ``x`` and ``y`` keywords can be used to distribute entries uniformly.
    When both are used ``x`` is the minimum and ``y`` is one greater than the the maximum.
    But still entries are never zero, even if the range contains zero. ::

        sage: random_matrix(ZZ, 4, 8, x=70, y=100)
        [81 82 70 81 78 71 79 94]
        [80 98 89 87 91 94 94 77]
        [86 89 85 92 95 94 72 89]
        [78 80 89 82 94 72 90 92]

        sage: random_matrix(ZZ, 3, 7, x=-5, y=5)
        [-3  3  1 -5  3  1  2]
        [ 3  3  0  3 -5 -2  1]
        [ 0 -2 -2  2 -3 -4 -2]

    If only ``x`` is given, then it is used as the upper bound of a range starting at 0. ::

        sage: random_matrix(ZZ, 5, 5, x=25)
        [20 16  8  3  8]
        [ 8  2  2 14  5]
        [18 18 10 20 11]
        [19 16 17 15  7]
        [ 0 24  3 17 24]

    To allow, and control, zero entries use the ``density`` keyword at a value
    strictly below the default of 1.0, even if distributing entries across an
    interval that does not contain zero already.  Note that for a square matrix it
    is only necessary to set a single dimension. ::

        sage: random_matrix(ZZ, 5, x=-10, y=10, density=0.75)
        [-6  1  0  0  0]
        [ 9  0  0  4  1]
        [-6  0  0 -8  0]
        [ 0  4  0  6  0]
        [ 1 -9  0  0 -8]

        sage: random_matrix(ZZ, 5, x=20, y=30, density=0.75)
        [ 0 28  0 27  0]
        [25 28 20  0  0]
        [ 0 21  0 21  0]
        [ 0 28 22  0  0]
        [ 0  0  0 26 24]

    It is possible to construct sparse matrices, where it may now be advisable
    (but not required) to control the density of nonzero entries. ::

        sage: A=random_matrix(ZZ, 5, 5)
        sage: A.is_sparse()
        False
        sage: A=random_matrix(ZZ, 5, 5, sparse=True)
        sage: A.is_sparse()
        True

        sage: random_matrix(ZZ, 5, 5, density=0.3, sparse=True)
        [ 4  0  0  0 -1]
        [ 0  0  0  0 -7]
        [ 0  0  2  0  0]
        [ 0  0  1  0 -4]
        [ 0  0  0  0  0]

    For algorithm testing you might want to control the number of bits,
    say 10,000 entries, each limited to 16 bits.  ::

        sage: A = random_matrix(ZZ, 100, 100, x=2^16); A
        100 x 100 dense matrix over Integer Ring (type 'print A.str()' to see all of the entries)

    Random rational matrices.  Now ``num_bound`` and ``den_bound`` control the
    generation of random elements, by specifying limits on the absolute value of
    numerators and denominators (respectively).  Entries will be positive and
    negative (map the absolute value function through the entries to get all
    positive values), and zeros are avoided unless the density is set.  If either
    the numerator or denominator bound (or both) is not used, then the values
    default to the distribution for `ZZ` described above that is most frequently
    positive or negative one. ::

        sage: random_matrix(QQ, 2, 8, num_bound=20, den_bound=4)
        [ -1/2     6    13   -12  -2/3  -1/4     5     5]
        [ -9/2   5/3    19  15/2  19/2  20/3 -13/4     0]

        sage: random_matrix(QQ, 4, density = 0.5, sparse=True)
        [    0    71     0  -1/2]
        [    0     0     0     0]
        [31/85     0 -31/2     0]
        [    1  -1/4     0     0]

        sage: A = random_matrix(QQ, 3, 10, num_bound = 99, den_bound = 99)
        sage: positives = map(abs, A.list())
        sage: matrix(QQ, 3, 10, positives)
        [61/18 47/41  1/22   1/2 75/68   6/7     1   1/2 72/41   7/3]
        [33/13   9/2 40/21 45/46 17/22     1 70/79 97/71  7/24  12/5]
        [ 13/8  8/25   1/3 61/14 92/45  4/85  3/38 95/16 82/71   1/5]

        sage: random_matrix(QQ, 4, 10, den_bound = 10)
        [  -1    0  1/8  1/6  2/9 -1/6  1/5 -1/8  1/5 -1/5]
        [ 1/9  1/5   -1  2/9  1/4 -1/7  1/8 -1/9    0    2]
        [ 2/3    2  1/8   -2    0    0   -2    2    0 -1/2]
        [   0    2    1 -2/3    0    0  1/6    0 -1/3 -2/9]

    Random matrices over other rings.  Several classes of matrices have specialized
    ``randomize()`` methods.  You can locate these with the Sage command::

        search_def('randomize')

    The default implementation of :meth:`~sage.matrix.matrix2.randomize` relies
    on the ``random_element()`` method for the base ring.  The ``density`` and
    ``sparse`` keywords behave as described above. ::

        sage: K.<a>=FiniteField(3^2)
        sage: random_matrix(K, 2, 5)
        [      1     2*a       1   a + 2       2]
        [      a 2*a + 1       0       2       1]

        sage: random_matrix(RR, 3, 4, density=0.66)
        [ 0.000000000000000  0.566500636438206 0.0870635178173962  0.000000000000000]
        [-0.662290145671671  0.000000000000000  0.475667133865666  0.000000000000000]
        [-0.276405104068647  0.000000000000000  0.000000000000000 -0.636689607643642]

        sage: A = random_matrix(ComplexField(32), 3, density=0.8, sparse=True); A
        [                             0   -0.193242896 + 0.460619713*I                              0]
        [   0.909948633 + 0.611092515*I 0.00128001347 + 0.0659103142*I   0.199796580 + 0.0955426861*I]
        [                             0                              0    0.354729790 - 0.184624095*I]
        sage: A.is_sparse()
        True

    Random matrices in echelon form.  The ``algorithm='echelon_form'`` keyword,
    along with a requested number of non-zero rows (``num_pivots``) will return
    a random matrix in echelon form.  When the base ring is ``QQ`` the result has integer
    entries.  Other exact rings may be also specified. ::

        sage: A=random_matrix(QQ, 4, 8, algorithm='echelon_form', num_pivots=3); A # random
        [ 1 -5  0 -2  0  1  1 -2]
        [ 0  0  1 -5  0 -3 -1  0]
        [ 0  0  0  0  1  2 -2  1]
        [ 0  0  0  0  0  0  0  0]
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (4, 8)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True
        sage: A.rank()
        3
        sage: A==A.rref()
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_rref_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_rref_matrix

    Random matrices with predictable echelon forms.  The ``algorithm='echelonizable'``
    keyword, along with a requested rank (``rank``) and optional size control
    (``upper_bound``) will return a random matrix in echelon form.  When the
    base ring is ``ZZ`` or ``QQ`` the result has integer entries, whose magnitudes
    can be limited by the value of ``upper_bound``, and the echelon form of the
    matrix also has integer entries.  Other exact rings may be also
    specified, but there is no notion of controlling the size.  Square matrices
    of full rank generated by this function always have determinant one, and
    can be constructed with the ``unimodular`` keyword. ::

        sage: A=random_matrix(QQ, 4, 8, algorithm='echelonizable', rank=3, upper_bound=60); A # random
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (4, 8)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True
        sage: A.rank()
        3
        sage: all([abs(x)<60 for x in A.list()])
        True
        sage: A.rref() in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_echelonizable_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_echelonizable_matrix

    Random diagonalizable matrices.  The ``algorithm='diagonalizable'`` keyword,
    along with a requested matrix size (``size``) and optional lists of
    eigenvalues (``eigenvalues``) and the corresponding eigenspace
    dimensions (``dimensions``) will return a random diagonalizable matrix.
    When the eigenvalues and dimensions are not specified the result will have
    randomly generated values for both that fit with the designated size. ::

        sage: A=random_matrix(QQ, 5, algorithm='diagonalizable', eigenvalues=[2,3,-1], dimensions=[1,2,2]); A # random
        sage: all([x in ZZ for x in (A-(2*identity_matrix(5))).rref().list()])
        True
        sage: all([x in ZZ for x in (A-(3*identity_matrix(5))).rref().list()])
        True
        sage: all([x in ZZ for x in (A-(-1*identity_matrix(5))).rref().list()])
        True
        sage: A.jordan_form()
        [ 2| 0| 0| 0| 0]
        [--+--+--+--+--]
        [ 0| 3| 0| 0| 0]
        [--+--+--+--+--]
        [ 0| 0| 3| 0| 0]
        [--+--+--+--+--]
        [ 0| 0| 0|-1| 0]
        [--+--+--+--+--]
        [ 0| 0| 0| 0|-1]

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_diagonalizable_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_diagonalizable_matrix

    Random matrices with predictable subspaces.  The ``algorithm='subspaces'``
    keyword, along with an optional rank (``rank``) will return
    a matrix whose natural basis vectors for its four fundamental subspaces, if computed as
    described in the documentation of the :func:`~sage.matrix.constructor.random_subspaces_matrix`
    contain only integer entries.  If ``rank``, is not set, the
    rank of the matrix will be generated randomly. ::

        sage: B=random_matrix(QQ, 5, 6, algorithm='subspaces', rank=3); B #random
        sage: B_expanded=B.augment(identity_matrix(5)).rref()
        sage: (B.nrows(), B.ncols())
        (5, 6)
        sage: all([x in ZZ for x in B_expanded.list()])
        True
        sage: C=B_expanded.submatrix(0,0,B.nrows()-B.nullity(),B.ncols())
        sage: L=B_expanded.submatrix(B.nrows()-B.nullity(),B.ncols())
        sage: B.right_kernel()==C.right_kernel()
        True
        sage: B.row_space()==C.row_space()
        True
        sage: B.column_space()==L.right_kernel()
        True
        sage: B.left_kernel()==L.row_space()
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_subspaces_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_subspaces_matrix

    Random unimodular matrices.  The ``algorithm='unimodular'``
    keyword, along with an optional entry size control (``upper_bound``)
    will return a matrix of determinant 1. When the base ring is ``ZZ``
    or ``QQ`` the result has integer entries, whose magnitudes
    can be limited by the value of ``upper_bound``. ::

        sage: C=random_matrix(QQ, 5, algorithm='unimodular', upper_bound=70); C # random
        sage: det(C)
        1
        sage: C.base_ring()
        Rational Field
        sage: (C.nrows(), C.ncols())
        (5, 5)
        sage: all([abs(x)<70 for x in C.list()])
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_unimodular_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_unimodular_matrix

    TESTS:

    We return an error for a bogus value of ``algorithm``::

        sage: random_matrix(ZZ, 5, algorithm = 'bogus')
        Traceback (most recent call last):
        ...
        ValueError: random matrix algorithm "bogus" is not recognized

    AUTHOR:

    - William Stein (2007-02-06)

    - Rob Beezer (2010-08-25) Documentation, code to allow additional types of output
    """
    if ncols is None:
        ncols = nrows
    sparse = kwds.pop('sparse', False)
    # Construct the parent of the desired matrix
    parent = matrix_space.MatrixSpace(ring, nrows, ncols, sparse=sparse)
    if algorithm == 'randomize':
        density = kwds.pop('density', None)
        # zero matrix is immutable, copy is mutable
        A = copy(parent.zero_matrix())
        if density is None:
            A.randomize(density=float(1), nonzero=False, *args, **kwds)
        else:
            A.randomize(density=density, nonzero=True, *args, **kwds)
        return A
    elif algorithm == 'echelon_form':
        return random_rref_matrix(parent, *args, **kwds)
    elif algorithm == 'echelonizable':
        return random_echelonizable_matrix(parent, *args, **kwds)
    elif algorithm == 'diagonalizable':
        return random_diagonalizable_matrix(parent, *args, **kwds)
    elif algorithm == 'subspaces':
        return random_subspaces_matrix(parent, *args, **kwds)
    elif algorithm == 'unimodular':
        return random_unimodular_matrix(parent, *args, **kwds)
    else:
        raise ValueError('random matrix algorithm "%s" is not recognized' % algorithm)


def diagonal_matrix(arg0=None, arg1=None, arg2=None, sparse=None):
    """
    INPUT:

    Supported formats

    1. diagonal_matrix(diagonal_entries, [sparse=True]):
       diagonal matrix with flat list of entries

    2. diagonal_matrix(nrows, diagonal_entries, [sparse=True]):
       diagonal matrix with flat list of entries and the rest zeros

    3. diagonal_matrix(ring, diagonal_entries, [sparse=True]):
       diagonal matrix over specified ring with flat list of entries

    4. diagonal_matrix(ring, nrows, diagonal_entries, [sparse=True]):
       diagonal matrix over specified ring with flat
       list of entries and the rest zeros

    5. diagonal_matrix(vect, [sparse=True]):
       diagonal matrix with entries taken from a vector

    The sparse option is optional, must be explicitly named (i.e.,
    sparse=True), and may be either True or False.

    EXAMPLES:

    Input format 1::

        sage: diagonal_matrix([1,2,3])
        [1 0 0]
        [0 2 0]
        [0 0 3]

    Input format 2::

        sage: diagonal_matrix(3, [1,2])
        [1 0 0]
        [0 2 0]
        [0 0 0]

    Input format 3::

        sage: diagonal_matrix(GF(3), [1,2,3])
        [1 0 0]
        [0 2 0]
        [0 0 0]

    Input format 4::

        sage: diagonal_matrix(GF(3), 3, [8,2])
        [2 0 0]
        [0 2 0]
        [0 0 0]

    Input format 5::

        sage: diagonal_matrix(vector(GF(3),[1,2,3]))
        [1 0 0]
        [0 2 0]
        [0 0 0]

    """
    ring = None
    if isinstance(arg0, (list, tuple)):
        # Format 1
        v = arg0
        nrows = len(v)
    elif isinstance(arg0, (int, long, rings.Integer)):
        # Format 2
        nrows = arg0
        v = arg1
    elif rings.is_Ring(arg0):
        ring = arg0
        if isinstance(arg1, (list, tuple)):
            # Format 3
            v = arg1
            nrows = len(v)
        else:
            # Format 4
            nrows = arg1
            v = arg2
    elif is_Vector(arg0):
        # Format 5
        v = list(arg0)
        nrows = len(v)
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


def identity_matrix(ring, n=0, sparse=False):
    r"""
    Return the `n \times n` identity matrix over the given
    ring.

    The default ring is the integers.

    EXAMPLES::

        sage: M = identity_matrix(QQ, 2); M
        [1 0]
        [0 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: M = identity_matrix(2); M
        [1 0]
        [0 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: M.is_mutable()
        True
        sage: M = identity_matrix(3, sparse=True); M
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: M.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
        sage: M.is_mutable()
        True
    """
    if isinstance(ring, (int, long, rings.Integer)):
        n = ring
        ring = rings.ZZ
    return matrix_space.MatrixSpace(ring, n, n, sparse)(1)


def zero_matrix(ring, nrows, ncols=None, sparse=False):
    r"""
    Return the `nrows \times ncols` zero matrix over the given
    ring.

    The default ring is the integers.

    EXAMPLES::

        sage: M = zero_matrix(QQ, 2); M
        [0 0]
        [0 0]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: M = zero_matrix(2, 3); M
        [0 0 0]
        [0 0 0]
        sage: M.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: M.is_mutable()
        True
        sage: M = zero_matrix(3, 1, sparse=True); M
        [0]
        [0]
        [0]
        sage: M.parent()
        Full MatrixSpace of 3 by 1 sparse matrices over Integer Ring
        sage: M.is_mutable()
        True
    """
    if isinstance(ring, (int, long, rings.Integer)):
        nrows, ncols = (ring, nrows)
        ring = rings.ZZ
    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse)(0)

def ones_matrix(ring, nrows, ncols=None, sparse=False):
    r"""
    Return the `nrows \times ncols` all-ones matrix over the given
    ring.

    The default ring is the integers.

    EXAMPLES::

        sage: M = ones_matrix(QQ, 2); M
        [1 1]
        [1 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: M = ones_matrix(2, 3); M
        [1 1 1]
        [1 1 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: M = ones_matrix(3, 1, sparse=True); M
        [1]
        [1]
        [1]
        sage: M.parent()
        Full MatrixSpace of 3 by 1 sparse matrices over Integer Ring
    """
    if isinstance(ring, (int, long, rings.Integer)):
        nrows, ncols = (ring, nrows)
        ring = rings.ZZ
    if ncols is None:
        nents = nrows**2
    else:
        nents = nrows*ncols
    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse).matrix([1]*nents)


def block_matrix(sub_matrices, nrows=None, ncols=None, subdivide=True):
    """
    Returns a larger matrix made by concatenating the sub_matrices
    (rows first, then columns). For example, the matrix

    ::

        [ A B ]
        [ C D ]

    is made up of submatrices A, B, C, and D.

    INPUT:


    -  ``sub_matrices`` - matrices (must be of the correct
       size, or constants)

    -  ``nrows`` - (optional) the number of block rows

    -  ``ncols`` - (optional) the number of block cols

    -  ``subdivide`` - boolean, whether or not to add
       subdivision information to the matrix


    EXAMPLES::

        sage: A = matrix(QQ, 2, 2, [3,9,6,10])
        sage: block_matrix([A, -A, ~A, 100*A])
        [    3     9|   -3    -9]
        [    6    10|   -6   -10]
        [-----------+-----------]
        [-5/12   3/8|  300   900]
        [  1/4  -1/8|  600  1000]

    One can use constant entries::

        sage: block_matrix([1, A, 0, 1])
        [ 1  0| 3  9]
        [ 0  1| 6 10]
        [-----+-----]
        [ 0  0| 1  0]
        [ 0  0| 0  1]

    A zero entry may represent any square or non-square zero matrix::

        sage: B = matrix(QQ, 1, 1, [ 1 ] )
        sage: C = matrix(QQ, 2, 2, [ 2, 3, 4, 5 ] )
        sage: block_matrix([B, 0, 0, C])
        [1|0 0]
        [-+---]
        [0|2 3]
        [0|4 5]

    One can specify the number of rows or columns (optional for square
    number of matrices)::

        sage: block_matrix([A, -A, ~A, 100*A], ncols=4)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

    ::

        sage: block_matrix([A, -A, ~A, 100*A], nrows=1)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

    It handle base rings nicely too::

        sage: R.<x> = ZZ['x']
        sage: block_matrix([1/2, A, 0, x-1])
        [  1/2     0|    3     9]
        [    0   1/2|    6    10]
        [-----------+-----------]
        [    0     0|x - 1     0]
        [    0     0|    0 x - 1]
        sage: block_matrix([1/2, A, 0, x-1]).parent()
        Full MatrixSpace of 4 by 4 dense matrices over Univariate Polynomial Ring in x over Rational Field

    Subdivisions are optional::

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
        raise ValueError, "Given number of rows (%s), columns (%s) incompatible with number of submatrices (%s)" % (nrows, ncols, n)

    # empty matrix
    if n == 0:
      ans = Matrix(ZZ,0,0,[])
      if subdivide:
        ans.subdivide([0]*(nrows-1),[0]*(ncols-1))
      return ans

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

    # finally concatenate
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
    Create a block matrix whose diagonal block entries are given by
    sub_matrices, with zero elsewhere.

    See also ``block_matrix``.

    EXAMPLES::

        sage: A = matrix(ZZ, 2, [1,2,3,4])
        sage: block_diagonal_matrix(A, A)
        [1 2|0 0]
        [3 4|0 0]
        [---+---]
        [0 0|1 2]
        [0 0|3 4]

    The sub-matrices need not be square::

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

def jordan_block(eigenvalue, size, sparse=False):
    r"""
    Form the Jordan block with the specified size associated with the
    eigenvalue.

    INPUT:


    -  ``eigenvalue`` - eigenvalue for the diagonal entries
       of the block

    -  ``size`` - size of the Jordan block

    -  ``sparse`` - (default False) if True, return a
       sparse matrix


    EXAMPLE::

        sage: jordan_block(5, 3)
        [5 1 0]
        [0 5 1]
        [0 0 5]
    """
    block = diagonal_matrix([eigenvalue]*size, sparse=sparse)
    for i in xrange(size-1):
        block[i,i+1]=1
    return block

def random_rref_matrix(parent, num_pivots):
    r"""
    Generate a matrix in reduced row-echelon form with a specified number of non-zero rows.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``num_pivots`` - The number of non-zero rows in the result, i.e. the rank.

    OUTPUT:

    A matrix in reduced row echelon form with ``num_pivots`` non-zero rows. If the
    base ring is `ZZ` or `QQ` then the entries are all integers.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='echelon_form'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    Matrices generated are in reduced row-echelon form with specified rank. If the
    base ring is `QQ` the result has only integer entries.  ::

        sage: from sage.matrix.constructor import random_rref_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5, 6)
        sage: A=random_rref_matrix(matrix_space, num_pivots=4); A # random
        [ 1  0  0 -6  0 -3]
        [ 0  1  0  2  0  3]
        [ 0  0  1 -4  0 -2]
        [ 0  0  0  0  1  3]
        [ 0  0  0  0  0  0]
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (5, 6)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 5, 6)
        True
        sage: A.rank()
        4
        sage: A==A.rref()
        True

    Matrices can be generated over other exact rings. ::

        sage: B=random_matrix(FiniteField(7), 4, 4, algorithm='echelon_form', num_pivots=3); B # random
        [1 0 0 0]
        [0 1 0 6]
        [0 0 1 4]
        [0 0 0 0]
        sage: B.rank() == 3
        True
        sage: B.base_ring()
        Finite Field of size 7
        sage: B==B.rref()
        True

    TESTS:

    Rank of a matrix must be an integer. ::

        sage: random_matrix(QQ, 120, 56, algorithm='echelon_form', num_pivots=61/2)
        Traceback (most recent call last):
        ...
        TypeError: the number of pivots must be an integer.

    Matrices must be generated over exact fields. ::

        sage: random_matrix(RR, 40, 88, algorithm='echelon_form', num_pivots=39)
        Traceback (most recent call last):
        ...
        TypeError: the base ring must be exact.

    Matrices must have the number of pivot columns be less than or equal to the number of rows. ::

        sage: C=random_matrix(ZZ, 6,4, algorithm='echelon_form', num_pivots=7); C
        Traceback (most recent call last):
        ...
        ValueError: number of pivots cannot exceed the number of rows or columns.

    Matrices must have the number of pivot columns be less than or equal to the number of columns. ::

        sage: D=random_matrix(QQ, 1,3, algorithm='echelon_form', num_pivots=5); D
        Traceback (most recent call last):
        ...
        ValueError: number of pivots cannot exceed the number of rows or columns.

    Matrices must have the number of pivot columns be greater than zero. ::

        sage: random_matrix(QQ, 5, 4, algorithm='echelon_form', num_pivots=-1)
        Traceback (most recent call last):
        ...
        ValueError: the number of pivots must be zero or greater.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    import sage.gsl.probability_distribution as pd
    from sage.misc.prandom import randint

    try:
        num_pivots=ZZ(num_pivots)
    except TypeError:
        raise TypeError("the number of pivots must be an integer.")
    if num_pivots<0:
        raise ValueError("the number of pivots must be zero or greater.")
    ring = parent.base_ring()
    if not ring.is_exact():
        raise TypeError("the base ring must be exact.")
    num_row = parent.nrows()
    num_col = parent.ncols()
    if num_pivots>num_row or num_pivots>num_col:
        raise ValueError("number of pivots cannot exceed the number of rows or columns.")
    else:
        one=ring.one()
        # Create a matrix of the desired size to be modified and then returned.
        return_matrix=copy(parent.zero_matrix())
        pivots=[0] #Force first column to be a pivot. No harm if no pivots at all.
        # Probability distribution for the placement of leading one's.
        pivot_generator=pd.RealDistribution("beta",[1.6,4.3])
        while len(pivots)<num_pivots:
            pivot_column=int(pivot_generator.get_random_element()*num_col)
            if pivot_column not in pivots:
                pivots.append(pivot_column)
        pivots.sort()
        pivot_row=0
        # Use the list of pivot columns to set the pivot entries of the return_matrix to leading ones.
        while pivot_row<num_pivots:
            return_matrix[pivot_row,pivots[pivot_row]]=one
            pivot_row+=1
        if ring==QQ or ring==ZZ:
            # Keep track of the non-pivot columns by using the pivot_index, start at the first column to
            # the right of the initial pivot column, go until the first column to the left of the next
            # pivot column.
            for pivot_index in range(num_pivots-1):
                for non_pivot_column_index in range(pivots[pivot_index]+1,pivots[pivot_index+1]):
                    entry_generator1=pd.RealDistribution("beta",[6,4])
                    # Experimental distribution used to generate the values.
                    for non_pivot_column_entry in range(pivot_index+1):
                        sign1=(2*randint(0,1)-1)
                        return_matrix[non_pivot_column_entry,non_pivot_column_index]=sign1*int(entry_generator1.get_random_element()*((1-non_pivot_column_entry/return_matrix.ncols())*7))
            # Use index to fill entries of the columns to the right of the last pivot column.
            for rest_non_pivot_column in range(pivots[num_pivots-1]+1,num_col):
                entry_generator2=pd.RealDistribution("beta",[2.6,4])
                # experimental distribution to generate small values.
                for rest_entries in range(num_pivots):
                    sign2=(2*randint(0,1)-1)
                    return_matrix[rest_entries,rest_non_pivot_column]=sign2*int(entry_generator2.get_random_element()*5)
        else:
            for pivot_index in range(num_pivots-1):
                for non_pivot_column_index in range(pivots[pivot_index]+1,pivots[pivot_index+1]):
                    for non_pivot_column_entry in range(pivot_index+1):
                            return_matrix[non_pivot_column_entry,non_pivot_column_index]=ring.random_element()
            for rest_non_pivot_column in range(pivots[num_pivots-1]+1,num_col):
                for rest_entries in range(num_pivots):
                    return_matrix[rest_entries,rest_non_pivot_column]=ring.random_element()
    return return_matrix

def random_echelonizable_matrix(parent, rank, upper_bound=None):
    r"""
    Generate a matrix of a desired size and rank, over a desired ring, whose reduced
    row-echelon form has only integral values.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``rank`` - Rank of result, i.e the number of non-zero rows in the
      reduced row echelon form.

    - ``upper_bound`` - If designated, size control of the matrix entries is desired.
      Set ``upper_bound`` to 1 more than the maximum value entries can achieve.
      If None, no size control occurs. (default: None)

    OUTPUT:

    A matrix not in reduced row-echelon form with the desired dimensions and properties.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='echelonizable'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    Generated matrices have the desired dimensions, rank and entry size. The
    matrix in reduced row-echelon form has only integer entries. ::

        sage: from sage.matrix.constructor import random_echelonizable_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5, 6)
        sage: A=random_echelonizable_matrix(matrix_space, rank=4, upper_bound=40); A # random
        [  1  -1   1  -3  -4   6]
        [  5  -4   0   8   4  19]
        [ -3   3  -2   4   7 -16]
        [ -4   5  -7  26  31 -31]
        [  2  -3   4 -11 -14  17]
        sage: A.rank()
        4
        sage: max(map(abs,A.list()))<40
        True
        sage: A.rref()==A.rref().change_ring(ZZ)
        True

    An example with default settings (i.e. no entry size control). ::

        sage: C=random_matrix(QQ, 6, 7, algorithm='echelonizable', rank=5); C # random
        [  1   0   5  -2 -26 -16   0]
        [ -3   1 -19   6  97  61   1]
        [  0   4 -15  -1  71  50   3]
        [  2   4  -9   0  39  25   8]
        [  2   2   3  -3 -18  -9   3]
        [ -3  -4  -2  14  14  -6   4]
        sage: C.rank()
        5
        sage: C.rref()==C.rref().change_ring(ZZ)
        True

    A matrix without size control may have very large entry sizes. ::

        sage: D=random_matrix(ZZ, 7, 8, algorithm='echelonizable', rank=6); D # random
        [    9   -53  -255    45 -1519  4043  9819  3324]
        [    3   -14   -64     8  -369   972  2350   810]
        [    2   -14   -65     9  -377  1000  2420   829]
        [    4   -24  -116    21  -693  1846  4485  1516]
        [   -3    14    68   -16   426 -1134 -2767  -919]
        [   -5    21    92   -13   548 -1432 -3466 -1183]
        [    1    -9   -42     7  -254   670  1624   547]

    Matrices can be generated over any exact ring. ::

        sage: F.<a>=GF(2^3)
        sage: B=random_matrix(F, 4, 5, algorithm='echelonizable', rank=4, upper_bound=None); B # random
        [      a + 1 a^2 + a + 1         a^2           0           1]
        [          1           a       a + 1         a^2     a^2 + a]
        [          a         a^2 a^2 + a + 1       a + 1           1]
        [          1           0 a^2 + a + 1           0     a^2 + 1]
        sage: B.rank()
        4

    Square matrices over ZZ or QQ with full rank are always unimodular. ::

        sage: E=random_matrix(QQ, 7, 7, algorithm='echelonizable', rank=7); E # random
        [  1  -1   5  12 -24 -41  47]
        [  0   1  -1   3   0 -11  40]
        [  1  -1   6   6 -19 -20 -11]
        [ -2   1 -10 -12  35  44  -4]
        [  3  -1   9   7 -35 -40 -18]
        [  0   0   0  -4   4  13 -32]
        [  3  -3  11   6 -33 -31 -35]
        sage: det(E)
        1

    TESTS:

    Matrices must have a rank zero or greater zero. ::

        sage: random_matrix(QQ, 3, 4, algorithm='echelonizable', rank=-1)
        Traceback (most recent call last):
        ...
        ValueError: matrices must have rank zero or greater.

    The base ring must be exact. ::

        sage: random_matrix(RR, 3, 3, algorithm='echelonizable', rank=2)
        Traceback (most recent call last):
        ...
        TypeError: the base ring must be exact.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    from sage.misc.prandom import randint

    ring = parent.base_ring()
    rows = parent.nrows()
    if rank<0:
        raise ValueError("matrices must have rank zero or greater.")
    matrix = random_rref_matrix(parent, rank)
    # Entries of matrices over the ZZ or QQ can get large, entry size is regulated by finding the largest
    # entry of the resultant matrix after addition of scalar multiple of a row.
    if ring==QQ or ring==ZZ:
        # If upper_bound is not set, don't control entry size.
        if upper_bound==None:
            upper_bound=50
            size=False
        else:
            size=True
        if size:
            for pivots in range(len(matrix.pivots())-1,-1,-1):
            # keep track of the pivot column positions from the pivot column with the largest index to
            # the one with the smallest.
                row_index=0
                while row_index<rows:
                    # To each row in a pivot column add a scalar multiple of the pivot row.
                    # for full rank, square matrices, using only this row operation preserves the determinant of 1.
                    if pivots!=row_index:
                    # To ensure a leading one is not removed by the addition of the pivot row by its
                    # additive inverse.
                        matrix_copy=matrix.with_added_multiple_of_row(row_index,matrix.pivot_rows()[pivots],randint(-5,5))
                        # Range for scalar multiples determined experimentally.
                    if max(map(abs,matrix_copy.list()))<upper_bound:
                    # Continue if the the largest entry after a row operation is within the bound.
                        matrix=matrix_copy
                        row_index+=1
            # The leading one in row one has not been altered, so add a scalar multiple of a random row
            # to row one.
            row1=0
            if rows>1:
                while row1<1:
                    matrix_copy=matrix.with_added_multiple_of_row(0,randint(1,rows-1),randint(-3,3))
                    if max(map(abs,matrix_copy.list()))<upper_bound:
                        matrix=matrix_copy
                        row1+=1
        # If size control is not desired, the routine will run slightly faster, particularly with large matrices.
        else:
            for pivots in range(rank-1,-1,-1):
                row_index=0
                while row_index<rows:
                    if pivots==row_index:
                        row_index+=1
                    if pivots!=row_index and row_index!=rows:
                        matrix.add_multiple_of_row(row_index,matrix.pivot_rows()[pivots],randint(-5,5))
                        row_index+=1
            if rows>1:
                matrix.add_multiple_of_row(0,randint(1,rows-1),randint(-3,3))
    # If the matrix generated over a different ring, random elements from the designated ring are used as and
    # the routine is run similarly to the size unchecked version for rationals and integers.
    else:
        for pivots in range(rank-1,-1,-1):
            row_index=0
            while row_index<rows:
                if pivots==row_index:
                    row_index+=1
                if pivots!=row_index and row_index!=rows:
                    matrix.add_multiple_of_row(row_index,matrix.pivot_rows()[pivots],ring.random_element())
                    row_index+=1
        if rows>1:
            matrix.add_multiple_of_row(0,randint(1,rows-1),ring.random_element())
    return matrix

def random_subspaces_matrix(parent, rank=None):
    """
    Create a matrix of the designated size and rank whose right and
    left null spaces, column space, and row space have desirable
    properties that simplify the subspaces.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions, and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``rank`` - The desired rank of the return matrix (default: None).

    OUTPUT:

    A matrix whose natrual basis vectors for its four subspaces, when
    computed, have reasonably sized, integral valued, entries.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='subspaces'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A 6x8 matrix with designated rank of 3.  The four subspaces are
    determined using one simple routine in which we augment the
    original matrix with the equal row dimension identity matrix.  The
    resulting matrix is then put in reduced row-echelon form and the
    subspaces can then be determined by analyzing subdivisions of this
    matrix. See the four subspaces routine in [BEEZER]_ for more. ::

        sage: from sage.matrix.constructor import random_subspaces_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 6, 8)
        sage: B=random_subspaces_matrix(matrix_space, rank=3); B # random
        [ 113  339  -46  218 -243  641 -269 -306]
        [ -33  -99   13  -63   69 -185   77   90]
        [  35  105  -14   67  -74  197  -82  -95]
        [ -18  -54    7  -34   37 -100   41   49]
        [ -26  -78   10  -49   53 -144   59   71]
        [  -8  -24    3  -15   16  -44   18   22]
        sage: B.rank()
        3
        sage: B.nullity()
        3
        sage: (B.nrows(), B.ncols())
        (6, 8)
        sage: all([x in ZZ for x in B.list()])
        True
        sage: B_expanded=B.augment(identity_matrix(6)).rref()
        sage: all([x in ZZ for x in B_expanded.list()])
        True
        sage: B_expanded # random
        [  1   3   0   0   1   1   3  -2   0   0  -3   0  -9  16]
        [  0   0   1   0   3  -2  -1  -3   0   0   2   0  11 -27]
        [  0   0   0   1  -1   2  -3  -1   0   0   2   0   7 -14]
        [  0   0   0   0   0   0   0   0   1   0  -5   0  -3   2]
        [  0   0   0   0   0   0   0   0   0   1   1   0   1  -3]
        [  0   0   0   0   0   0   0   0   0   0   0   1  -1   1]
        sage: B_expanded.subdivide(B.nrows()-B.nullity(),B.ncols());B_expanded # random
        [  1   3   0   0   1   1   3  -2|  0   0  -3   0  -9  16]
        [  0   0   1   0   3  -2  -1  -3|  0   0   2   0  11 -27]
        [  0   0   0   1  -1   2  -3  -1|  0   0   2   0   7 -14]
        [-------------------------------+-----------------------]
        [  0   0   0   0   0   0   0   0|  1   0  -5   0  -3   2]
        [  0   0   0   0   0   0   0   0|  0   1   1   0   1  -3]
        [  0   0   0   0   0   0   0   0|  0   0   0   1  -1   1]
        sage: C=B_expanded.subdivision(0,0)
        sage: C # random
        [ 1  3  0  0  1  1  3 -2]
        [ 0  0  1  0  3 -2 -1 -3]
        [ 0  0  0  1 -1  2 -3 -1]
        sage: L=B_expanded.subdivision(1,1)
        sage: L # random
        [ 1  0 -5  0 -3  2]
        [ 0  1  1  0  1 -3]
        [ 0  0  0  1 -1  1]
        sage: B.right_kernel()==C.right_kernel()
        True
        sage: B.row_space()==C.row_space()
        True
        sage: B.column_space()==L.right_kernel()
        True
        sage: B.left_kernel()==L.row_space()
        True

    A matrix to show that the null space of the L matrix is the column space of the starting matrix. ::

        sage: A=random_matrix(QQ, 5, 7, algorithm='subspaces', rank=None); A # random
        [-31  12  -9 -27  21   2 -15]
        [105 -24   6 103 -30 -34  79]
        [ 29  -9   5  26 -14  -5  17]
        [233 -55  16 228 -71 -73 173]
        [-42  10  -3 -41  13  13 -31]
        sage: (A.nrows(), A.ncols())
        (5, 7)
        sage: all([x in ZZ for x in A.list()])
        True
        sage: A.nullity() # random
        1
        sage: A_expanded=A.augment(identity_matrix(5)).rref()
        sage: A_expanded # random
        [  1   0   0   0   0   1   0   0   3   7  25 151]
        [  0   1   0   0   1   2   1   0   5  21  84 493]
        [  0   0   1   0  -1   2   0   0   2  13  53 308]
        [  0   0   0   1   0  -1   1   0  -2  -3  -9 -57]
        [  0   0   0   0   0   0   0   1  -3   1   1  -2]
        sage: all([x in ZZ for x in A_expanded.list()])
        True
        sage: C=A_expanded.submatrix(0,0,A.nrows()-A.nullity(),A.ncols())
        sage: L=A_expanded.submatrix(A.nrows()-A.nullity(),A.ncols())
        sage: A.right_kernel()==C.right_kernel()
        True
        sage: A.row_space()==C.row_space()
        True
        sage: A.column_space()==L.right_kernel()
        True
        sage: A.left_kernel()==L.row_space()
        True

    TESTS:

    The designated rank of the L matrix cannot be greater than the number of desired rows. ::

        sage: random_matrix(QQ, 19, 20, algorithm='subspaces', rank=21)
        Traceback (most recent call last):
        ...
        ValueError: rank cannot exceed the number of rows or columns.

    REFERENCES:

        .. [BEEZER] `A First Course in Linear Algebra <http://linear.ups.edu/>`_.
           Robert A. Beezer, accessed 15 July 2010.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    import sage.gsl.probability_distribution as pd

    ring = parent.base_ring()
    rows = parent.nrows()
    columns = parent.ncols()

    if rank>rows or rank>columns:
        raise ValueError("rank cannot exceed the number of rows or columns.")
    # If rank is not designated, generate using probability distribution skewing to smaller numbers, always at least 1.
    if rank==None:
        left_nullity_generator=pd.RealDistribution("beta",[1.4,5.5])
        nullity=int(left_nullity_generator.get_random_element()*(rows-1)+1)
        rank=rows-nullity
    else:
        rank=rank
    nullity=rows-rank
    B=random_matrix(ring, rows, columns, algorithm='echelon_form', num_pivots=rank)
    # Create a nonsingular matrix whose columns will be used to stack a matrix over the L matrix, forming a nonsingular matrix.
    K_nonzero_columns=random_matrix(ring, rank, rank, algorithm='echelonizable', rank=rank)
    K=matrix(QQ,rank,rows)
    L=random_matrix(ring, nullity, rows, algorithm='echelon_form', num_pivots=nullity)
    for column in range(len(L.nonpivots())):
        for entry in range(rank):
            K[entry,L.nonpivots()[column]]=K_nonzero_columns[entry,column]
    J=K.stack(L)
    # by multiplying the B matrix by J.inverse() we hide the B matrix
    # of the solution using row operations required to change the solution
    # K matrix to the identity matrix.
    return J.inverse()*B

def random_unimodular_matrix(parent, upper_bound=None):
    """
    Generate a random unimodular (determinant 1) matrix of a desired size over a desired ring.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions
      and representation (dense/sparse) for the result.  The base ring
      must be exact.

    - ``upper_bound`` - For large matrices over QQ or ZZ,
      ``upper_bound`` is the largest value matrix entries can achieve.

    OUTPUT:

    An invertible matrix with the desired properties and determinant 1.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='unimodular'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A matrix size 5 over QQ. ::

        sage: from sage.matrix.constructor import random_unimodular_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5)
        sage: A=random_unimodular_matrix(matrix_space); A # random
        [  -8   31   85  148 -419]
        [   2   -9  -25  -45  127]
        [  -1   10   30   65 -176]
        [  -3   12   33   58 -164]
        [   5  -21  -59 -109  304]
        sage: det(A)
        1

    A matrix size 6 with entries no larger than 50. ::

        sage: B=random_matrix(ZZ, 7, algorithm='unimodular', upper_bound=50);B # random
        [ -1   0   3   1  -2   2   9]
        [  1   2  -5   0  14  19 -49]
        [ -3  -2  12   5  -6  -4  24]
        [  1   2  -9  -3   3   4  -7]
        [ -2  -1   7   2  -8  -5  31]
        [  2   2  -6  -3   8  16 -32]
        [  1   2  -9  -2   5   6 -12]
        sage: det(B)
        1

    A matrix over the number Field in `y` with defining polynomial `y^2-2y-2`. ::

        sage: y = var('y')
        sage: K=NumberField(y^2-2*y-2,'y')
        sage: C=random_matrix(K, 3, algorithm='unimodular');C # random
        [   2*y - 33 681*y - 787   31*y - 37]
        [      y + 6 -155*y + 83    -7*y + 4]
        [         -y   24*y + 51       y + 3]
        sage: det(C)
        1

    TESTS:

    Unimodular matrices are square. ::

        sage: random_matrix(QQ, 5, 6, algorithm='unimodular')
        Traceback (most recent call last):
        ...
        TypeError: a unimodular matrix must be square.

    Only matrices over ZZ and QQ can have size control. ::

        sage: F.<a>=GF(5^7)
        sage: random_matrix(F, 5, algorithm='unimodular', upper_bound=20)
        Traceback (most recent call last):
        ...
        TypeError: only matrices over ZZ or QQ can have size control.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    ring=parent.base_ring()
    size=parent.nrows()
    if parent.nrows()!=parent.ncols():
        raise TypeError("a unimodular matrix must be square.")
    if upper_bound!=None and (ring!=ZZ and ring!=QQ):
        raise TypeError("only matrices over ZZ or QQ can have size control.")
    if upper_bound==None:
        # random_echelonizable_matrix() always returns a determinant one matrix if given full rank.
        return random_matrix(ring, size, algorithm='echelonizable', rank=size)
    elif upper_bound!=None and (ring==ZZ or ring==QQ):
        return random_matrix(ring, size,algorithm='echelonizable',rank=size, upper_bound=upper_bound)


def random_diagonalizable_matrix(parent,eigenvalues=None,dimensions=None):
    """
    Create a random matrix that diagonalizes nicely.

    To be used as a teaching tool.  Return matrices have only real
    eigenvalues.

    INPUT:

    If eigenvalues and dimensions are not specified in a list,
    they will be assigned randomly.

    - ``parent`` - the desired size of the square matrix.

    - ``eigenvalues`` - the list of desired eigenvalues (default=None).

    - ``dimensions`` - the list of dimensions corresponding to each
      eigenspace (default=None).

    OUTPUT:

    A square, diagonalizable, matrix with only integer entries. The
    eigenspaces of this matrix, if computed by hand, give basis
    vectors with only integer entries.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='diagonalizable'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A diagonalizable matrix, size 5. ::

        sage: from sage.matrix.constructor import random_diagonalizable_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5)
        sage: A=random_diagonalizable_matrix(matrix_space); A # random
        [ 10  18   8   4 -18]
        [ 20  10   8   4 -16]
        [-60 -54 -22 -12  18]
        [-60 -54 -24  -6   6]
        [-20 -18  -8  -4   8]
        sage: A.eigenvalues() # random
        [10,6,2,-8,-10]
        sage: S=A.right_eigenmatrix()[1]; S # random
        [ 1  1  1  1  0]
        [ 1  1  1  0  1]
        [-3 -3 -4 -3 -3]
        [-3 -4 -3 -3 -3]
        [-1 -1 -1 -1 -1]
        sage: S_inverse=S.inverse(); S_inverse # random
        [ 1  1  1  1 -5]
        [ 0  0  0 -1  3]
        [ 0  0 -1  0  3]
        [ 0 -1  0  0 -1]
        [-1  0  0  0 -1]
        sage: S_inverse*A*S # random
        [ 10   0   0   0   0]
        [  0   6   0   0   0]
        [  0   0   2   0   0]
        [  0   0   0  -8   0]
        [  0   0   0   0 -10]

    A diagonalizable matrix with eigenvalues and dimensions designated,
    with a check that if eigenvectors were calculated by hand
    entries would all be integers. ::

        sage: B=random_matrix(QQ, 6, algorithm='diagonalizable', eigenvalues=[-12,4,6],dimensions=[2,3,1]); B # random
        [ -52   32  240 -464  -96 -520]
        [   6    4  -48   72   36   90]
        [  46  -32 -108  296  -12  274]
        [  24  -16  -64  164    0  152]
        [  18  -16    0   72  -48   30]
        [   2    0  -16   24   12   34]
        sage: all([x in ZZ for x in (B-(-12*identity_matrix(6))).rref().list()])
        True
        sage: all([x in ZZ for x in (B-(4*identity_matrix(6))).rref().list()])
        True
        sage: all([x in ZZ for x in (B-(6*identity_matrix(6))).rref().list()])
        True
        sage: S=B.right_eigenmatrix()[1]; S_inverse=S.inverse(); S_inverse*B*S # random
        [  6   0   0   0   0   0]
        [  0 -12   0   0   0   0]
        [  0   0 -12   0   0   0]
        [  0   0   0   4   0   0]
        [  0   0   0   0   4   0]
        [  0   0   0   0   0   4]

    TESTS:

    Eigenvalues must all be integers. ::

        sage: random_matrix(QQ,3,algorithm='diagonalizable', eigenvalues=[2+I,2-I,2],dimensions=[1,1,1])
        Traceback (most recent call last):
        ...
        TypeError: eigenvalues must be integers.

    Diagonal matrices must be square. ::

        sage: random_matrix(QQ, 5, 7, algorithm='diagonalizable', eigenvalues=[-5,2,-3], dimensions=[1,1,3])
        Traceback (most recent call last):
        ...
        TypeError: a diagonalizable matrix must be square.

    A list of eigenvalues must be accompanied with a list of dimensions. ::

        sage: random_matrix(QQ,10,algorithm='diagonalizable',eigenvalues=[4,8])
        Traceback (most recent call last):
        ...
        ValueError: the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.

    A list of dimensions must be accompanied with a list of eigenvalues. ::

        sage: random_matrix(QQ, 10,algorithm='diagonalizable',dimensions=[2,2,4,2])
        Traceback (most recent call last):
        ...
        ValueError: the list of dimensions must have a list of corresponding eigenvalues.

    The sum of the eigenvalue dimensions must equal the size of the matrix. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6,-1],dimensions=[2,3,5,1])
        Traceback (most recent call last):
        ...
        ValueError: the size of the matrix must equal the sum of the dimensions.

    Each eigenspace dimension must be at least 1. ::

        sage: random_matrix(QQ,9,algorithm='diagonalizable',eigenvalues=[-15,22,8,-4,90,12],dimensions=[4,2,2,4,-3,0])
        Traceback (most recent call last):
        ...
        ValueError: eigenspaces must have a dimension of at least 1.

    Each eigenvalue must have a corresponding eigenspace dimension. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6,-1],dimensions=[4,3,5])
        Traceback (most recent call last):
        ...
        ValueError: each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.

    Each dimension must have an eigenvalue paired to it. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6],dimensions=[2,3,5,2])
        Traceback (most recent call last):
        ...
        ValueError: each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.

    TODO:

    Modify the routine to allow for complex eigenvalues.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    from sage.misc.prandom import randint

    size=parent.nrows()
    if parent.nrows()!=parent.ncols():
        raise TypeError("a diagonalizable matrix must be square.")
    if eigenvalues!=None and dimensions==None:
        raise ValueError("the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.")
    if eigenvalues==None and dimensions!=None:
        raise ValueError("the list of dimensions must have a list of corresponding eigenvalues.")
    if eigenvalues==None and dimensions==None:
        values=[]
        #create a list with "size" number of entries
        for eigen_index in range(size):
            eigenvalue=randint(-10,10)
            values.append(eigenvalue)
        values.sort()
        dimensions=[]
        eigenvalues=[]
        #create a list with no duplicate values to be the eigenvalues
        for eigenvalue in range(size):
            if values[eigenvalue] not in eigenvalues:
                eigenvalues.append(values[eigenvalue])
        for dimension in range(len(eigenvalues)):
            #dimension is equal to how many times an eigenvalue was generated in the 'values' list
            dimensions.append(values.count(eigenvalues[dimension]))
    size_check=0
    for check in range(len(dimensions)):
        size_check=size_check+dimensions[check]
    if not map(lambda x: x in ZZ, eigenvalues)==[True]*len(eigenvalues):
        raise TypeError("eigenvalues must be integers.")
    if size!=size_check:
        raise ValueError("the size of the matrix must equal the sum of the dimensions.")
    if min(dimensions)<1:
        raise ValueError("eigenspaces must have a dimension of at least 1.")
    if len(eigenvalues)!=len(dimensions):
        raise ValueError("each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.")
    #sort the dimensions in order of increasing size, and sort the eigenvalues list in an identical fashion, to maintain corresponding values.
    dimensions_sort=zip(dimensions,eigenvalues)
    dimensions_sort.sort()
    dimensions=[x[0] for x in dimensions_sort]
    eigenvalues=[x[1] for x in dimensions_sort]
    #Create the matrix of eigenvalues on the diagonal.  Use a lower limit and upper limit determined by the eigenvalue dimensions.
    diagonal_matrix=matrix(QQ,size)
    up_bound=0
    low_bound=0
    for row_index in range(len(dimensions)):
        up_bound=up_bound+dimensions[row_index]
        for entry in range(low_bound,up_bound):
            diagonal_matrix[entry,entry]=eigenvalues[row_index]
        low_bound=low_bound+dimensions[row_index]
    # Create a matrix to hold each of the eigenvectors as its columns, begin with an identity matrix so that after row and column
    # operations the resulting matrix will be unimodular.
    eigenvector_matrix=matrix(QQ,size,size,1)
    upper_limit=0
    lower_limit=0
    #run the routine over the necessary number of columns corresponding eigenvalue dimension.
    for dimension_index in range(len(dimensions)-1):
        upper_limit=upper_limit+dimensions[dimension_index]
        lowest_index_row_with_one=size-dimensions[dimension_index]
        #assign a one to the row that is the eigenvalue dimension rows up from the bottom row then assign ones diagonally down to the right.
        for eigen_ones in range(lower_limit,upper_limit):
            eigenvector_matrix[lowest_index_row_with_one,eigen_ones]=1
            lowest_index_row_with_one+=1
        lower_limit=lower_limit+dimensions[dimension_index]
    #Create a list to give the eigenvalue dimension corresponding to each column.
    dimension_check=[]
    for i in range(len(dimensions)):
        for k in range(dimensions[i]):
            dimension_check.append(dimensions[i])
    #run routine over the rows that are in the range of the protected ones.  Use addition of column multiples to fill entries.
    for dimension_multiplicity in range(max(dimensions),min(dimensions),-1):
        highest_one_row=size-dimension_multiplicity
        highest_one_column=0
        #find the column with the protected one in the lowest indexed row.
        while eigenvector_matrix[highest_one_row,highest_one_column]==0:
            highest_one_column+=1
        #dimension_check determines if column has a low enough eigenvalue dimension to take a column multiple.
        for bottom_entry_filler in range(len(dimension_check)):
            if dimension_check[bottom_entry_filler]<dimension_multiplicity and eigenvector_matrix[highest_one_row,bottom_entry_filler]==0:
                # randint range determined experimentally to keep entries manageable.
                eigenvector_matrix.add_multiple_of_column(bottom_entry_filler,highest_one_column,randint(-4,4))
    #Fill remaining rows using scalar row addition.
    for row in range(size-max(dimensions),size):
        for upper_row in range(size-max(dimensions)):
            # range of multiplier determined experimentally so that entries stay manageable for small matrices
            eigenvector_matrix.add_multiple_of_row(upper_row,row,randint(-4,4))
    return eigenvector_matrix*diagonal_matrix*(eigenvector_matrix.inverse())
